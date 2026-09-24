"""peptidehlacpp command-line tool.

Subcommands
-----------
predict     Peptide-HLA class-I binding: per-allele PSSM log10(IC50) +
            allele-conditioned CNN binder probability + the paper's
            z-scored ensemble score (batch mode).
screen      Score sequences for cell-penetrating-peptide character with the
            trained CPP CNN classifier plus biophysical descriptors.
design-cpp  Sample novel CPP candidates from the trained GRU generator and
            run the in-silico screening cascade (classifier + charge +
            hydrophobic moment + hydropathy + novelty vs training set).

All models load from committed artifacts (results/). Use --weights-dir to
point at another checkout. Run from the repository root or install with
`pip install .` and pass --weights-dir.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

from .design.generator import CPPGenerator, sample, screen_candidates
from .features import AA_INDEX, MAX_LEN_PHLA, hydrophobic_moment, mean_hydropathy, net_charge
from .models.cnn import PHLACNN
from .training.train_cpp import CPPCNN, enc
from .training.train_phla import enc_batch

CANONICAL = set(AA_INDEX)


def default_weights_dir() -> Path:
    here = Path(__file__).resolve()
    for up in here.parents:
        cand = up / "results"
        if (cand / "phla_cnn.pt").exists():
            return cand
    return Path("results")


def read_fasta(path: str) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0] or f"seq{len(out)}", []
        elif line:
            buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def check_seq(s: str) -> str:
    s = s.upper()
    if not s or not set(s) <= CANONICAL:
        raise ValueError(f"non-canonical or empty sequence: {s!r}")
    return s


class PHLAPredictor:
    """PSSM + CNN + z-scored ensemble, reproducing the paper pipeline."""

    def __init__(self, wdir: Path):
        d = np.load(wdir / "pssm_weights.npz", allow_pickle=False)
        self.pssm_alleles = list(d["alleles"])
        self.W, self.b = d["W"], d["b"]
        ck = torch.load(wdir / "phla_cnn.pt", map_location="cpu", weights_only=False)
        self.cnn_alleles = list(ck["alleles"])
        self.amap = {a: i for i, a in enumerate(self.cnn_alleles)}
        self.cnn = PHLACNN(n_alleles=len(self.cnn_alleles))
        self.cnn.load_state_dict(ck["state"])
        self.cnn.eval()

    @property
    def alleles(self) -> list[str]:
        return sorted(set(self.pssm_alleles) & set(self.cnn_alleles))

    def pssm_log10ic50(self, allele: str, seqs: list[str]) -> np.ndarray:
        i = self.pssm_alleles.index(allele)
        X = np.zeros((len(seqs), MAX_LEN_PHLA * 20))
        for r, s in enumerate(seqs):
            for p, a in enumerate(s[:MAX_LEN_PHLA]):
                X[r, p * 20 + AA_INDEX[a]] = 1.0
        return X @ self.W[i].ravel() + self.b[i]

    def cnn_p_binder(self, allele: str, seqs: list[str]) -> np.ndarray:
        X, m, al = enc_batch(seqs, np.array([allele] * len(seqs)), self.amap, False)
        with torch.no_grad():
            logits = self.cnn(X, m, al)[1].numpy()
        return 1.0 / (1.0 + np.exp(-logits))

    def predict(self, allele: str, seqs: list[str]) -> list[dict]:
        if allele not in self.alleles:
            raise ValueError(f"allele {allele!r} not covered; run `peptidehlacpp predict --list-alleles`")
        log_ic50 = self.pssm_log10ic50(allele, seqs)
        p_cnn = self.cnn_p_binder(allele, seqs)
        pssm_score = -log_ic50  # higher = stronger binder
        if len(seqs) > 1:
            z1 = (pssm_score - pssm_score.mean()) / (pssm_score.std() + 1e-9)
            z2 = (p_cnn - p_cnn.mean()) / (p_cnn.std() + 1e-9)
            ens = z1 + z2
        else:
            ens = np.full(len(seqs), np.nan)
        return [{
            "sequence": s,
            "allele": allele,
            "pssm_log10_ic50": round(float(log_ic50[i]), 4),
            "pssm_ic50_nM": round(float(10 ** log_ic50[i]), 1),
            "cnn_p_binder": round(float(p_cnn[i]), 4),
            "ensemble_z": None if np.isnan(ens[i]) else round(float(ens[i]), 4),
            "binder_500nM": bool(p_cnn[i] >= 0.5),
        } for i, s in enumerate(seqs)]


class CPPScreener:
    def __init__(self, wdir: Path):
        self.clf = CPPCNN()
        self.clf.load_state_dict(torch.load(wdir / "cpp_cnn.pt", map_location="cpu", weights_only=False))
        self.clf.eval()

    def p_cpp(self, seqs: list[str]) -> np.ndarray:
        X, m = enc(seqs)
        with torch.no_grad():
            logits = self.clf(X, m).numpy()
        return 1.0 / (1.0 + np.exp(-logits))

    def screen(self, seqs: list[str]) -> list[dict]:
        p = self.p_cpp(seqs)
        rows = [{
            "sequence": s,
            "p_cpp": round(float(p[i]), 4),
            "net_charge": round(net_charge(s), 2),
            "hydrophobic_moment": round(hydrophobic_moment(s), 3),
            "mean_hydropathy": round(mean_hydropathy(s), 3),
            "predicted_cpp": bool(p[i] >= 0.5),
        } for i, s in enumerate(seqs)]
        rows.sort(key=lambda r: -r["p_cpp"])
        return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="peptidehlacpp", description=__doc__.splitlines()[0])
    ap.add_argument("--weights-dir", default=None,
                    help="directory containing model artifacts (default: auto-detect results/)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_pred = sub.add_parser("predict", help="peptide-HLA binding prediction")
    p_pred.add_argument("--allele", default="HLA-A*02:01")
    p_pred.add_argument("--peptides", nargs="*", default=[])
    p_pred.add_argument("--fasta", default=None)
    p_pred.add_argument("--list-alleles", action="store_true")
    p_pred.add_argument("--json", action="store_true")

    p_scr = sub.add_parser("screen", help="CPP classification screen")
    p_scr.add_argument("--peptides", nargs="*", default=[])
    p_scr.add_argument("--fasta", default=None)
    p_scr.add_argument("--json", action="store_true")

    p_des = sub.add_parser("design-cpp", help="generate + screen novel CPP candidates")
    p_des.add_argument("--n", type=int, default=10, help="candidates to return")
    p_des.add_argument("--sample", type=int, default=1000, help="raw samples to draw")
    p_des.add_argument("--temperature", type=float, default=0.9)
    p_des.add_argument("--seed", type=int, default=17)
    p_des.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)
    wdir = Path(args.weights_dir) if args.weights_dir else default_weights_dir()
    if not wdir.exists():
        ap.error(f"weights dir {wdir} not found")

    if args.cmd == "predict":
        pred = PHLAPredictor(wdir)
        if args.list_alleles:
            print("\n".join(pred.alleles))
            return 0
        seqs = [check_seq(s) for s in args.peptides]
        if args.fasta:
            seqs += [check_seq(s) for _, s in read_fasta(args.fasta)]
        if not seqs:
            ap.error("predict needs --peptides and/or --fasta")
        rows = pred.predict(args.allele, seqs)
        if len(seqs) == 1:
            print("note: ensemble_z needs a batch of >=2 peptides (z-scored within the batch); "
                  "single-peptide runs report pssm/cnn only", file=sys.stderr)
    elif args.cmd == "screen":
        seqs = [check_seq(s) for s in args.peptides]
        if args.fasta:
            seqs += [check_seq(s) for _, s in read_fasta(args.fasta)]
        if not seqs:
            ap.error("screen needs --peptides and/or --fasta")
        rows = CPPScreener(wdir).screen(seqs)
    else:  # design-cpp
        gen = CPPGenerator()
        gen.load_state_dict(torch.load(wdir / "cpp_generator.pt", map_location="cpu", weights_only=False))
        scr = CPPScreener(wdir)
        train_seqs = set((wdir / "cpp_train_sequences.txt").read_text().split())
        cand = sample(gen, args.sample, temperature=args.temperature, seed=args.seed)
        passed = screen_candidates(sorted(set(cand)), lambda s: float(scr.p_cpp([s])[0]), train_seqs)
        rows = passed[: args.n]
        if not rows:
            print("no candidates passed the cascade; raise --sample or --temperature", file=sys.stderr)
            return 1

    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        if not rows:
            print("(no rows)")
            return 0
        keys = list(rows[0])
        print("\t".join(keys))
        for r in rows:
            print("\t".join(str(r[k]) for k in keys))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
