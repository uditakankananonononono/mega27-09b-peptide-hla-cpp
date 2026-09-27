"""B3 eluted-ligand benchmark (PREREG_B3_EL_BENCHMARK_2026-09-27.md, locked
before scoring). Phases (resume-safe, results/b3_partial/):
  parse     stream the 5 c00?_el files -> per-allele {peptide: target} json shards
  score     score covered alleles (full ensemble) + uncovered (B1 arm1 mode)
  finalize  falsifier + overlap audit + bootstrap + results/b3_el_benchmark.json
"""
import argparse, glob, json, sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, "src")
import numpy as np

from peptidehlacpp.eval import metrics as M

PART = Path("results/b3_partial")
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
SEED_BOOT, SEED_SHUFFLE = 23, 29


def norm_allele(raw):
    # HLA-C14:02 -> HLA-C*14:02
    if not raw.startswith("HLA-"):
        return None
    body = raw[4:]
    if "*" in body or "-" in body:
        return None  # unexpected format, treat as cell line
    # split locus letter(s) from digits: e.g. C14:02, A02:01, B44:02
    i = 0
    while i < len(body) and body[i].isalpha():
        i += 1
    if i == 0 or i == len(body):
        return None
    return f"HLA-{body[:i]}*{body[i:]}"


def phase_parse():
    out_p = PART / "el_monoallelic.json"
    if out_p.exists():
        print("parse done, skip"); return
    per_allele = defaultdict(dict)
    n_total = n_ma_excluded = n_bad = 0
    for f in sorted(glob.glob("data/external/NetMHCpan_train/c00?_el")):
        for line in open(f):
            n_total += 1
            parts = line.split()
            if len(parts) < 3:
                n_bad += 1
                continue
            pep, tgt, raw = parts[0], parts[1], parts[2]
            a = norm_allele(raw)
            if a is None:
                n_ma_excluded += 1
                continue
            if not (8 <= len(pep) <= 14) or not set(pep) <= AA20:
                n_bad += 1
                continue
            t = 1 if tgt == "1" else 0
            # ligand overrides decoy on (peptide, allele) conflict
            per_allele[a][pep] = max(per_allele[a].get(pep, 0), t)
        print(f"{f} parsed; total {n_total}", flush=True)
    out = {a: d for a, d in per_allele.items()}
    meta = {"n_lines": n_total, "n_ma_excluded": n_ma_excluded, "n_bad": n_bad,
            "n_alleles": len(out), "n_rows": sum(len(d) for d in out.values())}
    json.dump(out, open(out_p, "w"))
    json.dump(meta, open(PART / "parse_meta.json", "w"), indent=1)
    print(json.dumps(meta, indent=1), flush=True)


def phase_score():
    import torch
    from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv
    from peptidehlacpp.models.pssm import AllelePSSM
    from peptidehlacpp.models.cnn import PHLACNN
    from peptidehlacpp.features import stacked_enc, MAX_LEN_PHLA
    train_alleles = sorted({r["allele"] for r in json.load(open("results/per_allele_analysis.json"))["per_allele"]})
    tset = set(train_alleles)
    el = json.load(open(PART / "el_monoallelic.json"))
    covered = {a: d for a, d in el.items() if a in tset}
    uncovered = {a: d for a, d in el.items() if a not in tset}
    print(f"covered {len(covered)} alleles / {sum(len(d) for d in covered.values())} rows; "
          f"uncovered {len(uncovered)} / {sum(len(d) for d in uncovered.values())}", flush=True)
    # fit per-allele PSSMs on IEDB BA rows (production recipe; no EL row touches any fit)
    df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
    examples = aggregate(df)
    by_al = defaultdict(list)
    for e in examples:
        if e.allele in tset:
            by_al[e.allele].append(e)
    pssms = {a: AllelePSSM(ridge=1.0).fit([e.sequence for e in exs],
             np.array([e.log_ic50 for e in exs])) for a, exs in by_al.items()}
    amap = {a: i for i, a in enumerate(train_alleles)}
    cnn = PHLACNN(n_alleles=len(train_alleles))
    cnn.load_state_dict(torch.load("results/b1_partial/b1_cnn_best.pt")["state"])
    cnn.eval()
    with torch.no_grad():
        mean_vec = cnn.allele_emb.weight.mean(dim=0, keepdim=True)

    def cnn_scores(peps, allele=None):
        out = []
        with torch.no_grad():
            for s in range(0, len(peps), 2048):
                chunk = peps[s:s + 2048]
                X = np.stack([stacked_enc(p, MAX_LEN_PHLA) for p in chunk])
                mask = np.zeros((len(chunk), MAX_LEN_PHLA), dtype=np.float32)
                for r, p in enumerate(chunk):
                    mask[r, : len(p)] = 1.0
                Xt, mt = torch.from_numpy(X), torch.from_numpy(mask)
                h = Xt.transpose(1, 2)
                for b in cnn.blocks:
                    h = b(h)
                m1 = mt.unsqueeze(1)
                h = h * m1
                mean_pool = h.sum(dim=2) / m1.sum(dim=2).clamp(min=1.0)
                max_pool = h.masked_fill(m1 == 0, -1e4).amax(dim=2)
                if allele is not None:
                    emb = cnn.allele_emb.weight[amap[allele]].expand(len(chunk), -1)
                else:
                    emb = mean_vec.expand(len(chunk), -1)
                z = cnn.head(torch.cat([mean_pool, max_pool, emb], dim=1))
                out.append(cnn.cls_out(z).squeeze(-1).numpy())
        return np.concatenate(out)

    def z(x):
        x = np.asarray(x, dtype=np.float64)
        return (x - x.mean()) / (x.std() + 1e-9)

    res = {"covered": {}, "uncovered": {}}
    for group, scorer in (("covered", covered), ("uncovered", uncovered)):
        ck = PART / f"scores_{group}.json"
        done = json.load(open(ck)) if ck.exists() else {}
        for a, d in sorted(scorer.items()):
            if a in done:
                continue
            peps = sorted(d)
            labels = np.array([d[p] for p in peps])
            if group == "covered":
                s_pssm = -pssms[a].predict(peps)
                s = z(s_pssm) + z(cnn_scores(peps, allele=a))
            else:
                s = cnn_scores(peps, allele=None)
            done[a] = {"n": len(peps), "n_pos": int(labels.sum()),
                       "auroc": M.auc(labels, s) if 0 < labels.sum() < len(labels) else None,
                       "auc0.1": M.auc_top_frac(labels, s, frac=0.1) if 0 < labels.sum() < len(labels) else None,
                       "scores": [round(float(v), 5) for v in s],
                       "labels": labels.tolist()}
            json.dump(done, open(ck, "w"))
            print(f"{group} {a} auroc={done[a]['auroc']}", flush=True)
        res[group] = done
    print("score done", flush=True)


def phase_finalize():
    import glob as g
    cov = json.load(open(PART / "scores_covered.json"))
    unc = json.load(open(PART / "scores_uncovered.json"))
    meta = json.load(open(PART / "parse_meta.json"))
    # falsifier: seed-29 label shuffle within largest covered allele
    big = max(cov, key=lambda a: cov[a]["n"])
    f = cov[big]
    rng = np.random.default_rng(SEED_SHUFFLE)
    fals = M.auc(rng.permutation(np.array(f["labels"])), np.array(f["scores"]))
    # overlap audit: EL ligands also present as IEDB BA training peptides
    from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv
    examples = aggregate(load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv"))
    ba_peps = {e.sequence for e in examples}
    audit = {}
    for group, res in (("covered", cov), ("uncovered", unc)):
        for a, v in res.items():
            peps = sorted({p for p in v.get("scores", [])})  # placeholder, see below
    # per-allele overlap needs peptide lists; recompute from el_monoallelic
    el = json.load(open(PART / "el_monoallelic.json"))
    for a, d in el.items():
        if not d:
            continue
        hits = sum(1 for p in d if d[p] == 1 and p in ba_peps)
        np_ = sum(1 for p in d if d[p] == 1)
        audit[a] = {"ligand_overlap_frac": round(hits / np_, 4) if np_ else None, "n_ligands": np_}
    rngb = np.random.default_rng(SEED_BOOT)
    def summ(res):
        aus = np.array([v["auroc"] for v in res.values() if v["auroc"] is not None])
        a01 = np.array([v["auc0.1"] for v in res.values() if v["auc0.1"] is not None])
        boots = [float(rngb.choice(aus, size=len(aus), replace=True).mean()) for _ in range(10000)]
        return {"n_alleles": len(aus),
                "mean_auroc": float(aus.mean()),
                "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                "mean_auc0.1": float(a01.mean()),
                "per_allele": {a: {k: v[k] for k in ("n", "n_pos", "auroc", "auc0.1")}
                               for a, v in res.items()}}
    sc, su = summ(cov), summ(unc)
    def decision(m, m01):
        if m >= 0.90 and m01 >= 0.50:
            return "SUCCESS"
        if m >= 0.75:
            return "PARTIAL"
        return "FAILURE"
    out = {
        "prereg": "docs/PREREG_B3_EL_BENCHMARK_2026-09-27.md",
        "data": {"source": "NetMHCpan-4.1 training release (EL partitions, mono-allelic rows)",
                 "sha256": "06f2c9f20bb959238bf5d601fca0489a0ed3f17648b952f30640205afca8f9b4",
                 "parse_meta": meta},
        "falsifier": {"allele": big, "seed": SEED_SHUFFLE, "shuffled_auroc": fals,
                      "in_0.40_0.60": bool(0.40 <= fals <= 0.60)},
        "ligand_overlap_audit_vs_iedb_ba": audit,
        "covered_arm": sc, "uncovered_arm": su,
        "decision": {"covered": decision(sc["mean_auroc"], sc["mean_auc0.1"]),
                     "uncovered": decision(su["mean_auroc"], su["mean_auc0.1"])},
    }
    json.dump(out, open("results/b3_el_benchmark.json", "w"), indent=1)
    print(json.dumps({"falsifier": out["falsifier"], "decision": out["decision"],
                      "covered_mean_auroc": sc["mean_auroc"], "covered_mean_auc0.1": sc["mean_auc0.1"],
                      "uncovered_mean_auroc": su["mean_auroc"]}, indent=1))


def main():
    PART.mkdir(exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True, choices=["parse", "score", "finalize"])
    a = ap.parse_args()
    {"parse": phase_parse, "score": phase_score, "finalize": phase_finalize}[a.phase]()


if __name__ == "__main__":
    main()
