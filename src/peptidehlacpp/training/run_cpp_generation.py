"""Train the CPP generator on redundancy-filtered CPPsite positives, sample
candidates, and screen them through the trained CNN classifier + biophysical
filters. Writes results/cpp_designs.json."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from ..data.cppsite import (dedup_exact, natural_only, parse_fasta,
                            redundancy_filter)
from ..design.generator import sample, screen_candidates, train_generator
from ..training.train_cpp import CPPCNN, enc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pos", default="data/raw/cppsite2_natural.fa")
    ap.add_argument("--classifier", default="results/cpp_cnn.pt")
    ap.add_argument("--out", default="results/cpp_designs.json")
    ap.add_argument("--n-sample", type=int, default=4000)
    ap.add_argument("--gen-epochs", type=int, default=40)
    args = ap.parse_args()

    pos = redundancy_filter(dedup_exact(natural_only(parse_fasta(args.pos))))
    pos = [p for p in pos if 8 <= len(p.sequence) <= 40]
    seqs = [p.sequence for p in pos]
    print(f"training generator on {len(seqs)} CPPs", flush=True)
    gen = train_generator(seqs, epochs=args.gen_epochs)
    torch.save(gen.state_dict(), "results/cpp_generator.pt")

    clf = CPPCNN()
    clf.load_state_dict(torch.load(args.classifier))
    clf.eval()

    def clf_score(s: str) -> float:
        X, m = enc([s])
        with torch.no_grad():
            return float(torch.sigmoid(clf(X, m))[0])

    cand = sample(gen, args.n_sample, temperature=0.9)
    uniq = sorted(set(cand))
    print(f"sampled {len(cand)} ({len(uniq)} unique)", flush=True)
    passed = screen_candidates(uniq, clf_score, set(seqs))
    print(f"passed cascade: {len(passed)}", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps({
        "n_training_cpp": len(seqs), "n_sampled": len(cand),
        "n_unique": len(uniq), "n_passed": len(passed),
        "top_candidates": passed[:50]}, indent=2))
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
