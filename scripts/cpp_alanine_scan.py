"""Computational falsification test for the basic-run design rule.

Prediction (paper, CPP 3-mer weights section): mutating one basic 3-run
(Arg/Lys) to Ala in a 9B-CPP named candidate should reduce its 3-mer LR
log-odds by ~the top-weight magnitude and drop p below the 0.7 cascade
threshold. Here we actually run that scan with the refit LR (same split
protocol as the benchmark) on the 18 named candidates and record outcomes
either way. Saves results/cpp_alanine_scan.json.
"""
from __future__ import annotations

import json
import sys
sys.path.insert(0, "src")

import numpy as np
from sklearn.linear_model import LogisticRegression

from peptidehlacpp.data.cppsite import (dedup_exact, natural_only,
                                        parse_fasta, redundancy_filter,
                                        sample_length_matched_windows,
                                        kmer_jaccard)
from peptidehlacpp.eval import metrics as M
from peptidehlacpp.training.train_cpp import build_kmer_features

MAX_LEN = 40
pos = redundancy_filter(dedup_exact(natural_only(parse_fasta("data/raw/cppsite2_natural.fa"))))
pos = [p for p in pos if 8 <= len(p.sequence) <= MAX_LEN]
neg_pool = []
for np_ in ["data/raw/uniprot_reviewed_len8_35.fasta",
            "data/raw/uniprot_human_reviewed_len60_400.fasta"]:
    neg_pool.extend(p for p in natural_only(parse_fasta(np_)) if 8 <= len(p.sequence))
neg_pool = dedup_exact(neg_pool)
neg = sample_length_matched_windows(neg_pool, pos, n_per_pos=1)
rng = np.random.default_rng(11)
order = rng.permutation(len(pos))
n_test = int(0.15 * len(pos)); n_val = int(0.1 * len(pos))
idx_te, idx_va = order[:n_test], order[n_test:n_test + n_val]
idx_tr = order[n_test + n_val:]
te_pos = [pos[i].sequence for i in idx_te]
idx_tr = [i for i in idx_tr
          if all(kmer_jaccard(pos[i].sequence, t) < 0.6 for t in te_pos)]
def mk(idxs_p):
    n_same = int(round(len(idxs_p) * len(neg) / len(pos)))
    neg_sel = [neg[j] for j in rng.choice(len(neg), size=n_same, replace=False)]
    seqs = [pos[i].sequence for i in idxs_p] + [n.sequence for n in neg_sel]
    labels = np.array([1] * len(idxs_p) + [0] * len(neg_sel), dtype=np.float32)
    return seqs, labels
tr_s, tr_y = mk(idx_tr)
Xtr, vocab = build_kmer_features(tr_s)
lr = LogisticRegression(max_iter=2000, C=1.0, n_jobs=2).fit(Xtr, tr_y)

def logodds(seq):
    X, _ = build_kmer_features([seq], vocab=vocab)
    return float(X @ lr.coef_[0] + lr.intercept_[0])

def p(seq):
    z = logodds(seq)
    return 1.0 / (1.0 + np.exp(-z))

named = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
rows = []
for c in named:
    s = c["sequence"]
    z0 = logodds(s)
    # all maximal basic runs of length >= 3
    runs = []
    i = 0
    while i < len(s):
        if s[i] in "RK":
            j = i
            while j < len(s) and s[j] in "RK":
                j += 1
            if j - i >= 3:
                runs.append((i, j))
            i = j
        else:
            i += 1
    best = None
    for (i, j) in runs:
        mut = s[:i] + "A" * (j - i) + s[j:]
        dz = z0 - logodds(mut)
        if best is None or dz > best["drop"]:
            best = {"run": s[i:j], "span": [i, j], "mutant": mut, "drop": dz,
                    "p_mutant": p(mut)}
    rows.append({"name": c["name"], "sequence": s, "z0": z0, "p0": p(s),
                 "n_basic_runs_ge3": len(runs), "worst_run": best})

# second variant: strip ALL basic residues (R/K -> A everywhere)
for r in rows:
    s = r["sequence"]
    mut = "".join("A" if a in "RK" else a for a in s)
    r["allRK_mutant"] = mut
    r["p_allRK"] = p(mut)
    r["drop_allRK"] = r["z0"] - logodds(mut)
drops = np.array([r["worst_run"]["drop"] for r in rows if r["worst_run"]])
drops_all = np.array([r["drop_allRK"] for r in rows])
frac_below = float(np.mean([r["worst_run"]["p_mutant"] < 0.7
                            for r in rows if r["worst_run"]]))
summary = {
    "n_candidates": len(rows),
    "n_with_basic_run_ge3": int((drops.size)),
    "mean_logodds_drop": float(drops.mean()),
    "min_drop": float(drops.min()), "max_drop": float(drops.max()),
    "frac_mutants_below_cascade_0.7": frac_below,
    "frac_mutants_below_0.5": float(np.mean([r["worst_run"]["p_mutant"] < 0.5
                                             for r in rows if r["worst_run"]])),
    "mean_drop_allRK": float(drops_all.mean()),
    "frac_allRK_below_0.5": float(np.mean([r["p_allRK"] < 0.5 for r in rows])),
    "frac_allRK_below_0.1": float(np.mean([r["p_allRK"] < 0.1 for r in rows])),
    "candidates": rows,
}
json.dump(summary, open("results/cpp_alanine_scan.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "candidates"}, indent=1))
