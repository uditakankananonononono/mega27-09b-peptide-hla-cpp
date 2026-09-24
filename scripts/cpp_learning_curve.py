"""CPP 3-mer LR learning curve: is the classifier data-limited at n~600?

Same split protocol as the benchmark (seed 11 prep); refit LR on train
fractions [0.1..1.0] x 4 seeds, evaluate on the fixed homology-guarded
test set. Saves results/cpp_learning_curve.json.
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
def mk(idxs_p, rng_):
    n_same = int(round(len(idxs_p) * len(neg) / len(pos)))
    neg_sel = [neg[j] for j in rng_.choice(len(neg), size=n_same, replace=False)]
    seqs = [pos[i].sequence for i in idxs_p] + [n.sequence for n in neg_sel]
    labels = np.array([1] * len(idxs_p) + [0] * len(neg_sel), dtype=np.float32)
    return seqs, labels
tr_s, tr_y = mk(idx_tr, rng)
va_s, va_y = mk(idx_va, rng)  # consumed for rng fidelity with the benchmark
te_s, te_y = mk(idx_te, rng)
Xtr_full, vocab = build_kmer_features(tr_s)
Xte, _ = build_kmer_features(te_s, vocab=vocab)

FRACS = [0.1, 0.25, 0.5, 0.75, 1.0]
SEEDS = [0, 1, 2, 3]
out = {"fractions": FRACS, "seeds": SEEDS, "n_train_full": len(tr_s),
       "n_test": len(te_s), "runs": []}
for frac in FRACS:
    aucs = []
    for sd in SEEDS:
        r2 = np.random.default_rng(100 + sd)
        k = max(2, int(round(frac * len(tr_s))))
        sel = r2.choice(len(tr_s), size=k, replace=False) if frac < 1.0 else np.arange(len(tr_s))
        lr = LogisticRegression(max_iter=2000, C=1.0, n_jobs=2).fit(
            Xtr_full[sel], tr_y[sel])
        aucs.append(float(M.auc(te_y, Xte @ lr.coef_[0] + lr.intercept_[0])))
    out["runs"].append({"fraction": frac, "n_train": int(round(frac * len(tr_s))),
                        "auc_mean": float(np.mean(aucs)),
                        "auc_sd": float(np.std(aucs)),
                        "aucs": aucs})
    print(frac, round(np.mean(aucs), 4), "+/-", round(np.std(aucs), 4), flush=True)
json.dump(out, open("results/cpp_learning_curve.json", "w"), indent=1)
