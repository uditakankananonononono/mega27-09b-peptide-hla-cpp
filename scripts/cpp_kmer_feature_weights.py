"""CPP k-mer LR feature-weight analysis: which 3-mers drive CPP classification.

Mirrors train_cpp.py data prep exactly (same rng seed 11 and call order) so the
refit LR reproduces the reported benchmark split; coefficients are then read off
and grouped by amino-acid chemistry. Saves results/cpp_kmer_feature_weights.json.
"""
from __future__ import annotations

import json

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
tr_s, tr_y = mk(idx_tr); va_s, va_y = mk(idx_va); te_s, te_y = mk(idx_te)

Xtr, vocab = build_kmer_features(tr_s)
lr = LogisticRegression(max_iter=2000, C=1.0, n_jobs=2).fit(Xtr, tr_y)
Xte, _ = build_kmer_features(te_s, vocab=vocab)
scores = Xte @ lr.coef_[0] + lr.intercept_[0]
test_auc = M.auc(te_y, scores)
test_auc
inv = {i: k for k, i in vocab.items()}
w = lr.coef_[0]
order_w = np.argsort(-w)
top_pos = [{"kmer": inv[i], "weight": float(w[i])} for i in order_w[:25]]
top_neg = [{"kmer": inv[i], "weight": float(w[i])} for i in order_w[::-1][:15]]

def contains_any(s, chars): return any(c in s for c in chars)

BASIC, AROM, HPHOB = "RK", "WYF", "LIV"
n50 = 50
top50 = [inv[i] for i in order_w[:n50]]
groups = {
    "top50_containing_R_or_K": sum(contains_any(k, BASIC) for k in top50),
    "top50_containing_W_Y_F": sum(contains_any(k, AROM) for k in top50),
    "top50_containing_L_I_V": sum(contains_any(k, HPHOB) for k in top50),
}
kmers = np.array([inv[i] for i in range(len(inv))])
def mean_weight(mask_chars):
    sel = np.array([contains_any(k, mask_chars) for k in kmers])
    return {"n_kmers": int(sel.sum()), "mean_weight": float(w[sel].mean()),
            "mean_weight_rest": float(w[~sel].mean())}
chem = {"contains_R": mean_weight("R"), "contains_K": mean_weight("K"),
        "contains_RK": mean_weight(BASIC), "contains_W": mean_weight("W"),
        "contains_WYF": mean_weight(AROM), "contains_LIV": mean_weight(HPHOB),
        "contains_DE": mean_weight("DE")}

json.dump({
    "refit_test_auc": float(test_auc),
    "benchmark_test_auc": 0.9265232974910395,
    "n_train": len(tr_s), "n_test": len(te_s), "vocab_size": len(vocab),
    "intercept": float(lr.intercept_[0]),
    "top_positive": top_pos, "top_negative": top_neg,
    "top50_composition": groups, "chemistry_effects": chem,
}, open("results/cpp_kmer_feature_weights.json", "w"), indent=1)
print(json.dumps({"refit_test_auc": float(test_auc), "vocab": len(vocab),
                  "top5": [d["kmer"] for d in top_pos[:5]],
                  "top50_RK": groups["top50_containing_R_or_K"]}, indent=1))
