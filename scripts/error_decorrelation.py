"""Direct measurement of error decorrelation between the ensemble and
MHCflurry on the 6,000-peptide head-to-head subset.

If the ensemble's AUC win came from modeling something MHCflurry knows,
the two models would misorder the same positive-negative pairs. If it
comes from decorrelated errors, the ensemble should order a substantial
fraction of MHCflurry's misordered pairs correctly (and vice versa).
Output: results/error_decorrelation.json
"""
import json
import numpy as np
from scipy.stats import spearmanr

d = np.load("results/h2h_scores.npz")
y = d["labels"]
mfl, ens, cnn, pssm = d["mhcflurry"], d["ensemble"], d["cnn"], d["pssm"]

def auc(scores):
    pos, neg = scores[y], scores[~y]
    c = (pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()
    return float(c)

def pair_stats(s_a, s_b):
    """Among (pos,neg) pairs that model A misorders, fraction model B
    orders correctly; plus the reverse. Ties count as half."""
    pos_a, neg_a = s_a[y], s_a[~y]
    pos_b, neg_b = s_b[y], s_b[~y]
    wrong_a = (pos_a[:, None] < neg_a[None, :])
    right_b = (pos_b[:, None] > neg_b[None, :])
    tie_b = (pos_b[:, None] == neg_b[None, :])
    n_wrong_a = int(wrong_a.sum())
    frac_b_right_given_a_wrong = float((right_b[wrong_a].sum() + 0.5 * tie_b[wrong_a].sum()) / n_wrong_a)
    return n_wrong_a, frac_b_right_given_a_wrong

n_wrong_mfl, ens_right = pair_stats(mfl, ens)
n_wrong_ens, mfl_right = pair_stats(ens, mfl)
n_pairs = int(y.sum()) * int((~y).sum())

out = {
    "n_peptides": int(len(y)),
    "n_pairs": n_pairs,
    "auc_mhcflurry": auc(mfl),
    "auc_ensemble": auc(ens),
    "mfl_misordered_pairs": n_wrong_mfl,
    "mfl_misordered_frac": n_wrong_mfl / n_pairs,
    "ens_correct_given_mfl_wrong": ens_right,
    "ens_misordered_pairs": n_wrong_ens,
    "mfl_correct_given_ens_wrong": mfl_right,
    "spearman_mfl_ens_all": float(spearmanr(mfl, ens).statistic),
    "spearman_mfl_ens_pos": float(spearmanr(mfl[y], ens[y]).statistic),
    "spearman_mfl_ens_neg": float(spearmanr(mfl[~y], ens[~y]).statistic),
    "spearman_cnn_pssm_pos": float(spearmanr(cnn[y], pssm[y]).statistic),
    "spearman_cnn_pssm_neg": float(spearmanr(cnn[~y], pssm[~y]).statistic),
}
json.dump(out, open("results/error_decorrelation.json", "w"), indent=1)
for k, v in out.items():
    print(k, round(v, 4) if isinstance(v, float) else v)
