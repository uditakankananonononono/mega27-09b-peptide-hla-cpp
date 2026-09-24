"""Does anchor-epistasis strength predict where the ensemble beats MHCflurry?

Joins per-allele head-to-head win margins (ensemble AUC - MHCflurry AUC,
results/per_allele_headtohead.json) with anchor-epistasis summaries
(results/anchor_epistasis.json) on shared alleles; Spearman + Pearson with a
peptide-count note. Saves results/epistasis_winmargin_link.json.
"""
import json
import numpy as np
from scipy.stats import spearmanr, pearsonr

epi = json.load(open("results/anchor_epistasis.json"))
h2h = json.load(open("results/per_allele_headtohead.json"))
rows_h = h2h["per_allele"]
margin = {}
ntest = {}
for a, r in rows_h.items():
    if r.get("ensemble_auc") is not None and r.get("mhcflurry_auc") is not None:
        margin[a] = r["ensemble_auc"] - r["mhcflurry_auc"]
        ntest[a] = r.get("n", 0)
common = sorted(set(margin) & set(epi))
x = np.array([epi[a]["rms_interaction"] for a in common])
g = np.array([abs(epi[a]["hh_contrast"]) for a in common])
y = np.array([margin[a] for a in common])
n = np.array([ntest[a] for a in common])
res = {
    "n_alleles": len(common),
    "alleles": common,
    "spearman_rms_vs_margin": {"rho": float(spearmanr(x, y).statistic),
                               "p": float(spearmanr(x, y).pvalue)},
    "pearson_rms_vs_margin": {"r": float(pearsonr(x, y).statistic),
                              "p": float(pearsonr(x, y).pvalue)},
    "spearman_abshh_vs_margin": {"rho": float(spearmanr(g, y).statistic),
                                 "p": float(spearmanr(g, y).pvalue)},
    "spearman_ntest_vs_margin": {"rho": float(spearmanr(n, y).statistic),
                                 "p": float(spearmanr(n, y).pvalue)},
    "margin_mean": float(y.mean()), "margin_sd": float(y.std()),
}
json.dump(res, open("results/epistasis_winmargin_link.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "alleles"}, indent=1))
