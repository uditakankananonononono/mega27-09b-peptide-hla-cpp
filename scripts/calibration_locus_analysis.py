"""Murphy Brier-score decomposition from the calibration bins + locus-level AUC.
Outputs: results/calibration_decomposition.json, results/locus_auc.json
"""
import json
import numpy as np

cal = json.load(open("results/calibration.json"))
bins = cal["prob_bins"]
N = sum(b["n"] for b in bins)
obar = sum(b["n"] * b["empirical"] for b in bins) / N
rel = sum(b["n"] * (b["mean_pred"] - b["empirical"]) ** 2 for b in bins) / N
res = sum(b["n"] * (b["empirical"] - obar) ** 2 for b in bins) / N
unc = obar * (1 - obar)
brier = rel - res + unc
json.dump({
    "n": N, "base_rate": obar,
    "brier_binned": brier, "reliability": rel, "resolution": res, "uncertainty": unc,
    "brier_baseline": unc,
    "skill_vs_climatology": 1 - brier / unc,
}, open("results/calibration_decomposition.json", "w"), indent=1)
print("brier", round(brier, 5), "rel", round(rel, 5), "res", round(res, 5),
      "unc", round(unc, 5), "skill", round(1 - brier / unc, 4))

d = np.load("results/ensemble_test_predictions.npz")
labels, ens, alleles = d["labels"], d["ensemble"], d["alleles"]
sys_auc = None
from sklearn.metrics import roc_auc_score
out = {"loci": []}
for locus in ["HLA-A", "HLA-B", "HLA-C"]:
    sel = np.array([a.startswith(locus) for a in alleles])
    if sel.sum() > 50 and len(set(labels[sel])) > 1:
        out["loci"].append({
            "locus": locus, "n": int(sel.sum()),
            "n_binders": int(labels[sel].sum()),
            "ensemble_auc": float(roc_auc_score(labels[sel], ens[sel])),
        })
json.dump(out, open("results/locus_auc.json", "w"), indent=1)
print(json.dumps(out, indent=1))
