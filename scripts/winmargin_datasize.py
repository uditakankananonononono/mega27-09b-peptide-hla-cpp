"""Is the ensemble win over MHCflurry carried by data-rich alleles?

Joins per-allele head-to-head margins with per-allele train sizes and
correlates them. Output: results/winmargin_datasize.json
"""
import json
import numpy as np
from scipy.stats import spearmanr

h2h = json.load(open("results/per_allele_headtohead.json"))["per_allele"]
ana = {r["allele"]: r for r in json.load(open("results/per_allele_analysis.json"))["per_allele"]}

rows = []
for a, v in h2h.items():
    if a in ana:
        rows.append({
            "allele": a,
            "n_train": ana[a]["n_train"],
            "n_test": v["n"],
            "delta_auc": v["ensemble_auc"] - v["mhcflurry_auc"],
        })
rows.sort(key=lambda r: -r["n_train"])
d = np.array([r["delta_auc"] for r in rows])
ntr = np.array([r["n_train"] for r in rows])
nte = np.array([r["n_test"] for r in rows])

half = len(rows) // 2
out = {
    "n_alleles": len(rows),
    "spearman_delta_vs_log_ntrain": float(spearmanr(d, np.log10(ntr)).statistic),
    "spearman_delta_vs_log_ntrain_p": float(spearmanr(d, np.log10(ntr)).pvalue),
    "spearman_delta_vs_log_ntest": float(spearmanr(d, np.log10(nte)).statistic),
    "mean_delta_top_half_train": float(d[:half].mean()),
    "mean_delta_bottom_half_train": float(d[half:].mean()),
    "wins_top_half": int((d[:half] > 0).sum()),
    "wins_bottom_half": int((d[half:] > 0).sum()),
    "half": half,
    "n_alleles_odd": len(rows) - 2 * half,
    "top_half_alleles": [r["allele"] for r in rows[:half]],
    "bottom_half_alleles": [r["allele"] for r in rows[half:]],
}
json.dump(out, open("results/winmargin_datasize.json", "w"), indent=1)
for k, v in out.items():
    if not k.endswith("alleles"):
        print(k, round(v, 4) if isinstance(v, float) else v)
