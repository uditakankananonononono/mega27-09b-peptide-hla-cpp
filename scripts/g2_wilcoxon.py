"""G2 (locked addendum v2): Wilcoxon signed-rank over per-allele AUC deltas
(ensemble - MHCflurry) from results/per_allele_headtohead.json. Pure stats on
committed numbers; no refit, no rescoring."""
import json
import numpy as np
from scipy import stats

d = json.load(open("results/per_allele_headtohead.json"))
rows = [(a, v["ensemble_auc"] - v["mhcflurry_auc"], v["n"])
        for a, v in d["per_allele"].items()]
deltas = np.array([r[1] for r in rows])
w = stats.wilcoxon(deltas, alternative="greater")
out = {
    "n_alleles": len(rows),
    "ensemble_wins": int((deltas > 0).sum()),
    "median_delta_auc": float(np.median(deltas)),
    "mean_delta_auc": float(deltas.mean()),
    "wilcoxon_stat": float(w.statistic),
    "wilcoxon_p_one_sided": float(w.pvalue),
    "per_allele_deltas": {a: round(dlt, 4) for a, dlt, _ in rows},
    "source": "results/per_allele_headtohead.json (committed)",
}
json.dump(out, open("results/g2_wilcoxon.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "per_allele_deltas"}, indent=1))
