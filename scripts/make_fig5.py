import json, sys
sys.path.insert(0, "src")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

pa = json.load(open("results/per_allele_analysis.json"))["per_allele"]
n = [r["n_train"] for r in pa]
delta = [r["cnn_auc"] - r["pssm_auc"] for r in pa]
fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.scatter(n, delta, s=22, alpha=0.75, edgecolor="k", linewidth=0.3)
coef = np.polyfit(np.log10(n), delta, 1)
xs = np.linspace(min(np.log10(n)), max(np.log10(n)), 50)
ax.plot(10**xs, coef[0]*xs + coef[1], "r--", lw=1.2,
        label=f"log-linear fit (slope {coef[0]:.3f})")
ax.axhline(0, color="gray", lw=0.8)
ax.set_xscale("log")
ax.set_xlabel("training pairs (allele, log scale)")
ax.set_ylabel("CNN AUC - PSSM AUC")
ax.set_title("Non-additive gain grows with per-allele data")
ax.legend()
fig.tight_layout(); fig.savefig("paper/figures/fig5_datadependence.png")
print("fig5 written")
