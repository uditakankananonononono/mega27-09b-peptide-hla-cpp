"""Fig 23: error decorrelation between ensemble and MHCflurry."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("results/h2h_scores.npz")
dec = json.load(open("results/error_decorrelation.json"))
y = d["labels"]
mfl, ens = d["mhcflurry"], d["ensemble"]

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
ax = axes[0]
rng = np.random.default_rng(7)
idx = rng.choice(len(y), 4000, replace=False)
ax.scatter(mfl[idx][~y[idx]], ens[idx][~y[idx]], s=3, alpha=0.25, c="steelblue", label="non-binder")
ax.scatter(mfl[idx][y[idx]], ens[idx][y[idx]], s=3, alpha=0.25, c="firebrick", label="binder")
ax.set_xlabel("MHCflurry score")
ax.set_ylabel("ensemble score")
ax.set_title(f"(a) per-peptide scores (Spearman {dec['spearman_mfl_ens_all']:.3f})")
ax.legend(markerscale=3, fontsize=8, frameon=False)

ax = axes[1]
labels = ["MHCflurry wrong:\nensemble right", "ensemble wrong:\nMHCflurry right"]
vals = [dec["ens_correct_given_mfl_wrong"], dec["mfl_correct_given_ens_wrong"]]
bars = ax.bar(labels, vals, color=["seagreen", "gray"], width=0.5)
ax.axhline(0.5, ls="--", c="k", lw=0.8)
ax.text(1.28, 0.505, "chance", fontsize=8)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.008, f"{v:.3f}", ha="center", fontsize=9)
ax.set_ylim(0, 0.68)
ax.set_ylabel("fraction of misordered pairs\nordered correctly by the other model")
ax.set_title("(b) asymmetric error correction")
fig.tight_layout()
fig.savefig("paper/figures/fig23_error_decorrelation.pdf")
print("fig23 saved")
