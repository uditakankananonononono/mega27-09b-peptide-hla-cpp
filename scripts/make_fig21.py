"""fig21: CNN vs 3-mer LR scores for the 18 named candidates (cross-model
corroboration), with the 0.7 cascade threshold lines."""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

scan = json.load(open("results/cpp_alanine_scan.json"))["candidates"]
named = {c["name"]: c for c in json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]}
x = np.array([named[r["name"]]["p_cpp"] for r in scan])   # CNN cascade score
y = np.array([r["p0"] for r in scan])                     # LR score
fig, ax = plt.subplots(figsize=(5.6, 5.2))
ax.axvline(0.7, color="gray", lw=0.8, ls="--"); ax.axhline(0.7, color="gray", lw=0.8, ls="--")
ax.plot([0.4, 1.02], [0.4, 1.02], color="lightgray", lw=0.8)
colors = ["#b2182b" if yy >= 0.7 else "#2166ac" for yy in y]
ax.scatter(x, y, c=colors, s=42, zorder=3)
for xi, yi, r in zip(x, y, scan):
    if yi >= 0.7 or xi < 0.995:
        ax.annotate(r["name"].replace("9B-CPP-", ""), (xi, yi),
                    textcoords="offset points", xytext=(4, 4), fontsize=7)
ax.set_xlabel("CNN classifier p (cascade model)")
ax.set_ylabel("3-mer LR p (independent model family)")
ax.set_title("Cross-model corroboration of the 18 named candidates")
ax.text(0.42, 0.93, "corroborated by both\n(upper right of 0.7/0.7)", fontsize=8,
        transform=ax.transAxes, va="top")
fig.tight_layout()
fig.savefig("paper/figures/fig21_crossmodel.png", dpi=150)
fig.savefig("paper/figures/fig21_crossmodel.pdf")
both = [r["name"] for r in scan if r["p0"] >= 0.7 and named[r["name"]]["p_cpp"] >= 0.99]
print("corroborated:", both, "mean CNN p:", x.mean().round(4), "mean LR p:", y.mean().round(4))
