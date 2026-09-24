"""fig22: CPP 3-mer LR learning curve (log-x) with conservative last-octave
slope extrapolation to a 2x-data point."""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("results/cpp_learning_curve.json"))
mu = np.array([r["auc_mean"] for r in d["runs"]])
sd = np.array([r["auc_sd"] for r in d["runs"]])
n = np.array([r["n_train"] for r in d["runs"]], dtype=float)
fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.errorbar(n, mu, yerr=sd, marker="o", capsize=3, color="#b2182b", label="CPP 3-mer LR (4 seeds)")
ax.set_xscale("log"); ax.set_xlabel("training peptides (log scale)")
ax.set_ylabel("held-out test AUC"); ax.set_ylim(0.75, 0.99)
ax.axhline(0.9265, color="gray", lw=0.7, ls=":")
# conservative: last-octave slope (75%->100%), extrapolate one doubling
slope_last = (mu[-1] - mu[-2]) / np.log2(n[-1] / n[-2])
x2 = 2 * n[-1]; y2 = mu[-1] + slope_last
ax.plot([n[-1], x2], [mu[-1], y2], ls="--", color="#b2182b", alpha=0.6)
ax.scatter([x2], [y2], marker="o", facecolors="none", edgecolors="#b2182b")
ax.annotate(f"2x data ~ {y2:.3f}\n(last-octave slope {slope_last:.3f}/doubling)",
            (x2, y2), textcoords="offset points", xytext=(-118, 8), fontsize=8)
ax.set_title("CPP classifier is data-limited: still climbing at 100% of CPPsite")
ax.legend(loc="lower right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig("paper/figures/fig22_cpp_learning_curve.png", dpi=150)
fig.savefig("paper/figures/fig22_cpp_learning_curve.pdf")
print("last-octave slope/doubling:", round(slope_last, 4), "2x estimate:", round(y2, 4))
