"""Fig 24: per-length head-to-head bar chart (ensemble vs MHCflurry)."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("results/per_length_headtohead.json"))
strata = d["strata"]
Ls = sorted(strata.keys(), key=int)
ens = [strata[L]["auc_ensemble"] for L in Ls]
mfl = [strata[L]["auc_mhcflurry"] for L in Ls]
ns = [strata[L]["n"] for L in Ls]
x = np.arange(len(Ls)); w = 0.36
fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.bar(x - w/2, ens, w, label="ensemble (ours)", color="#2166ac")
ax.bar(x + w/2, mfl, w, label="MHCflurry 2.2.1", color="#b2182b")
for i, L in enumerate(Ls):
    ax.text(x[i], 0.62, f"n={ns[i]}", ha="center", fontsize=8)
    dlt = strata[L]["delta"]
    ax.text(x[i], max(ens[i], mfl[i]) + 0.012, f"$\\Delta$={dlt:+.3f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([f"{L}-mer" for L in Ls])
ax.set_ylim(0.55, 1.0); ax.set_ylabel("AUC (held-out h2h subset)")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("paper/figures/fig24.pdf")
print("fig24 written")
