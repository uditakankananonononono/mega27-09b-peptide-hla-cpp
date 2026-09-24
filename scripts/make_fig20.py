"""fig20: CPP k-mer LR feature weights - top 20 positive + top 10 negative 3-mers,
colored by chemistry (basic R/K, aromatic W/Y/F, hydrophobic L/I/V, other)."""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("results/cpp_kmer_feature_weights.json"))
pos, neg = d["top_positive"][:20], d["top_negative"][:10]
items = [(x["kmer"], x["weight"]) for x in pos] + [(x["kmer"], x["weight"]) for x in neg]

def color(k):
    if any(c in k for c in "RK"): return "#b2182b"
    if any(c in k for c in "WYF"): return "#2166ac"
    if any(c in k for c in "LIV"): return "#1b7837"
    return "#7f7f7f"

labels = [k for k, _ in items][::-1]
vals = [w for _, w in items][::-1]
cols = [color(k) for k in labels]
fig, ax = plt.subplots(figsize=(7.2, 6.4))
ax.barh(range(len(vals)), vals, color=cols)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(labels, fontsize=8)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("logistic-regression coefficient (log-odds per 3-mer occurrence)")
ax.set_title("CPP 3-mer classifier: top 20 positive and top 10 negative feature weights")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#b2182b", label="contains R/K (basic)"),
                   Patch(color="#2166ac", label="contains W/Y/F (aromatic)"),
                   Patch(color="#1b7837", label="contains L/I/V (hydrophobic)"),
                   Patch(color="#7f7f7f", label="other")],
          loc="lower right", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig("paper/figures/fig20_cpp_kmer_weights.pdf")
fig.savefig("paper/figures/fig20_cpp_kmer_weights.png", dpi=150)
print("fig20 saved;", sum(any(c in k for c in 'RK') for k in labels[:10]), "of top 10 negative contain R/K")
