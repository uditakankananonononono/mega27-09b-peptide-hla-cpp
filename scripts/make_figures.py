"""Generate paper figures from results JSONs (real outputs only)."""
import json
import sys
sys.path.insert(0, "src")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.features import hydrophobic_moment, net_charge

phla = json.load(open("results/phla_benchmark.json"))
cpp = json.load(open("results/cpp_benchmark.json"))
designs = json.load(open("results/cpp_designs.json"))

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

# Fig 1: model comparison bars (pHLA)
models = ["pssm", "cnn", "gnn"]
metrics = ["auc", "auc0.1", "ppv", "srcc"]
fig, ax = plt.subplots(figsize=(7, 3.6))
x = np.arange(len(metrics)); w = 0.25
for k, m in enumerate(models):
    ax.bar(x + (k - 1) * w, [phla[m][mm] for mm in metrics], w, label=m.upper())
ax.set_xticks(x); ax.set_xticklabels(["AUC", "AUC0.1", "PPV", "SRCC"])
ax.set_ylim(0, 1); ax.legend(); ax.set_title("pHLA-I binding: model comparison (held-out test)")
fig.tight_layout(); fig.savefig("paper/figures/fig1_phla_models.png"); plt.close(fig)

# Fig 2: CPP classifier + baseline
fig, ax = plt.subplots(figsize=(4.2, 3.4))
labels = ["3-mer LR", "CNN"]
aucs = [cpp["kmer_lr"]["auc"], cpp["cnn"]["auc"]]
ppvs = [cpp["kmer_lr"]["ppv"], cpp["cnn"]["ppv"]]
x = np.arange(2); w = 0.35
ax.bar(x - w/2, aucs, w, label="AUC"); ax.bar(x + w/2, ppvs, w, label="PPV")
ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylim(0, 1); ax.legend()
ax.set_title("CPP classification (held-out test)")
fig.tight_layout(); fig.savefig("paper/figures/fig2_cpp_classifier.png"); plt.close(fig)

# Fig 3: generated candidates property map
top = designs["top_candidates"]
charges = [c["net_charge"] for c in top]
moments = [c["hydrophobic_moment"] for c in top]
scores = [c["p_cpp"] for c in top]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
sc = ax.scatter(charges, moments, c=scores, cmap="viridis", s=28, edgecolor="k", linewidth=0.3)
fig.colorbar(sc, label="p(CPP)")
ax.set_xlabel("net charge (pH 7)"); ax.set_ylabel("hydrophobic moment (delta=100)")
ax.set_title(f"Generated CPP candidates passing cascade (top {len(top)})")
fig.tight_layout(); fig.savefig("paper/figures/fig3_cpp_designs.png"); plt.close(fig)

# Fig 4: allele coverage
import pandas as pd
df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
ex = aggregate(df)
from collections import Counter
cnt = Counter(e.allele for e in ex)
top15 = cnt.most_common(15)
fig, ax = plt.subplots(figsize=(7, 3.6))
ax.bar([a.replace("HLA-", "") for a, _ in top15], [c for _, c in top15])
ax.set_ylabel("unique peptide-allele pairs"); ax.set_title("IEDB class-I data coverage (top 15 alleles)")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
fig.tight_layout(); fig.savefig("paper/figures/fig4_coverage.png"); plt.close(fig)
print("figures written")
