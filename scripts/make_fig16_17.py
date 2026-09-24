"""fig16: motif-logo grid for 12 alleles (binder frequency x info content).
fig17: learned PSSM (ridge) weight heatmaps for 4 alleles - model introspection."""
import sys, json
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.models.pssm import AllelePSSM

AA = "ACDEFGHIKLMNPQRSTVWY"; AAI = {a:i for i,a in enumerate(AA)}
df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples: counts[e.allele] += 1
keep = {a for a,c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

# ---- fig16: logo grid ----
top12 = sorted(keep, key=lambda a:-counts[a])[:12]
def freq_matrix(seqs, L=9):
    F = np.ones((L,20))*0.5
    for s in seqs:
        for p,ch in enumerate(s[:L]):
            if ch in AAI: F[p,AAI[ch]] += 1
    return F/F.sum(axis=1,keepdims=True)
fig, axes = plt.subplots(3, 4, figsize=(11, 6.4))
for ax, a in zip(axes.ravel(), top12):
    bind = [e.sequence for e in train if e.allele==a and e.binder and len(e.sequence)==9]
    P = freq_matrix(bind)
    H = -(P*np.log2(P+1e-12)).sum(axis=1)
    IC = np.log2(20)-H
    # height per residue = p*IC; render top-3 letters per column as text stacks
    for pos in range(9):
        heights = P[pos]*IC[pos]
        order = np.argsort(heights)[::-1][:3]
        y = 0.0
        for r in order:
            ax.text(pos, y+heights[r]/2, AA[r], ha="center", va="center",
                    fontsize=4.5+11*P[pos,r], color="black")
            y += heights[r]
    ax.set_xlim(-0.6, 8.6); ax.set_ylim(0, 4.4)
    ax.set_title(f"{a.replace('HLA-','')} (n={len(bind):,})", fontsize=8)
    ax.set_xticks(range(9)); ax.set_xticklabels([f"P{p+1}" for p in range(9)], fontsize=6)
    ax.set_yticks([0,2,4]); ax.tick_params(labelsize=6)
axes.ravel()[0].set_ylabel("bits", fontsize=7)
fig.suptitle("Sequence logos of 9-mer binders, 12 best-covered alleles (train split)", fontsize=10)
fig.tight_layout(); fig.savefig("paper/figures/fig16_motif_logos.png", dpi=200); plt.close(fig)
print("fig16 done", flush=True)

# ---- fig17: learned PSSM weights ----
top4 = sorted(keep, key=lambda a:-counts[a])[:4]
fig, axes = plt.subplots(2, 2, figsize=(9.5, 7))
for ax, a in zip(axes.ravel(), top4):
    exs = [e for e in train if e.allele==a and len(e.sequence)==9]
    m = AllelePSSM().fit([e.sequence for e in exs], np.array([e.log_ic50 for e in exs]))
    blk = m.W[:9]  # (9, 20) learned weights, 9-mer block
    im = ax.imshow(-blk.T, aspect="auto", cmap="RdBu_r",
                   vmin=-np.abs(blk).max(), vmax=np.abs(blk).max())
    ax.set_xticks(range(9)); ax.set_xticklabels([f"P{p+1}" for p in range(9)], fontsize=7)
    ax.set_yticks(range(20)); ax.set_yticklabels(list(AA), fontsize=6)
    ax.set_title(f"{a.replace('HLA-','')} (n={len(exs):,})", fontsize=9)
fig.suptitle("Learned additive weights (ridge, 9-mers only): red = favorable (lower IC50)", fontsize=10)
fig.tight_layout(); fig.savefig("paper/figures/fig17_pssm_weights.png", dpi=200); plt.close(fig)
print("fig17 done", flush=True)
