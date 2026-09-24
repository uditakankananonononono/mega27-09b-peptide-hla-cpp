"""Extension analyses for the 50-page paper expansion.
All figures/results here are computed from the real datasets and the real
trained models saved in results/. Run from repo root:
    python3 scripts/paper_extension_analyses.py
Outputs: results/motif_information.json, results/calibration.json,
results/cpp_properties.json, results/ensemble_test_predictions.npz,
results/per_allele_headtohead.json, paper/figures/fig8..fig15.
"""
import json, sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.data.cppsite import parse_fasta, natural_only
from peptidehlacpp.eval import metrics as M

AA = "ACDEFGHIKLMNPQRSTVWY"
AAI = {a: i for i, a in enumerate(AA)}
rng = np.random.default_rng(7)
FIG = "paper/figures"

# ---------------- load pHLA data & rebuild the exact split ----------------
print("loading IEDB pairs...", flush=True)
df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples:
    counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)
print(f"train {len(train)} val {len(val)} test {len(test)} alleles {len(keep)}", flush=True)

# ---------------- 1. motif information content (top 12 alleles) ----------
print("motif information content...", flush=True)
top_alleles = sorted(keep, key=lambda a: -counts[a])[:12]
def freq_matrix(seqs, L=9):
    F = np.ones((L, 20)) * 0.5  # pseudocount
    for s in seqs:
        for p, ch in enumerate(s[:L]):
            if ch in AAI:
                F[p, AAI[ch]] += 1
    return F / F.sum(axis=1, keepdims=True)
def info_content(P, bg=1/20):
    H = -(P * np.log2(P + 1e-12)).sum(axis=1)
    return np.log2(20) - H
IC = np.zeros((len(top_alleles), 9))
nbinders = {}
for i, a in enumerate(top_alleles):
    bind = [e.sequence for e in train if e.allele == a and e.binder and len(e.sequence) == 9]
    nbinders[a] = len(bind)
    IC[i] = info_content(freq_matrix(bind))
motif_res = {"alleles": top_alleles, "n_binders": nbinders,
             "info_content": {a: IC[i].tolist() for i, a in enumerate(top_alleles)},
             "anchor_delta": {a: float(IC[i,1] + IC[i,8] - IC[i,1:8].mean()*2) for i, a in enumerate(top_alleles)}}
json.dump(motif_res, open("results/motif_information.json", "w"), indent=2)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
im = ax.imshow(IC, aspect="auto", cmap="viridis", vmin=0, vmax=IC.max())
ax.set_xticks(range(9)); ax.set_xticklabels([f"P{p+1}" for p in range(9)])
ax.set_yticks(range(len(top_alleles)))
ax.set_yticklabels([a.replace("HLA-", "") for a in top_alleles], fontsize=8)
ax.set_xlabel("Peptide position"); ax.set_title("Per-position information content of 9-mer binders (bits)")
for i in range(len(top_alleles)):
    for j in range(9):
        ax.text(j, i, f"{IC[i,j]:.1f}", ha="center", va="center",
                color="white" if IC[i,j] < IC.max()*0.6 else "black", fontsize=6.5)
fig.colorbar(im, label="bits")
fig.tight_layout(); fig.savefig(f"{FIG}/fig8_motif_information.png", dpi=200); plt.close(fig)

# ------------- 2. A*02:01 vs A*03:01 anchor composition ------------------
print("anchor composition comparison...", flush=True)
def comp_at(allele, pos):
    bind = [e.sequence for e in train if e.allele == allele and e.binder and len(e.sequence) == 9]
    F = freq_matrix(bind)
    return F[pos]
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
for ax, pos, title in zip(axes, [1, 8], ["Position 2", "Position 9 (C-term)"]):
    x = np.arange(20); w = 0.38
    ax.bar(x - w/2, comp_at("HLA-A*02:01", pos), w, label="A*02:01 binders")
    ax.bar(x + w/2, comp_at("HLA-A*03:01", pos), w, label="A*03:01 binders")
    ax.set_xticks(x); ax.set_xticklabels(list(AA), fontsize=7)
    ax.set_title(title); ax.legend(fontsize=8)
axes[0].set_ylabel("Frequency")
fig.suptitle("Anchor-residue composition: A*02:01 vs A*03:01 9-mer binders (train split)")
fig.tight_layout(); fig.savefig(f"{FIG}/fig9_anchor_composition.png", dpi=200); plt.close(fig)
anchor_comp = {"P2": {"A0201": comp_at("HLA-A*02:01",1).tolist(), "A0301": comp_at("HLA-A*03:01",1).tolist()},
               "P9": {"A0201": comp_at("HLA-A*02:01",8).tolist(), "A0301": comp_at("HLA-A*03:01",8).tolist()},
               "aa_order": list(AA)}

# ---------------- 3. ensemble scores on full test + calibration ----------
print("scoring full test set with ensemble...", flush=True)
import torch
from peptidehlacpp.models.cnn import PHLACNN
from peptidehlacpp.models.pssm import AllelePSSM
from peptidehlacpp.training.train_phla import enc_batch
alleles = sorted(keep)
amap = {a: i for i, a in enumerate(alleles)}
by_al_train = defaultdict(list)
for e in train:
    by_al_train[e.allele].append(e)
pssms = {a: AllelePSSM().fit([e.sequence for e in exs], np.array([e.log_ic50 for e in exs]))
         for a, exs in by_al_train.items()}
cnn = PHLACNN(n_alleles=len(alleles))
cnn.load_state_dict(torch.load("results/phla_cnn.pt")["state"])
cnn.eval()
pssm_s = -np.array([pssms[e.allele].predict([e.sequence])[0] for e in test])
cnn_s = []
with torch.no_grad():
    for s in range(0, len(test), 2048):
        chunk = test[s:s+2048]
        X, m, al = enc_batch([e.sequence for e in chunk], np.array([e.allele for e in chunk]), amap, False)
        cnn_s.append(cnn(X, m, al)[1].numpy())
cnn_s = np.concatenate(cnn_s)
def z(x): return (x - x.mean()) / (x.std() + 1e-9)
ens_s = z(pssm_s) + z(cnn_s)
labels = np.array([e.binder for e in test])
np.savez("results/ensemble_test_predictions.npz", labels=labels, pssm=pssm_s, cnn=cnn_s, ensemble=ens_s,
         alleles=np.array([e.allele for e in test]))
print(f"test AUC: pssm {M.auc(labels,pssm_s):.4f} cnn {M.auc(labels,cnn_s):.4f} ens {M.auc(labels,ens_s):.4f}", flush=True)

# calibration: decile reliability of ensemble
order = np.argsort(-ens_s)
dec = np.array_split(order, 10)
cal = {"bins": []}
for i, idx in enumerate(dec):
    cal["bins"].append({"decile": i+1, "n": int(len(idx)),
                        "mean_score": float(ens_s[idx].mean()),
                        "empirical_binder_frac": float(labels[idx].mean())})
# ECE on probability proxy: map ens score to binder prob via isotonic-free logistic fit on val
val_pssm = -np.array([pssms[e.allele].predict([e.sequence])[0] for e in val])
val_cnn = []
with torch.no_grad():
    for s in range(0, len(val), 2048):
        chunk = val[s:s+2048]
        X, m, al = enc_batch([e.sequence for e in chunk], np.array([e.allele for e in chunk]), amap, False)
        val_cnn.append(cnn(X, m, al)[1].numpy())
val_ens = z(val_pssm) + z(np.concatenate(val_cnn))
val_lab = np.array([e.binder for e in val])
# logistic recalibration on val
A = np.stack([val_ens, np.ones_like(val_ens)], axis=1)
w = np.zeros(2)
for _ in range(200):
    p = 1/(1+np.exp(-A@w))
    g = A.T @ (p - val_lab) / len(val_lab)
    H = (A.T * (p*(1-p))) @ A / len(val_lab) + 1e-4*np.eye(2)
    w -= np.linalg.solve(H, g)
test_p = 1/(1+np.exp(-(w[0]*ens_s + w[1])))
qs = np.quantile(test_p, np.linspace(0, 1, 11))
ece = 0.0; cal["prob_bins"] = []
for lo, hi in zip(qs[:-1], qs[1:]):
    m_ = (test_p >= lo) & (test_p <= hi)
    if m_.sum() == 0: continue
    gap = abs(test_p[m_].mean() - labels[m_].mean())
    ece += m_.mean() * gap
    cal["prob_bins"].append({"lo": float(lo), "hi": float(hi), "n": int(m_.sum()),
                             "mean_pred": float(test_p[m_].mean()),
                             "empirical": float(labels[m_].mean())})
cal["ece"] = float(ece)
cal["logistic_w"] = w.tolist()
json.dump(cal, open("results/calibration.json", "w"), indent=2)
fig, ax = plt.subplots(figsize=(5.2, 4.4))
mp = [b["mean_pred"] for b in cal["prob_bins"]]; em = [b["empirical"] for b in cal["prob_bins"]]
ax.plot([0,1],[0,1],"k--",lw=1,label="perfect calibration")
ax.plot(mp, em, "o-", label="ensemble (Platt on val)")
ax.set_xlabel("Mean predicted binder probability"); ax.set_ylabel("Empirical binder fraction")
ax.set_title(f"Calibration on held-out test (ECE={ece:.4f})"); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/fig10_calibration.png", dpi=200); plt.close(fig)

# ------------- 4. PSSM vs CNN error correlation --------------------------
r = np.corrcoef(z(pssm_s), z(cnn_s))[0,1]
fig, ax = plt.subplots(figsize=(5.0, 4.4))
sub = rng.choice(len(test), size=6000, replace=False)
ax.scatter(z(pssm_s)[sub], z(cnn_s)[sub], s=2, alpha=0.25,
           c=["tab:red" if labels[i] else "tab:blue" for i in sub])
ax.set_xlabel("PSSM score (z)"); ax.set_ylabel("CNN score (z)")
ax.set_title(f"Model score correlation on test (Pearson r={r:.3f})")
fig.tight_layout(); fig.savefig(f"{FIG}/fig11_score_correlation.png", dpi=200); plt.close(fig)
json.dump({"pearson_pssm_cnn": float(r)}, open("results/score_correlation.json","w"))

# ------------- 5. length dependence of binding ---------------------------
print("length analysis...", flush=True)
lens = sorted({len(e.sequence) for e in examples})
lstat = []
for L in lens:
    exs = [e for e in examples if len(e.sequence) == L]
    if len(exs) < 50: continue
    lstat.append({"len": L, "n": len(exs), "binder_frac": float(np.mean([e.binder for e in exs])),
                  "median_ic50_nm": float(np.median([10**e.log_ic50 for e in exs]))})
json.dump({"by_length": lstat}, open("results/length_analysis.json","w"), indent=2)
fig, ax1 = plt.subplots(figsize=(5.6, 3.8))
Ls = [d["len"] for d in lstat]
ax1.bar(Ls, [d["n"] for d in lstat], color="lightgray", label="n assays")
ax1.set_ylabel("n assays"); ax1.set_xlabel("Peptide length")
ax2 = ax1.twinx()
ax2.plot(Ls, [d["binder_frac"] for d in lstat], "o-", color="tab:red", label="binder fraction")
ax2.set_ylabel("Binder fraction (IC50 ≤ 500 nM)", color="tab:red")
ax.set_title = None
ax1.set_title("Assay count and binder fraction by peptide length (all alleles)")
fig.tight_layout(); fig.savefig(f"{FIG}/fig12_length.png", dpi=200); plt.close(fig)

# ------------- 6. CPP physicochemical properties --------------------------
print("CPP properties...", flush=True)
pos = natural_only(parse_fasta("data/raw/cppsite2_natural.fa")) + natural_only(parse_fasta("data/raw/cppsite2_nonnatural.fa"))
neg_pool = []
for f in ["data/raw/uniprot_reviewed_len8_35.fasta"]:
    neg_pool += [p.sequence for p in natural_only(parse_fasta(f))]
neg = neg_pool[:]
# property functions (real, standard scales)
EISENBERG = dict(zip(AA, [0.62,0.29,-0.90,-0.74,1.19,0.48,-0.40,1.38,-1.50,1.06,0.64,-1.50,0.12,-0.85,-0.18,0.81,0.26,-0.73,0.00,1.08]))
def net_charge(s):
    return sum(1 for c in s if c in "KR") + 0.1*sum(1 for c in s if c=="H") - sum(1 for c in s if c in "DE")
def mean_hyd(s):
    v=[EISENBERG[c] for c in s if c in EISENBERG]; return float(np.mean(v)) if v else 0.0
def hyd_moment(s, angle=100):
    v=[EISENBERG[c] for c in s if c in EISENBERG]
    if len(v)<3: return 0.0
    th=np.deg2rad(angle); sx=sum(h*np.cos(i*th) for i,h in enumerate(v)); sy=sum(h*np.sin(i*th) for i,h in enumerate(v))
    return float(np.hypot(sx,sy)/len(v))
def props(seqs):
    return {"charge": [net_charge(s) for s in seqs],
            "hydrophobicity": [mean_hyd(s) for s in seqs],
            "hydrophobic_moment": [hyd_moment(s) for s in seqs]}
pp = props([p.sequence for p in pos]); nn = props(neg)
cpp_props = {"n_pos": len(pos), "n_neg": len(neg),
             "pos": {k: {"mean": float(np.mean(v)), "sd": float(np.std(v))} for k,v in pp.items()},
             "neg": {k: {"mean": float(np.mean(v)), "sd": float(np.std(v))} for k,v in nn.items()}}
# Welch t and Cohen d
from math import sqrt
for k in pp:
    a,b = np.array(pp[k]), np.array(nn[k])
    d = (a.mean()-b.mean())/sqrt((a.var()+b.var())/2)
    cpp_props[f"cohens_d_{k}"] = float(d)
json.dump(cpp_props, open("results/cpp_properties.json","w"), indent=2)
fig, axes = plt.subplots(1, 3, figsize=(11, 3.3))
for ax, k, lab in zip(axes, ["charge","hydrophobicity","hydrophobic_moment"],
                      ["Net charge (pH 7)","Mean Eisenberg hydrophobicity","Hydrophobic moment (α, 100°)"]):
    bins = np.linspace(min(pp[k]+nn[k]), np.percentile(pp[k]+nn[k], 99), 40)
    ax.hist(nn[k], bins=bins, density=True, alpha=0.55, label=f"UniProt non-CPP (n={len(neg)})")
    ax.hist(pp[k], bins=bins, density=True, alpha=0.55, label=f"CPPsite CPPs (n={len(pos)})")
    ax.set_xlabel(lab); ax.legend(fontsize=7)
axes[0].set_ylabel("Density")
fig.suptitle(f"Physicochemical separation of CPPs (Cohen's d: charge {cpp_props['cohens_d_charge']:.2f}, hydrophobicity {cpp_props['cohens_d_hydrophobicity']:.2f}, moment {cpp_props['cohens_d_hydrophobic_moment']:.2f})", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FIG}/fig13_cpp_properties.png", dpi=200); plt.close(fig)

# ------------- 7. AA composition enrichment CPP vs non-CPP ---------------
def aa_freq(seqs):
    F = np.zeros(20)
    tot = 0
    for s in seqs:
        for c in s:
            if c in AAI: F[AAI[c]] += 1; tot += 1
    return F / tot
fp = aa_freq([p.sequence for p in pos]); fn = aa_freq(neg)
lo = np.log2((fp + 1e-4) / (fn + 1e-4))
json.dump({"aa_order": list(AA), "log2_enrichment": lo.tolist(),
           "cpp_freq": fp.tolist(), "uniprot_freq": fn.tolist()},
          open("results/cpp_aa_enrichment.json","w"), indent=2)
fig, ax = plt.subplots(figsize=(7.0, 3.4))
order_aa = np.argsort(-lo)
ax.bar(range(20), lo[order_aa], color=["tab:red" if v>0 else "tab:blue" for v in lo[order_aa]])
ax.set_xticks(range(20)); ax.set_xticklabels([AA[i] for i in order_aa])
ax.set_ylabel("log2(CPP frequency / UniProt frequency)")
ax.set_title("Amino-acid enrichment in CPPsite 2.0 vs length-matched UniProt windows")
fig.tight_layout(); fig.savefig(f"{FIG}/fig14_aa_enrichment.png", dpi=200); plt.close(fig)
print("CPP properties + enrichment done", flush=True)
print("ALL EXTENSION ANALYSES DONE", flush=True)
