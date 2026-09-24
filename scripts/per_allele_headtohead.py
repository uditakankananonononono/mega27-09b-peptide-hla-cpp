"""Per-allele head-to-head: MHCflurry vs our ensemble on the identical
6000-pair held-out subset (same rng seed as head_to_head_mhcflurry.py)."""
import json, sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples: counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

from mhcflurry import Class1AffinityPredictor
pred = Class1AffinityPredictor.load()
supported = set(pred.supported_alleles)
test_supported = [e for e in test if e.allele in supported]
rng = np.random.default_rng(31)
binders = [e for e in test_supported if e.binder]
nonbinders = [e for e in test_supported if not e.binder]
nb = rng.choice(len(nonbinders), size=min(len(nonbinders), 6000 - len(binders)), replace=False)
subset = binders + [nonbinders[i] for i in nb]
print("subset", len(subset), flush=True)
peps = [e.sequence for e in subset]; als = [e.sequence and e.allele for e in subset]
mfl = []
B = 512
for s in range(0, len(subset), B):
    r = pred.predict(peptides=peps[s:s+B], alleles=als[s:s+B])
    mfl.append(-np.log(np.clip(np.asarray(r), 1e-6, None)))
mfl = np.concatenate(mfl)

# our ensemble on the same subset
import torch
from peptidehlacpp.models.cnn import PHLACNN
from peptidehlacpp.models.pssm import AllelePSSM
from peptidehlacpp.training.train_phla import enc_batch
alleles = sorted(keep); amap = {a: i for i, a in enumerate(alleles)}
by_al_train = defaultdict(list)
for e in train: by_al_train[e.allele].append(e)
pssms = {a: AllelePSSM().fit([e.sequence for e in exs], np.array([e.log_ic50 for e in exs])) for a, exs in by_al_train.items()}
cnn = PHLACNN(n_alleles=len(alleles)); cnn.load_state_dict(torch.load("results/phla_cnn.pt")["state"]); cnn.eval()
pssm_s = -np.array([pssms[e.allele].predict([e.sequence])[0] for e in subset])
cnn_s = []
with torch.no_grad():
    for s in range(0, len(subset), 2048):
        chunk = subset[s:s+2048]
        X, m, al = enc_batch([e.sequence for e in chunk], np.array([e.allele for e in chunk]), amap, False)
        cnn_s.append(cnn(X, m, al)[1].numpy())
cnn_s = np.concatenate(cnn_s)
def z(x): return (x - x.mean()) / (x.std() + 1e-9)
ens_s = z(pssm_s) + z(cnn_s)
labels = np.array([e.binder for e in subset])

# per-allele AUCs
per = {}
for a in sorted(set(als)):
    idx = [i for i, x in enumerate(als) if x == a]
    if len(idx) < 40: continue
    lab = labels[idx]
    if lab.sum() < 8 or (~lab.astype(bool)).sum() < 8: continue
    per[a] = {"n": len(idx), "n_binders": int(lab.sum()),
              "mhcflurry_auc": M.auc(lab, mfl[idx]),
              "pssm_auc": M.auc(lab, pssm_s[idx]),
              "cnn_auc": M.auc(lab, cnn_s[idx]),
              "ensemble_auc": M.auc(lab, ens_s[idx])}
wins = sum(1 for v in per.values() if v["ensemble_auc"] > v["mhcflurry_auc"])
res = {"n_alleles_compared": len(per), "ensemble_wins": wins, "per_allele": per}
json.dump(res, open("results/per_allele_headtohead.json", "w"), indent=2)
print(f"alleles compared: {len(per)}, ensemble wins: {wins}", flush=True)

fig, ax = plt.subplots(figsize=(5.4, 4.8))
xs = [v["mhcflurry_auc"] for v in per.values()]; ys = [v["ensemble_auc"] for v in per.values()]
ax.plot([0.5, 1], [0.5, 1], "k--", lw=1)
ax.scatter(xs, ys, s=28, c=["tab:red" if y > x else "tab:blue" for x, y in zip(xs, ys)])
for a, v in list(per.items()):
    if abs(v["ensemble_auc"] - v["mhcflurry_auc"]) > 0.06:
        ax.annotate(a.replace("HLA-", ""), (v["mhcflurry_auc"], v["ensemble_auc"]), fontsize=6)
ax.set_xlabel("MHCflurry 2.2.1 AUC (per allele)"); ax.set_ylabel("Our ensemble AUC (per allele)")
ax.set_title(f"Per-allele head-to-head on identical held-out pairs\n(ensemble wins {wins}/{len(per)} alleles)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figures/fig15_per_allele_headtohead.png", dpi=200)
print("DONE", flush=True)
