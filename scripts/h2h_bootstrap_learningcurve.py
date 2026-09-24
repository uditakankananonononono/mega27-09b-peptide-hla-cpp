"""Bootstrap CIs for the head-to-head metrics + PSSM learning curve."""
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

# --- rebuild the identical h2h subset and score with MHCflurry ---
from mhcflurry import Class1AffinityPredictor
pred = Class1AffinityPredictor.load()
supported = set(pred.supported_alleles)
test_supported = [e for e in test if e.allele in supported]
rng = np.random.default_rng(31)
binders = [e for e in test_supported if e.binder]
nonbinders = [e for e in test_supported if not e.binder]
nb = rng.choice(len(nonbinders), size=min(len(nonbinders), 6000 - len(binders)), replace=False)
subset = binders + [nonbinders[i] for i in nb]
labels = np.array([e.binder for e in subset])
mfl = []
for s in range(0, len(subset), 512):
    r = pred.predict(peptides=[e.sequence for e in subset[s:s+512]],
                     alleles=[e.allele for e in subset[s:s+512]])
    mfl.append(-np.log(np.clip(np.asarray(r), 1e-6, None)))
mfl = np.concatenate(mfl)
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
np.savez("results/h2h_scores.npz", labels=labels, mhcflurry=mfl, pssm=pssm_s, cnn=cnn_s, ensemble=ens_s)

# --- bootstrap CIs ---
rngb = np.random.default_rng(11)
B = 2000
n = len(labels)
boot = {"ensemble": [], "mhcflurry": []}
diffs = []
for b in range(B):
    idx = rngb.integers(0, n, n)
    if labels[idx].sum() in (0, len(idx)): continue
    a_e = M.auc(labels[idx], ens_s[idx]); a_m = M.auc(labels[idx], mfl[idx])
    boot["ensemble"].append(a_e); boot["mhcflurry"].append(a_m)
    diffs.append(a_e - a_m)
ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
res = {"n_boot": len(diffs),
       "ensemble_auc_ci95": ci(boot["ensemble"]),
       "mhcflurry_auc_ci95": ci(boot["mhcflurry"]),
       "delta_auc_mean": float(np.mean(diffs)),
       "delta_auc_ci95": ci(diffs),
       "p_win": float(np.mean(np.array(diffs) > 0))}
json.dump(res, open("results/h2h_bootstrap.json", "w"), indent=2)
print("bootstrap:", json.dumps(res, indent=1), flush=True)

# --- PSSM learning curve (fractions of train, 3 seeds) ---
fracs = [0.05, 0.1, 0.25, 0.5, 1.0]
curve = {}
test_labels = np.array([e.binder for e in test])
for fr in fracs:
    aucs = []
    for seed in range(3):
        r = np.random.default_rng(seed)
        tr_sub = [e for e in train if r.random() < fr]
        by_al = defaultdict(list)
        for e in tr_sub: by_al[e.allele].append(e)
        pm = {a: AllelePSSM().fit([e.sequence for e in exs], np.array([e.log_ic50 for e in exs]))
              for a, exs in by_al.items()}
        sc = -np.array([pm[e.allele].predict([e.sequence])[0] if e.allele in pm else 0.0 for e in test])
        aucs.append(M.auc(test_labels, sc))
    curve[str(fr)] = {"mean": float(np.mean(aucs)), "sd": float(np.std(aucs))}
    print(fr, curve[str(fr)], flush=True)
json.dump(curve, open("results/pssm_learning_curve.json", "w"), indent=2)
fig, ax = plt.subplots(figsize=(5.4, 4))
xs = [float(k) for k in fracs]
ys = [curve[str(fr)]["mean"] for fr in fracs]
es = [curve[str(fr)]["sd"] for fr in fracs]
ax.errorbar([f * len(train) for f in fracs], ys, yerr=es, fmt="o-", capsize=4)
ax.set_xscale("log")
ax.set_xlabel("training pairs (log scale)"); ax.set_ylabel("held-out test AUC")
ax.set_title("PSSM learning curve (3 seeds); CNN point at full data for reference")
ax.axhline(0.9044, color="tab:orange", ls="--", lw=1, label="CNN (full data)")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("paper/figures/fig18_learning_curve.png", dpi=200)
print("DONE", flush=True)
