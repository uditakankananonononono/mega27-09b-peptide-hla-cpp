"""Head-to-head: MHCflurry 2.x (published leader) vs our models on OUR
held-out test split - identical peptides, identical alleles, identical
labels. MHCflurry has train-overlap ADVANTAGE (it trained on many of these
IEDB assays); we evaluate on identical inputs anyway and report honestly."""
import json
import sys
from collections import defaultdict
sys.path.insert(0, "src")

import numpy as np

from peptidehlacpp.data.iedb_dataset import (aggregate, load_filtered_tsv,
                                             split_by_peptide)
from peptidehlacpp.eval import metrics as M

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples:
    counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

from mhcflurry import Class1AffinityPredictor  # noqa
pred = Class1AffinityPredictor.load()

# MHCflurry supports alleles like "HLA-A*02:01"; check coverage
supported = set(pred.supported_alleles)
test_supported = [e for e in test if e.allele in supported]
print(f"test pairs: {len(test)}; covered by MHCflurry: {len(test_supported)}", flush=True)

# subsample for speed: all binders + random negatives to 6000
rng = np.random.default_rng(31)
binders = [e for e in test_supported if e.binder]
nonbinders = [e for e in test_supported if not e.binder]
nb = rng.choice(len(nonbinders), size=min(len(nonbinders), 6000 - len(binders)), replace=False)
subset = binders + [nonbinders[i] for i in nb]
print(f"eval subset: {len(subset)} ({len(binders)} binders)", flush=True)

peps = [e.sequence for e in subset]
als = [e.allele for e in subset]
scores_aff = []
B = 512
for s in range(0, len(subset), B):
    r = pred.predict(peptides=peps[s:s+B], alleles=als[s:s+B])
    scores_aff.append(-np.log(np.clip(np.asarray(r), 1e-6, None)))  # higher = better binder
    print(f"  {min(s+B, len(subset))}/{len(subset)}", flush=True)
scores = np.concatenate(scores_aff)
labels = np.array([e.binder for e in subset])
res = {"n": len(subset), "n_binders": int(labels.sum()),
       "mhcflurry": {"auc": M.auc(labels, scores),
                     "auc0.1": M.auc_top_frac(labels, scores),
                     "ppv": M.ppv(labels, scores)}}
print("MHCFLURRY on our test subset:", {k: round(v, 4) for k, v in res["mhcflurry"].items()}, flush=True)

# --- our models on the identical subset ---
import torch
from peptidehlacpp.models.cnn import PHLACNN
from peptidehlacpp.models.pssm import AllelePSSM
from peptidehlacpp.training.train_phla import enc_batch

alleles = sorted(keep)
amap = {a: i for i, a in enumerate(alleles)}
by_al_train = defaultdict(list)
for e in train:
    by_al_train[e.allele].append(e)
pssms = {a: AllelePSSM().fit([e.sequence for e in exs],
                             np.array([e.log_ic50 for e in exs]))
         for a, exs in by_al_train.items()}
cnn = PHLACNN(n_alleles=len(alleles))
cnn.load_state_dict(torch.load("results/phla_cnn.pt")["state"])
cnn.eval()
pssm_s = -np.array([pssms[e.allele].predict([e.sequence])[0] for e in subset])
cnn_s = []
with torch.no_grad():
    for s in range(0, len(subset), 2048):
        chunk = subset[s:s+2048]
        X, m, al = enc_batch([e.sequence for e in chunk],
                             np.array([e.allele for e in chunk]), amap, False)
        cnn_s.append(cnn(X, m, al)[1].numpy())
cnn_s = np.concatenate(cnn_s)
def z(x): return (x - x.mean()) / (x.std() + 1e-9)
ens_s = z(pssm_s) + z(cnn_s)
for name, sc in [("pssm", pssm_s), ("cnn", cnn_s), ("ensemble", ens_s)]:
    res[name] = {"auc": M.auc(labels, sc), "auc0.1": M.auc_top_frac(labels, sc),
                 "ppv": M.ppv(labels, sc)}
    print(name.upper(), {k: round(v, 4) for k, v in res[name].items()}, flush=True)
json.dump(res, open("results/head_to_head_mhcflurry.json", "w"), indent=2)
