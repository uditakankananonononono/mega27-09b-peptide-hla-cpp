"""Per-allele breakdown + PSSM/CNN ensemble evaluation on the held-out test."""
import json
import sys
from collections import defaultdict
sys.path.insert(0, "src")

import numpy as np
import torch

from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M
from peptidehlacpp.features import MAX_LEN_PHLA, backbone_graph, stacked_enc
from peptidehlacpp.models.cnn import PHLACNN
from peptidehlacpp.models.pssm import AllelePSSM
from peptidehlacpp.training.train_phla import MIN_PER_ALLELE, enc_batch

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples:
    counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= MIN_PER_ALLELE}
examples = [e for e in examples if e.allele in keep]
alleles = sorted(keep)
amap = {a: i for i, a in enumerate(alleles)}
train, val, test = split_by_peptide(examples)

by_al_train = defaultdict(list)
for e in train:
    by_al_train[e.allele].append(e)
pssms = {a: AllelePSSM().fit([e.sequence for e in exs],
                             np.array([e.log_ic50 for e in exs]))
         for a, exs in by_al_train.items()}

cnn = PHLACNN(n_alleles=len(alleles))
cnn.load_state_dict(torch.load("results/phla_cnn.pt")["state"])
cnn.eval()

def cnn_scores(exs):
    out = []
    with torch.no_grad():
        for s in range(0, len(exs), 2048):
            chunk = exs[s:s + 2048]
            X, m, al = enc_batch([e.sequence for e in chunk],
                                 np.array([e.allele for e in chunk]), amap, False)
            out.append(cnn(X, m, al)[1].numpy())
    return np.concatenate(out)

labels = np.array([e.binder for e in test])
reg_true = np.array([e.log_ic50 for e in test])
pssm_s = -np.array([pssms[e.allele].predict([e.sequence])[0] for e in test])
cnn_s = cnn_scores(test)

def z(x):
    return (x - x.mean()) / (x.std() + 1e-9)

ens_s = z(pssm_s) + z(cnn_s)
print("ENSEMBLE:", {"auc": round(M.auc(labels, ens_s), 4),
                    "ppv": round(M.ppv(labels, ens_s), 4)})

# per-allele AUCs for CNN vs PSSM, with train size
rows = []
by_al_test = defaultdict(list)
for i, e in enumerate(test):
    by_al_test[e.allele].append(i)
for a, idxs in by_al_test.items():
    idxs = np.array(idxs)
    if labels[idxs].sum() < 5 or (~labels[idxs]).sum() < 5:
        continue
    rows.append({"allele": a, "n_test": len(idxs),
                 "n_train": len(by_al_train[a]),
                 "pssm_auc": round(M.auc(labels[idxs], pssm_s[idxs]), 4),
                 "cnn_auc": round(M.auc(labels[idxs], cnn_s[idxs]), 4),
                 "ens_auc": round(M.auc(labels[idxs], ens_s[idxs]), 4)})
rows.sort(key=lambda r: -r["n_test"])
json.dump({"ensemble": {"auc": M.auc(labels, ens_s), "ppv": M.ppv(labels, ens_s)},
           "per_allele": rows}, open("results/per_allele_analysis.json", "w"), indent=2)
for r in rows[:12]:
    print(r)
wins = sum(1 for r in rows if r["cnn_auc"] > r["pssm_auc"])
print(f"CNN beats PSSM on {wins}/{len(rows)} alleles")
ens_wins = sum(1 for r in rows if r["ens_auc"] > max(r["pssm_auc"], r["cnn_auc"]))
print(f"Ensemble beats both on {ens_wins}/{len(rows)} alleles")
