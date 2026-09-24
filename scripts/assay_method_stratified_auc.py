"""Ensemble AUC stratified by IEDB assay method.

Replicates the exact benchmark split (seed 7, peptide-level) to attach the
dominant assay method of each test (peptide, allele) pair, then recomputes
AUC within method strata from the saved ensemble scores. Verifies the
replicated split against the saved npz before trusting the join.
Output: results/assay_method_auc.json
"""
import json, sys
from collections import defaultdict, Counter
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
# dominant method per (peptide, allele)
meth = defaultdict(Counter)
for seq, al, m in df[["sequence", "allele", "method"]].itertuples(index=False):
    meth[(seq, al)][m] += 1
dom_meth = {k: v.most_common(1)[0][0] for k, v in meth.items()}

examples = aggregate(df)
counts = defaultdict(int)
for e in examples:
    counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

npz = np.load("results/ensemble_test_predictions.npz")
assert len(npz["labels"]) == len(test), (len(npz["labels"]), len(test))
lab_regen = np.array([e.binder for e in test])
al_regen = np.array([e.allele for e in test])
assert (lab_regen == npz["labels"]).all(), "label mismatch - split replication failed"
assert (al_regen == npz["alleles"]).all(), "allele mismatch - split replication failed"
print("split replication verified: labels and alleles match the saved npz")

labels, ens = npz["labels"], npz["ensemble"]
methods = np.array([dom_meth.get((e.sequence, e.allele), "unknown") for e in test])
top = [m for m, _ in Counter(methods).most_common(4)]
out = {"n_test": int(len(test)), "strata": []}
for m in top:
    sel = methods == m
    if sel.sum() < 50 or len(set(labels[sel])) < 2:
        continue
    out["strata"].append({
        "method": m, "n": int(sel.sum()),
        "binder_frac": float(labels[sel].mean()),
        "ensemble_auc": float(M.auc(labels[sel], ens[sel])),
    })
rest = ~np.isin(methods, top)
if rest.sum() > 0 and len(set(labels[rest])) > 1:
    out["strata"].append({
        "method": "all other", "n": int(rest.sum()),
        "binder_frac": float(labels[rest].mean()),
        "ensemble_auc": float(M.auc(labels[rest], ens[rest])),
    })
out["overall_auc"] = float(M.auc(labels, ens))
json.dump(out, open("results/assay_method_auc.json", "w"), indent=1)
print(json.dumps(out, indent=1))
