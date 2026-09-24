"""Head-to-head vs MHCflurry within assay-method strata.

Replicates the exact 6,000-peptide h2h subset (rng seed 31 over the
MHCflurry-supported test peptides), verifies the replication against the
saved h2h score arrays, attaches each pair's dominant assay method, and
computes both models' AUC within each large stratum.
Output: results/assay_method_headtohead.json
"""
import json, sys
from collections import defaultdict, Counter
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
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

from mhcflurry import Class1AffinityPredictor
pred = Class1AffinityPredictor.load()
supported = set(pred.supported_alleles)
test_supported = [e for e in test if e.allele in supported]
rng = np.random.default_rng(31)
binders = [e for e in test_supported if e.binder]
nonbinders = [e for e in test_supported if not e.binder]
nb = rng.choice(len(nonbinders), size=min(len(nonbinders), 6000 - len(binders)), replace=False)
subset = binders + [nonbinders[i] for i in nb]

npz = np.load("results/h2h_scores.npz")
labels = np.array([e.binder for e in subset])
assert (labels == npz["labels"]).all(), "subset replication failed (labels)"
mfl_regen = []
for s in range(0, len(subset), 512):
    r = pred.predict(peptides=[e.sequence for e in subset[s:s+512]],
                     alleles=[e.allele for e in subset[s:s+512]])
    mfl_regen.append(-np.log(np.clip(np.asarray(r), 1e-6, None)))
mfl_regen = np.concatenate(mfl_regen)
max_diff = float(np.abs(mfl_regen - npz["mhcflurry"]).max())
assert max_diff < 1e-9, f"mhcflurry score mismatch {max_diff}"
print(f"subset replication verified (labels + MHCflurry scores, max diff {max_diff:.2e})")

ens = npz["ensemble"]
methods = np.array([dom_meth.get((e.sequence, e.allele), "unknown") for e in subset])
top = [m for m, _ in Counter(methods).most_common(3)]
out = {"n_subset": int(len(subset)), "replication_max_diff": max_diff, "strata": []}
for m in top:
    sel = methods == m
    if sel.sum() < 100 or len(set(labels[sel])) < 2:
        continue
    out["strata"].append({
        "method": m, "n": int(sel.sum()),
        "n_binders": int(labels[sel].sum()),
        "mhcflurry_auc": float(M.auc(labels[sel], npz["mhcflurry"][sel])),
        "ensemble_auc": float(M.auc(labels[sel], ens[sel])),
        "delta": float(M.auc(labels[sel], ens[sel]) - M.auc(labels[sel], npz["mhcflurry"][sel])),
    })
json.dump(out, open("results/assay_method_headtohead.json", "w"), indent=1)
print(json.dumps(out, indent=1))
