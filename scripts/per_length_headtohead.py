"""Per-length head-to-head: ensemble vs MHCflurry, stratified by peptide length.

Controls for length composition as a confound of the overall head-to-head win
(the 10-mer AUC gap is documented in results/per_length_auc.json). Rebuilds the
exact h2h subset (verified recipe), confirms labels match results/h2h_scores.npz
bit-exactly, then computes per-length AUCs and a paired bootstrap on the delta.
"""
import json, sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
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

# --- replication guard: labels must match the committed score file exactly ---
saved = np.load("results/h2h_scores.npz")
labels = np.array([e.binder for e in subset])
assert labels.shape == saved["labels"].shape and (labels == saved["labels"]).all(), \
    "rebuilt subset does not match results/h2h_scores.npz"

lengths = np.array([len(e.sequence) for e in subset])
mfl = saved["mhcflurry"]; ens = saved["ensemble"]

def auc(y, s):
    return M.auc(y, s)

def paired_boot(y, s1, s2, n_boot=10000, seed=123):
    r = np.random.default_rng(seed)
    n = len(y)
    deltas = []
    idx = np.arange(n)
    for _ in range(n_boot):
        b = r.choice(idx, size=n, replace=True)
        if y[b].sum() in (0, len(b)):
            continue
        deltas.append(auc(y[b], s1[b]) - auc(y[b], s2[b]))
    deltas = np.array(deltas)
    return float(deltas.mean()), float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5)), float((deltas > 0).mean())

out = {"n_total": int(len(labels)), "replication": "labels match results/h2h_scores.npz bit-exactly",
       "strata": {}}
for L in sorted(set(lengths.tolist())):
    m = lengths == L
    y = labels[m]; e = ens[m]; f = mfl[m]
    if y.sum() == 0 or y.sum() == len(y):
        continue
    d, lo, hi, pwin = paired_boot(y, e, f)
    out["strata"][str(L)] = {
        "n": int(m.sum()), "n_binders": int(y.sum()),
        "auc_ensemble": round(float(auc(y, e)), 4),
        "auc_mhcflurry": round(float(auc(y, f)), 4),
        "delta": round(d, 4), "ci95": [round(lo, 4), round(hi, 4)],
        "p_win": round(pwin, 4),
    }
json.dump(out, open("results/per_length_headtohead.json", "w"), indent=1)
print(json.dumps(out, indent=1))
