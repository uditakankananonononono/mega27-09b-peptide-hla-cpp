"""G1 (locked addendum v2): paired bootstrap CI for the IEDB-API
(NetMHCpan 4.1 EL) vs ensemble comparison on the SAME deterministic
400-example subset as scripts/iedb_api_benchmark.py. Rebuilds the subset
identically, re-fetches IEDB percentile ranks (polite client, fresh cache),
loads committed per-example ensemble scores, 10,000-replicate paired
bootstrap of the AUC delta. No refit; no new outcome data."""
import json, os, sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M
from iedb_api_benchmark import iedb_predict

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples:
    counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)
rng = np.random.default_rng(7)
by_al = defaultdict(list)
for i, e in enumerate(test):
    by_al[e.allele].append(i)
top_alleles = sorted(by_al, key=lambda a: -len(by_al[a]))[:8]
sub_idx = []
for a in top_alleles:
    idx = np.array(by_al[a])
    binders = idx[[test[i].binder for i in idx]]
    nonb = idx[[not test[i].binder for i in idx]]
    nb = min(35, len(binders)); nn = min(15, len(nonb))
    sub_idx += list(rng.choice(binders, nb, replace=False))
    sub_idx += list(rng.choice(nonb, nn, replace=False))
sub_idx = sorted(set(sub_idx))
print(f"subset: {len(sub_idx)}", flush=True)

groups = defaultdict(list)
for i in sub_idx:
    e = test[i]
    groups[(e.allele, len(e.sequence))].append((i, e.sequence))
CACHE = "/tmp/iedb_g1_cache.json"
pred = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
pred = {int(k): v for k, v in pred.items()}
todo = [(k, v) for k, v in sorted(groups.items()) if not all(str(i) in pred for i, _ in v)]
for (a, L), items in todo:
    r = iedb_predict([s for _, s in items], a, L)
    for i, s in items:
        if s in r:
            pred[i] = r[s]
    json.dump({str(k): v for k, v in pred.items()}, open(CACHE, "w"))
    print(f"  {a} len{L}: cached", flush=True)

z = np.load("results/ensemble_test_predictions.npz", allow_pickle=True)
have = sorted(pred)
y = z["labels"][have]
iedb_rank = np.array([pred[i]["percentile_rank"] for i in have])
ens = z["ensemble"][have]
auc_iedb = M.auc(y, -iedb_rank)   # lower percentile rank = stronger binder
auc_ens = M.auc(y, ens)
print(f"AUC iedb {auc_iedb:.4f} ensemble {auc_ens:.4f} (n={len(y)})", flush=True)

rngb = np.random.default_rng(123)
n = len(y)
deltas = []
for _ in range(10000):
    b = rngb.integers(0, n, n)
    if y[b].sum() in (0, len(b)):
        continue
    deltas.append(M.auc(y[b], ens[b]) - M.auc(y[b], -iedb_rank[b]))
deltas = np.array(deltas)
out = {
    "n": int(n), "n_boot": int(len(deltas)),
    "iedb_auc": float(auc_iedb), "ensemble_auc": float(auc_ens),
    "delta_auc_mean": float(deltas.mean()),
    "delta_auc_ci95": [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))],
    "p_win": float((deltas > 0).mean()),
    "subset": "same deterministic 400-example subset as results/iedb_api_benchmark.json (rng 7)",
    "note": "IEDB percentile rank negated so higher = stronger binder; paired bootstrap on identical examples",
}
json.dump(out, open("results/g1_iedb_bootstrap.json", "w"), indent=1)
print(json.dumps(out, indent=1))
