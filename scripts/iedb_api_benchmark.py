"""IEDB Analysis Resource tools API (MHC-I, method=recommended = NetMHCpan
4.1 EL): official-engine comparison and CPP immunogenicity screen.

Two analyses, one polite client (1 request per (allele, length) group,
<=100 sequences per request, 1.5 s between requests):

1. benchmark_subset: a seeded stratified 400-peptide subset of the held-out
   test split (rebuilt with the exact project data path), scored by IEDB;
   AUC of the IEDB percentile rank vs our ensemble on the same peptides.
2. cpp_immunogenicity: the designed CPPs (8-15-mers, MHC-I ligand range)
   against four common class-I alleles - a delivery-vehicle safety screen
   (low predicted binding = lower immunogenicity risk).

Usage: python3 scripts/iedb_api_benchmark.py
"""
import json, os, sys, time, subprocess
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M

API = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
OUT_B = "results/iedb_api_benchmark.json"
OUT_C = "results/cpp_iedb_immunogenicity.json"

def iedb_predict(seqs, allele, length):
    body = "\n".join(seqs)
    r = subprocess.run(
        ["curl", "-s", "-m", "120", "-X", "POST", API,
         "-d", f"method=recommended&allele={allele}&length={length}",
         "--data-urlencode", f"sequence_text={body}"],
        capture_output=True, text=True, timeout=130)
    if r.returncode != 0 or "peptide" not in r.stdout[:200]:
        raise RuntimeError(f"IEDB call failed: {r.stdout[:200]}")
    rows = {}
    for line in r.stdout.strip().split("\n")[1:]:
        f = line.split("\t")
        if len(f) >= 10:
            rows[f[5]] = {"score": float(f[8]), "percentile_rank": float(f[9])}
    return rows

# ---------- 1. benchmark subset ----------
print("rebuilding test split...", flush=True)
df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
from collections import defaultdict
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
print(f"subset: {len(sub_idx)} peptides, {len(top_alleles)} alleles", flush=True)

groups = defaultdict(list)
for i in sub_idx:
    e = test[i]
    groups[(e.allele, len(e.sequence))].append((i, e.sequence))
CACHE = "/tmp/iedb_bench_cache.json"
pred0 = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
pred = {int(k): v for k, v in pred0.items()}
todo = [(k, v) for k, v in sorted(groups.items()) if not all(str(i) in pred0 for i, _ in v)]
from concurrent.futures import ThreadPoolExecutor
def one(g):
    (a, L), items = g
    try:
        r = iedb_predict([s for _, s in items], a, L)
        got = {i: r[s] for i, s in items if s in r}
    except Exception as ex:
        print("group failed", a, L, ex, flush=True)
        got = {}
    return a, L, items, got
with ThreadPoolExecutor(max_workers=3) as ex:
    for a, L, items, got in ex.map(one, todo):
        pred.update(got)
        json.dump({str(k): v for k, v in pred.items()}, open(CACHE, "w"))
        print(f"  {a} len{L}: {sum(1 for i,_ in items if i in pred)}/{len(items)}", flush=True)

z = np.load("results/ensemble_test_predictions.npz", allow_pickle=True)
have = sorted(pred)
y = z["labels"][have]
iedb_rank = np.array([pred[i]["percentile_rank"] for i in have])
ens = z["ensemble"][have]
res = {
    "engine": "IEDB tools API mhci recommended (NetMHCpan 4.1 EL)",
    "n_requested": len(sub_idx), "n_scored": len(have),
    "alleles": top_alleles,
    "iedb_auc_on_rank": round(float(M.auc(y, -iedb_rank)), 4),
    "ensemble_auc_same_subset": round(float(M.auc(y, ens)), 4),
    "n_binders": int(y.sum()),
}
json.dump(res, open(OUT_B, "w"), indent=1)
print("BENCHMARK:", json.dumps(res, indent=1), flush=True)

# ---------- 2. CPP immunogenicity screen ----------
cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
cpps = {}
for c in cands:
    if isinstance(c, dict): cpps[c["name"]] = c["sequence"]
    else: cpps[c[0]] = c[1]
SCREEN_ALLELES = ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*35:01"]
cgroups = defaultdict(list)
skipped_len = []
for n, s in sorted(cpps.items()):
    if 8 <= len(s) <= 15:
        cgroups[len(s)].append((n, s))
    else:
        skipped_len.append((n, len(s)))
CCACHE = "/tmp/iedb_cpp_cache.json"
cres = json.load(open(CCACHE)) if os.path.exists(CCACHE) else {}
ctodo = []
for L, items in sorted(cgroups.items()):
    for a in SCREEN_ALLELES:
        if not all(n in cres and a in cres[n] for n, _ in items):
            ctodo.append((L, a, items))
def cone(g):
    L, a, items = g
    try:
        r = iedb_predict([s for _, s in items], a, L)
        return L, a, items, r
    except Exception as ex:
        print("cpp group failed", a, L, ex, flush=True)
        return L, a, items, {}
with ThreadPoolExecutor(max_workers=3) as ex:
    for L, a, items, r in ex.map(cone, ctodo):
        for n, s in items:
            cres.setdefault(n, {"sequence": s})[a] = r.get(s)
        json.dump(cres, open(CCACHE, "w"))
        print(f"  CPP len{L} {a} done", flush=True)
summary = {}
for n, d in cres.items():
    ranks = [v["percentile_rank"] for k, v in d.items() if k != "sequence" and v]
    summary[n] = {"sequence": d["sequence"],
                  "min_percentile_rank": min(ranks) if ranks else None,
                  "mean_percentile_rank": round(float(np.mean(ranks)), 3) if ranks else None,
                  "n_strong_binders_rank_lt0.5": sum(1 for r in ranks if r < 0.5),
                  "n_weak_binders_rank_lt2": sum(1 for r in ranks if r < 2.0)}
out = {"engine": "IEDB tools API mhci recommended (NetMHCpan 4.1 EL)",
       "alleles": SCREEN_ALLELES, "skipped_length_outside_8_15": skipped_len,
       "per_cpp": summary}
json.dump(out, open(OUT_C, "w"), indent=1)
print("CPP SCREEN:", json.dumps({n: (v['min_percentile_rank'], v['n_weak_binders_rank_lt2'])
                                  for n, v in summary.items()}, indent=1))
