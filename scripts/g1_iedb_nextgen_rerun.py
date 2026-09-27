"""G1 nextgen rerun (prereg docs/PREREG_G1_NEXTGEN_RERUN_2026-09-28.md).
Transport-only rerun of the locked G1 benchmark + addendum-v2 paired
bootstrap via the IEDB nextgen API. Subset, metrics, ensemble scores and
bootstrap resampling identical to the locked scripts. Resume-safe: raw
API JSON cached (and committed) under results/g1_nextgen_raw/.
Usage: python3 scripts/g1_iedb_nextgen_rerun.py
"""
import json, os, sys, time, urllib.request
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.eval import metrics as M

API = "https://api-nextgen-tools.iedb.org/api/v1"
RAW = "results/g1_nextgen_raw"

def post_json(url, payload, timeout=120):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())

def get_json(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.loads(r.read().decode())

def run_batch(allele, L, peps):
    fasta = "".join(f">p{i}\n{p}\n" for i, p in enumerate(peps))
    payload = {"pipeline_id": "", "run_stage_range": [1, 1], "stages": [
        {"stage_number": 1, "tool_group": "mhci",
         "input_sequence_text": fasta,
         "input_parameters": {"alleles": allele,
                              "peptide_length_range": [L, L],
                              "predictors": [{"type": "binding",
                                              "method": "netmhcpan_el"}]}}]}
    sub = post_json(f"{API}/pipeline", payload)
    if "results_uri" not in sub:
        raise RuntimeError(str(sub.get("errors", sub))[:300])
    uri = sub["results_uri"]
    for _ in range(90):
        time.sleep(4)
        res = get_json(uri)
        if res.get("status") == "done":
            return res
        if res.get("status") not in ("pending", "running", None):
            raise RuntimeError(f"bad status {res.get('status')}")
    raise TimeoutError(f"{allele} L{L} still pending")

def parse_percentile(res):
    out = {}
    for t in res.get("data", {}).get("results", []):
        if t.get("type") != "peptide_table":
            continue
        cols = [c["name"] for c in t["table_columns"]]
        ip, ipct = cols.index("peptide"), cols.index("netmhcpan_el_percentile")
        for row in t["table_data"]:
            out[row[ip]] = float(row[ipct])
    return out

# ---- rebuild the exact locked subset (verbatim code path, rng 7) ----
print("rebuilding test split + subset (rng 7)...", flush=True)
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
print(f"subset: {len(sub_idx)} peptides, {len(top_alleles)} alleles", flush=True)

# ---- fetch (resume-safe, committed raw cache) ----
groups = defaultdict(list)
for i in sub_idx:
    e = test[i]
    groups[(e.allele, len(e.sequence))].append((i, e.sequence))
os.makedirs(RAW, exist_ok=True)
failed = []
for (a, L), items in sorted(groups.items()):
    cache = f"{RAW}/{a.replace('*','').replace(':','')}_L{L}.json"
    if not os.path.exists(cache):
        peps = [s for _, s in items]
        try:
            res = run_batch(a, L, peps)
        except Exception as ex:
            for attempt in range(3):
                try:
                    time.sleep(10)
                    res = run_batch(a, L, peps)
                    break
                except Exception as ex2:
                    res = {"failed": True, "error": str(ex2)[:300]}
            if isinstance(res, dict) and res.get("failed"):
                pass
        json.dump(res, open(cache, "w"))
        time.sleep(2)
    res = json.load(open(cache))
    if res.get("failed"):
        failed.append((a, L))
        print(f"  {a} L{L}: FAILED (excluded symmetrically, disclosed)", flush=True)
    else:
        print(f"  {a} L{L}: cached", flush=True)

# ---- assemble per-example percentiles ----
pred = {}
for (a, L), items in sorted(groups.items()):
    cache = f"{RAW}/{a.replace('*','').replace(':','')}_L{L}.json"
    res = json.load(open(cache))
    if res.get("failed"):
        continue
    pct = parse_percentile(res)
    for i, s in items:
        if s in pct:
            pred[i] = pct[s]

z = np.load("results/ensemble_test_predictions.npz", allow_pickle=True)
have = sorted(pred)
y = z["labels"][have]
iedb_rank = np.array([pred[i] for i in have])
ens = z["ensemble"][have]
auc_iedb = float(M.auc(y, -iedb_rank))
auc_ens = float(M.auc(y, ens))
old = json.load(open("results/iedb_api_benchmark.json"))
bench = {
    "engine": "IEDB nextgen API mhci netmhcpan_el (NetMHCpan EL percentile)",
    "prereg": "docs/PREREG_G1_NEXTGEN_RERUN_2026-09-28.md",
    "n_requested": len(sub_idx), "n_scored": len(have),
    "failed_groups_excluded_symmetrically": [f"{a} L{L}" for a, L in failed],
    "alleles": top_alleles,
    "iedb_auc_on_rank": round(auc_iedb, 4),
    "ensemble_auc_same_subset": round(auc_ens, 4),
    "n_binders": int(y.sum()),
    "cross_engine_note": {
        "old_api_iedb_auc": old["iedb_auc_on_rank"],
        "old_api_ensemble_auc": old["ensemble_auc_same_subset"],
        "old_api_n_scored": old["n_scored"],
        "delta_iedb_auc_nextgen_minus_old": round(auc_iedb - old["iedb_auc_on_rank"], 4)},
}
json.dump(bench, open("results/g1_iedb_nextgen_rerun.json", "w"), indent=1)
print("BENCHMARK:", json.dumps(bench, indent=1), flush=True)

# ---- paired bootstrap (verbatim locked resampling, rng 123) ----
rngb = np.random.default_rng(123)
n = len(y)
deltas = []
for _ in range(10000):
    b = rngb.integers(0, n, n)
    if y[b].sum() in (0, len(b)):
        continue
    deltas.append(M.auc(y[b], ens[b]) - M.auc(y[b], -iedb_rank[b]))
deltas = np.array(deltas)
boot = {
    "engine": bench["engine"],
    "prereg": bench["prereg"],
    "n": int(n), "n_boot": int(len(deltas)),
    "iedb_auc": auc_iedb, "ensemble_auc": auc_ens,
    "delta_auc_mean": float(deltas.mean()),
    "delta_auc_ci95": [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))],
    "p_win": float((deltas > 0).mean()),
    "subset": "same deterministic 400-example subset as results/iedb_api_benchmark.json (rng 7)",
    "note": "IEDB percentile rank negated so higher = stronger binder; paired bootstrap on identical examples",
}
json.dump(boot, open("results/g1_bootstrap_nextgen.json", "w"), indent=1)
print("BOOTSTRAP:", json.dumps(boot, indent=1), flush=True)
print("G1 RERUN DONE", flush=True)
