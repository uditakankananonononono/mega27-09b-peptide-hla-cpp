"""Novelty verification for generated CPP candidates: exact + near-neighbor
search against CPPsite 2.0 (natural + non-natural) and the reviewed human
UniProt pool. Candidates with no near-neighbor (k-mer Jaccard >= 0.5 or
80% ungapped identity over any aligned window) are named 9B-CPP-N.
Also clusters the novel set (single-linkage at J>=0.5) to detect new
motif families absent from CPPsite."""
import json
import sys
sys.path.insert(0, "src")

from peptidehlacpp.data.cppsite import (dedup_exact, kmer_jaccard,
                                        natural_only, parse_fasta)

designs = json.load(open("results/cpp_designs.json"))["top_candidates"]
db_seqs = []
for f in ["data/raw/cppsite2_natural.fa", "data/raw/cppsite2_nonnatural.fa"]:
    db_seqs += [p.sequence for p in parse_fasta(f)]
for f in ["data/raw/uniprot_reviewed_len8_35.fasta",
          "data/raw/uniprot_human_reviewed_len60_400.fasta"]:
    db_seqs += [p.sequence for p in natural_only(parse_fasta(f))]
db_seqs = sorted(set(db_seqs))
print(f"database: {len(db_seqs)} sequences", flush=True)


def ungapped_identity(a: str, b: str) -> float:
    """Best ungapped identity of a against any same-length window of b."""
    if len(b) < len(a):
        a, b = b, a
    best = 0.0
    for i in range(len(b) - len(a) + 1):
        m = sum(1 for x, y in zip(a, b[i:i + len(a)]) if x == y) / len(a)
        if m > best:
            best = m
    return best


def nearest_hit(seq: str) -> tuple[float, float, str]:
    best_j, best_id, hit = 0.0, 0.0, ""
    for db in db_seqs:
        j = kmer_jaccard(seq, db)
        if j > best_j:
            best_j, hit = j, db
    best_id = max((ungapped_identity(seq, db) for db in db_seqs
                    if abs(len(db) - len(seq)) <= 3), default=0.0)
    return best_j, best_id, hit


novel = []
for c in designs:
    seq = c["sequence"]
    j, ident, hit = nearest_hit(seq)
    c["max_db_jaccard"] = round(j, 3)
    c["max_db_identity"] = round(ident, 3)
    if j < 0.5 and ident < 0.8:
        novel.append(c)
print(f"novel (J<0.5 and identity<0.8): {len(novel)} of {len(designs)}", flush=True)

# cluster novel set into motif families (single linkage at J>=0.5)
n = len(novel)
parent = list(range(n))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for i in range(n):
    for k in range(i + 1, n):
        if kmer_jaccard(novel[i]["sequence"], novel[k]["sequence"]) >= 0.5:
            pi, pk = find(i), find(k)
            if pi != pk:
                parent[pi] = pk
families = {}
for i in range(n):
    families.setdefault(find(i), []).append(novel[i]["sequence"])
families = sorted(families.values(), key=len, reverse=True)
print(f"motif families in novel set: {len(families)}; largest: {[len(f) for f in families[:5]]}")

named = []
for i, c in enumerate(novel, 1):
    c["name"] = f"9B-CPP-{i}"
    named.append(c)
json.dump({"n_screened": len(designs), "n_novel": len(novel),
           "named_candidates": named,
           "families": families[:10]},
          open("results/cpp_novel_candidates.json", "w"), indent=2)
for c in named[:12]:
    print(c["name"], c["sequence"], f"p={c['p_cpp']}", f"J={c['max_db_jaccard']}",
          f"id={c['max_db_identity']}")
