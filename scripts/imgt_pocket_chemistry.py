"""IMGT/HLA pocket chemistry across all covered alleles: parse hla_prot.fasta
(IPD-IMGT/HLA database flat file), take the reference (first) protein record
per covered 4-digit allele, and extract the B- and F-pocket residues in mature
protein numbering (signal peptide = 24 aa). Tests whether the pocket-identity
contrast found in three crystal structures (A*02: N66/Q70/D116; A*03/A*11:
K66/H70/Y116) holds across the full 51-allele panel, and computes the
epistasis-sign pocket table.
"""
import json

B_POCKET = [7, 9, 24, 34, 45, 63, 66, 67, 70, 99]
F_POCKET = [74, 77, 80, 81, 84, 95, 97, 114, 116, 123, 143, 146, 147]
LEADER = 24

# covered alleles from the benchmark per-allele analysis
cov = [r["allele"] for r in json.load(open("results/per_allele_analysis.json"))["per_allele"]]

seqs = {}
name, buf = None, []
for line in open("data/raw/hla_prot.fasta"):
    line = line.strip()
    if line.startswith(">"):
        if name:
            seqs[name] = "".join(buf)
        parts = line[1:].split()
        name = parts[1] if len(parts) > 1 else None  # e.g. A*02:01:01:01
        buf = []
    elif line:
        buf.append(line)
if name:
    seqs[name] = "".join(buf)

def pocket(seq, positions):
    out = []
    for p in positions:
        i = LEADER + p - 1
        out.append(seq[i] if i < len(seq) else "?")
    return "".join(out)

rows, missing = [], []
for allele in cov:
    short = allele.replace("HLA-", "")  # e.g. A*02:01
    keys = [k for k in seqs if ":".join(k.split(":")[:2]) == short]
    if not keys:
        missing.append(allele)
        continue
    key = sorted(keys)[0]  # reference (lowest-numbered) protein record
    s = seqs[key]
    if len(s) < LEADER + 147 or "X" in s[:LEADER + 147]:
        missing.append(allele + "(seq)")
        continue
    rows.append({"allele": allele, "imgt_record": key,
                 "b_pocket": pocket(s, B_POCKET), "f_pocket": pocket(s, F_POCKET),
                 "res66": s[LEADER + 65], "res70": s[LEADER + 69], "res116": s[LEADER + 115]})

json.dump({"n_alleles": len(rows), "missing": missing, "rows": rows},
          open("results/imgt_pocket_chemistry.json", "w"), indent=1)
print("covered:", len(rows), "missing:", missing)
from collections import Counter
print("res66:", Counter(r["res66"] for r in rows))
print("res70:", Counter(r["res70"] for r in rows))
print("res116:", Counter(r["res116"] for r in rows))
# link to epistasis sign: A*02 family favorable coupling vs A*03/A*11 unfavorable
for a in ["HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*15:01"]:
    r = next((x for x in rows if x["allele"] == a), None)
    if r:
        print(a, "B:", r["b_pocket"], "F:", r["f_pocket"], r["res66"], r["res70"], r["res116"])
