"""AAindex1 physicochemical-axis analysis of the designed CPPs (tool: aaindex).

For every AAindex1 index with complete 20-AA values, scores each sequence by
its mean per-residue index value and measures separation (Cohen's d) between
the 18 designed CPPs and 300 seeded natural 8-35-mer UniProt decoys (same
pool/seed as the ESM-2 novelty screen). Reports the top separating axes and
whether the designed CPPs exaggerate them relative to the 7 CPPsite
archetypes. Output: results/cpp_aaindex_axes.json
"""
import json, random
import numpy as np
from aaindex import aaindex1

OUT = "results/cpp_aaindex_axes.json"
random.seed(7)

def load_cpp():
    cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
    return {c["name"] if isinstance(c, dict) else c[0]:
            c["sequence"] if isinstance(c, dict) else c[1] for c in cands}

def load_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line.replace("-", "")
    return seqs

cpps = load_cpp()
arch = {k: v for k, v in load_fasta("results/ebi_mafft_cpp_alignment.fa").items() if k not in cpps}
decoys = random.sample(list(load_fasta("data/raw/uniprot_reviewed_len8_35.fasta").values()), 300)

AA = "ARNDCQEGHILKMFPSTWYV"
def scorer(vals):
    v = np.array([vals[a] for a in AA])
    def f(seq):
        x = [v[AA.index(c)] for c in seq if c in AA]
        return float(np.mean(x)) if x else None
    return f

rows, skipped = [], 0
for code in aaindex1.record_codes():
    rec = aaindex1[code]
    vals = rec["values"]
    if any(vals.get(a) in (None, "-", []) for a in AA):
        skipped += 1
        continue
    f = scorer(vals)
    c = np.array([f(s) for s in cpps.values()])
    a = np.array([f(s) for s in arch.values()])
    d = np.array([x for x in (f(s) for s in decoys) if x is not None])
    sd = np.sqrt((c.var(ddof=1) + d.var(ddof=1)) / 2)
    cohen = (c.mean() - d.mean()) / sd if sd > 0 else 0.0
    rows.append({"accession": code, "description": rec["description"],
                 "cpp_mean": round(float(c.mean()), 3), "arch_mean": round(float(a.mean()), 3),
                 "decoy_mean": round(float(d.mean()), 3), "cohens_d": round(float(cohen), 3)})

rows.sort(key=lambda r: -abs(r["cohens_d"]))
out = {"n_indices_used": len(rows), "n_indices_skipped_incomplete": skipped,
       "n_cpps": len(cpps), "n_archetypes": len(arch), "n_decoys": len(decoys),
       "top_axes": rows[:20]}
json.dump(out, open(OUT, "w"), indent=1)
print(f"{len(rows)} indices used, {skipped} skipped")
for r in rows[:12]:
    print(f'{r["cohens_d"]:+.2f}  {r["accession"]}  {r["description"][:70]}  cpp={r["cpp_mean"]} arch={r["arch_mean"]} decoy={r["decoy_mean"]}')
