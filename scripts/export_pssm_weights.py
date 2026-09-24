"""Export per-allele PSSM ridge weights for the shipped CLI.

Refits AllelePSSM on the paper's exact training split (alleles with >=200
pairs, split_by_peptide defaults) and saves weights to results/pssm_weights.npz
so the peptidehlacpp CLI can reproduce the paper's ensemble from a fresh clone.
"""
import sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide
from peptidehlacpp.models.pssm import AllelePSSM
from peptidehlacpp.features import MAX_LEN_PHLA

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
examples = aggregate(df)
counts = defaultdict(int)
for e in examples: counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

by_al = defaultdict(list)
for e in train: by_al[e.allele].append(e)
alleles = sorted(by_al)
W = np.zeros((len(alleles), MAX_LEN_PHLA, 20), dtype=np.float64)
b = np.zeros(len(alleles), dtype=np.float64)
for i, a in enumerate(alleles):
    exs = by_al[a]
    m = AllelePSSM().fit([e.sequence for e in exs], np.array([e.log_ic50 for e in exs]))
    W[i] = m.W; b[i] = m.b
np.savez("results/pssm_weights.npz", alleles=np.array(alleles), W=W, b=b,
         max_len=MAX_LEN_PHLA, ridge=1.0)
print(f"saved results/pssm_weights.npz: {len(alleles)} alleles")

# smoke check against the live refit used by the paper pipeline
chk = AllelePSSM().fit([e.sequence for e in by_al["HLA-A*02:01"]],
                       np.array([e.log_ic50 for e in by_al["HLA-A*02:01"]]))
i = alleles.index("HLA-A*02:01")
assert np.allclose(W[i], chk.W) and np.isclose(b[i], chk.b)
print("smoke check ok")
