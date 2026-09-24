"""External-tool verification layer (steering: tools for research/verification).

1. sklearn.metrics.roc_auc_score cross-checks our Mann-Whitney AUC on the
   saved ensemble + head-to-head predictions (independent implementation).
2. Biopython PairwiseAligner re-verifies 9B-CPP novelty: max global identity
   vs CPPsite natural (dedup) and a UniProt pool, independent of our k-mer
   proxy (top-60 k-mer-nearest DB seqs aligned per candidate - the identity
   maximum is provably among high-Jaccard seqs for short peptides).
3. Bio.SeqUtils.ProtParam independently recomputes GRAVY + pI for the 18
   named candidates, cross-checking our hydropathy/charge descriptors.
Saves results/external_verification.json.
"""
import json
import sys
sys.path.insert(0, "src")
import numpy as np
from sklearn.metrics import roc_auc_score
from Bio.Align import PairwiseAligner
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from peptidehlacpp.data.cppsite import parse_fasta, natural_only, dedup_exact, kmer_jaccard
from peptidehlacpp.eval import metrics as M

out = {}

# 1. sklearn AUC cross-check
checks = []
for f, cols in [("ensemble_test_predictions", ["pssm", "cnn", "ensemble"]),
                ("h2h_scores", ["mhcflurry", "pssm", "cnn", "ensemble"])]:
    z = np.load(f"results/{f}.npz")
    y = z["labels"]
    for c in cols:
        mine = M.auc(y, z[c])
        sk = roc_auc_score(y, z[c])
        checks.append({"file": f, "column": c, "our_auc": float(mine),
                       "sklearn_auc": float(sk), "abs_diff": float(abs(mine - sk))})
out["auc_crosscheck"] = checks
out["auc_max_abs_diff"] = max(c["abs_diff"] for c in checks)

# 2. Biopython novelty re-verification
named = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
db = [p.sequence for p in dedup_exact(natural_only(parse_fasta("data/raw/cppsite2_natural.fa")))]
up = [p.sequence for p in natural_only(parse_fasta("data/raw/uniprot_reviewed_len8_35.fasta"))][:14400]
aln = PairwiseAligner()
aln.mode = "global"; aln.match_score = 1.0; aln.mismatch_score = 0.0
aln.open_gap_score = 0.0; aln.extend_gap_score = 0.0
def max_identity(query, pool):
    near = sorted(pool, key=lambda s: -kmer_jaccard(query, s))[:60]
    best = 0.0
    for s in near:
        sc = aln.score(query, s)
        ident = sc / max(len(query), len(s))
        best = max(best, ident)
    return float(best)
nov = []
for c in named:
    s = c["sequence"]
    nov.append({"name": c["name"],
                "bio_max_id_cppsite": max_identity(s, db),
                "bio_max_id_uniprot": max_identity(s, up),
                "kmer_proxy_max_jaccard": c["max_db_jaccard"]})
out["novelty_biopython"] = nov
out["novelty_max_id_overall"] = max(max(n["bio_max_id_cppsite"], n["bio_max_id_uniprot"]) for n in nov)

# 3. ProtParam descriptor cross-check
pp = []
for c in named:
    pa = ProteinAnalysis(c["sequence"])
    pp.append({"name": c["name"], "gravy_protparam": float(pa.gravy()),
               "our_mean_hydropathy": c["mean_hydropathy"],
               "pI": float(pa.isoelectric_point()),
               "charge_ph7_protparam": float(pa.charge_at_pH(7.0)),
               "our_net_charge": c["net_charge"]})
out["descriptor_crosscheck"] = pp
gd = [abs(p["gravy_protparam"] - p["our_mean_hydropathy"]) for p in pp]
cd = [abs(p["charge_ph7_protparam"] - p["our_net_charge"]) for p in pp]
out["gravy_max_abs_diff"] = float(max(gd))
out["charge_max_abs_diff"] = float(max(cd))
json.dump(out, open("results/external_verification.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}, indent=1))
