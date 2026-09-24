"""Independent cross-implementation verification of the CPP descriptors.

The screening cascade and the named-candidate table use features.py
(net charge, Eisenberg hydrophobic moment at 100 deg, Kyte-Doolittle mean
hydropathy). Here the same quantities are recomputed for all 18 named
candidates with two independent published packages - modlamp (charge,
pKa-based) and the peptides package (EMBOSS-scale charge, KD hydropathy,
hydrophobic moment) - and compared. Output: results/descriptor_crosscheck.json
"""
import json
import numpy as np
import peptides
import modlamp.descriptors as md

cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
ours_charge = np.array([c["net_charge"] for c in cands])
ours_moment = np.array([c["hydrophobic_moment"] for c in cands])
ours_gravy = np.array([c["mean_hydropathy"] for c in cands])

ml_charge, pt_charge, pt_gravy, pt_moment = [], [], [], []
for c in cands:
    s = c["sequence"]
    g = md.GlobalDescriptor([s])
    g.calculate_charge(ph=7.0, amide=False)
    ml_charge.append(float(g.descriptor[0][0]))
    p = peptides.Peptide(s)
    pt_charge.append(float(p.charge(pH=7.0, pKscale="EMBOSS")))
    pt_gravy.append(float(p.hydrophobicity(scale="KyteDoolittle")))
    w = min(11, len(s))
    pt_moment.append(float(p.hydrophobic_moment(angle=100, window=w)))

ml_charge = np.array(ml_charge); pt_charge = np.array(pt_charge)
pt_gravy = np.array(pt_gravy); pt_moment = np.array(pt_moment)

def pearson(a, b): return float(np.corrcoef(a, b)[0, 1])
def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)

out = {
    "n_candidates": len(cands),
    "mean_hydropathy_vs_peptides_KD": {
        "max_abs_diff": float(np.abs(ours_gravy - pt_gravy).max()),
        "pearson": pearson(ours_gravy, pt_gravy),
        "verdict": "identical scale, must match to rounding"},
    "net_charge_vs_modlamp_pka": {
        "max_abs_diff": float(np.abs(ours_charge - ml_charge).max()),
        "mean_signed_diff": float((ml_charge - ours_charge).mean()),
        "pearson": pearson(ours_charge, ml_charge),
        "spearman": spearman(ours_charge, ml_charge),
        "verdict": "integer sidechain count vs pKa-partials at pH 7: small "
                   "systematic offset expected (His ~0.1, termini), rank must agree"},
    "net_charge_vs_peptides_EMBOSS": {
        "max_abs_diff": float(np.abs(ours_charge - pt_charge).max()),
        "pearson": pearson(ours_charge, pt_charge),
        "spearman": spearman(ours_charge, pt_charge)},
    "hydrophobic_moment_vs_peptides_w11": {
        "pearson": pearson(ours_moment, pt_moment),
        "spearman": spearman(ours_moment, pt_moment),
        "verdict": "different windowing convention (full-length vs sliding "
                   "window <=11); correlation expected, equality not"},
}
json.dump(out, open("results/descriptor_crosscheck.json", "w"), indent=1)
print(json.dumps(out, indent=1))
