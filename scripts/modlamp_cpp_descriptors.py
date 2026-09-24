"""modlamp descriptor panel: 18 designed CPPs vs 7 CPPsite archetypes.

Tool 28 (modlamp) deepened: full GlobalDescriptor panel (length, MW, charge,
charge density, pI, instability, aromaticity, aliphatic index, Boman index,
hydrophobic ratio) plus Eisenberg mean and windowed hydrophobic moment
(angle 100 = alpha helix) for the designed CPPs and the archetype set.
Output: results/cpp_modlamp_descriptors.json
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor

cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
cpps = {c["name"]: c["sequence"] if isinstance(c, dict) else None for c in cands}
cpps = {c["name"] if isinstance(c, dict) else c[0]:
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
arch = {k: v for k, v in load_fasta("results/ebi_mafft_cpp_alignment.fa").items() if k not in cpps}

def panel(seqs):
    g = GlobalDescriptor(seqs)
    feats = {}
    for name, fn in [("Length", g.length), ("MW", g.calculate_MW), ("Charge", g.calculate_charge),
                     ("ChargeDensity", g.charge_density), ("IsoelectricPoint", g.isoelectric_point),
                     ("InstabilityIndex", g.instability_index), ("Aromaticity", g.aromaticity),
                     ("AliphaticIndex", g.aliphatic_index), ("BomanIndex", g.boman_index),
                     ("HydrophobicRatio", g.hydrophobic_ratio)]:
        fn(); feats[name] = g.descriptor[:, -1].copy()
    e = PeptideDescriptor(seqs, "eisenberg"); e.calculate_global()
    feats["EisenbergMean"] = e.descriptor[:, -1].copy()
    e2 = PeptideDescriptor(seqs, "eisenberg"); e2.calculate_moment(window=11, angle=100)
    feats["EisenbergMoment11"] = e2.descriptor[:, -1].copy()
    return {s: {k: round(float(v[i]), 4) for k, v in feats.items()} for i, s in enumerate(seqs)}

cpp_rows = panel([cpps[n] for n in sorted(cpps)])
arch_rows = panel([arch[n] for n in sorted(arch)])
def rng(rows, k):
    v = [rows[s][k] for s in rows]
    return {"mean": round(float(np.mean(v)), 3), "min": round(min(v), 3), "max": round(max(v), 3)}
keys = list(next(iter(cpp_rows.values())).keys())
out = {"tool": "modlamp 4.3.0 (Muller et al. 2017) GlobalDescriptor + PeptideDescriptor(Eisenberg)",
       "summary": {k: {"designed_18": rng(cpp_rows, k), "archetypes_7": rng(arch_rows, k)} for k in keys},
       "by_name": {n: cpp_rows[cpps[n]] for n in sorted(cpps)},
       "archetypes": {n: arch_rows[arch[n]] for n in sorted(arch)}}
json.dump(out, open("results/cpp_modlamp_descriptors.json", "w"), indent=1)
print("wrote results/cpp_modlamp_descriptors.json:", len(cpp_rows), "designed,", len(arch_rows), "archetypes")
