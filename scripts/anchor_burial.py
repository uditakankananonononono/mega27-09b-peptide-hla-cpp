"""Anchor burial in pHLA crystal structures, measured with FreeSASA.

For 1DUZ / 8RNI / 7OW3: compute per-residue solvent-accessible surface
area of the bound peptide, then rank P2 and P-omega against every other
peptide position. Prediction from the pocket analysis: the two groove
anchors are the MOST buried positions; the central bulge is the most
exposed. Output: results/anchor_burial.json
"""
import json
import freesasa
from Bio.PDB import PDBParser

THREE2ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
STRUCTS = {"1DUZ": ("data/raw/pdb/1duz.pdb", "C"),
           "8RNI": ("data/raw/pdb/8rni.pdb", None),   # peptide chain auto
           "7OW3": ("data/raw/pdb/7ow3.pdb", None)}

out = {}
for pdb_id, (path, pep_chain) in STRUCTS.items():
    st = PDBParser(QUIET=True).get_structure(pdb_id, path)[0]
    chains = {}
    for ch in st:
        res = [r for r in ch if r.id[0] == ' ' and r.resname in THREE2ONE]
        if res: chains[ch.id] = res
    if pep_chain is None:
        cands = [c for c, r in chains.items() if 8 <= len(r) <= 12]
        pep_chain = min(cands, key=lambda c: -len(chains[c]))
    pep = chains[pep_chain]
    pep_seq = ''.join(THREE2ONE[r.resname] for r in pep)

    fstruct = freesasa.Structure(path)
    result = freesasa.calc(fstruct)
    residue_areas = {}
    for i in range(result.nAtoms()):
        rn = fstruct.residueNumber(i); cn = fstruct.chainLabel(i)
        residue_areas[(cn, rn)] = residue_areas.get((cn, rn), 0.0) + result.atomArea(i)
    # map peptide residue order to SASA (residue numbers are strings)
    sasas = []
    for r in pep:
        rn = str(r.id[1]) + (r.id[2].strip() if r.id[2].strip() else '')
        sasas.append(residue_areas.get((pep_chain, str(r.id[1])), None))
    # freesasa residueNumber returns the pdb serial as string; robust fallback:
    if any(s is None for s in sasas):
        keys = [k for k in residue_areas if k[0] == pep_chain]
        keys.sort(key=lambda k: int(k[1]))
        sasas = [residue_areas[k] for k in keys[:len(pep)]]
    n = len(sasas)
    p2, pO = sasas[1], sasas[-1]
    order = sorted(sasas)
    rank = {v: i for i, v in enumerate(order)}
    out[pdb_id] = {
        "peptide_seq": pep_seq, "peptide_chain": pep_chain,
        "sasa_by_position": [round(s, 2) for s in sasas],
        "p2_sasa": round(p2, 2), "pomega_sasa": round(pO, 2),
        "p2_burial_rank": f"{n - rank[p2]}/{n}",
        "pomega_burial_rank": f"{n - rank[pO]}/{n}",
        "max_exposed_position": int(sasas.index(max(sasas))) + 1,
        "max_exposed_sasa": round(max(sasas), 2),
        "p2_in_top3_buried": bool(p2 <= order[2]),
        "pomega_in_top3_buried": bool(pO <= order[2]),
        "most_exposed_is_central": int(sasas.index(max(sasas))) + 1 in range(4, n - 2),
    }
json.dump(out, open("results/anchor_burial.json", "w"), indent=1)
print(json.dumps(out, indent=1))
