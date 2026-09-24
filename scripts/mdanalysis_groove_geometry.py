"""Groove geometry of the 72-structure pHLA set with MDAnalysis.

For each peptide-HLA class I structure, identify the heavy (alpha) chain,
B2M and the peptide chain, then measure antigen-binding-groove geometry
with MDAnalysis atom selections and MDAnalysis.lib.distances:

  width_B_A  = |CA(63) - CA(163)|   alpha1-alpha2 helix opening over B pocket
  width_F_A  = |CA(81) - CA(146)|   helix opening over F pocket
  pep_span_A = |CA(P1) - CA(P-omega)|  peptide end-to-end span (bulge metric)
  depth_P2_A / depth_Pomega_A = signed distance of the anchor CA below its
    pocket rim plane (B-plane through CA 59,63,159,163; F-plane through
    CA 77,81,143,146); positive = buried toward the beta-sheet floor
    (sign fixed by floor residue CA 8)

Chain roles: the heavy chain is chosen among 200-320-CA chains by an
IMGT-numbering test (conserved Tyr at positions 7/59/84/159/171); this
rejects identically sequenced copies with author offset numbering
(e.g. 7k80/7k81 chain G, resids 6-294). Altloc duplicates take the first
CA. Single-chain trimers and complexes without a free 6-15-residue peptide
chain are recorded as honest skips.

MDAnalysis 2.9 has no mmCIF topology parser, so mmCIF -> PDB conversion is
done with Biopython (already in the tool inventory); all selections and
geometry are MDAnalysis.
"""
import json, os, glob, warnings, statistics as stats
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

OUT = "results/structure_mdanalysis_groove_geometry.json"
CONV = "/tmp/mda_conv"
os.makedirs(CONV, exist_ok=True)

_mf = json.load(open("data/raw/pdb_structures/search_manifest.json"))
ALLELE_OF = {p.lower(): a for a, ids in _mf.items() for p in ids}

TYR_TEST = (7, 59, 84, 159, 171)          # conserved IMGT tyrosines
B_PLANE = (59, 63, 159, 163)
F_PLANE = (77, 81, 143, 146)
FLOOR_RES = 8

def ca_pos(u, chain, resid):
    sel = u.select_atoms(f"chainID {chain} and resid {resid} and name CA")
    return sel.positions[0] if len(sel) else None

def pick_heavy(u):
    best, best_key = None, None
    for ch in np.unique(u.atoms.chainIDs):
        ca = u.select_atoms(f"chainID {ch} and name CA")
        if not (200 <= len(ca) <= 320):
            continue
        tyr = sum(1 for r in TYR_TEST
                  if u.select_atoms(f"chainID {ch} and resid {r} and name CA").resnames[:1].tolist() == ["TYR"])
        key = (tyr, -int(ca.resids.min()))
        if best is None or key > best_key:
            best, best_key = ch, key
    return best

def plane_fit(points):
    P = np.asarray(points)
    c = P.mean(axis=0)
    _, _, vh = np.linalg.svd(P - c)
    n = vh[-1]
    return c, n / np.linalg.norm(n)

parser = MMCIFParser(QUIET=True)
io = PDBIO()
rows, skipped = [], []
for fn in sorted(glob.glob("data/raw/pdb_structures/*.cif")):
    pdb = os.path.basename(fn)[:4]
    allele = ALLELE_OF.get(pdb)
    try:
        st = parser.get_structure(pdb, fn)
        pdb_path = os.path.join(CONV, pdb + ".pdb")
        io.set_structure(st)
        io.save(pdb_path)
        u = mda.Universe(pdb_path)
        heavy = pick_heavy(u)
        if heavy is None:
            skipped.append((pdb, "no 200-320 CA heavy-chain candidate"))
            continue
        pep_cands = []
        for ch in np.unique(u.atoms.chainIDs):
            if ch == heavy:
                continue
            n = len(u.select_atoms(f"chainID {ch} and name CA"))
            if 6 <= n <= 15:
                pep_cands.append(ch)
        if not pep_cands:
            skipped.append((pdb, "no free 6-15 residue peptide chain"))
            continue
        need = set(B_PLANE) | set(F_PLANE) | {63, 163, 81, 146, FLOOR_RES}
        lm = {r: ca_pos(u, heavy, r) for r in need}
        if any(v is None for v in lm.values()):
            skipped.append((pdb, "missing groove landmarks"))
            continue
        groove_center = (lm[63] + lm[163]) / 2.0
        def pep_dist(ch):
            cas = u.select_atoms(f"chainID {ch} and name CA").positions
            return float(distance_array(cas, groove_center.reshape(1, 3)).min())
        pep = min(pep_cands, key=pep_dist)
        if pep_dist(pep) > 15.0:
            skipped.append((pdb, "no peptide chain within 15 A of the groove"))
            continue
        d = lambda a, b: float(distance_array(lm[a].reshape(1, 3), lm[b].reshape(1, 3))[0, 0])
        width_B, width_F = d(63, 163), d(81, 146)
        pep_resids = sorted(set(u.select_atoms(f"chainID {pep} and name CA").resids))
        p1 = ca_pos(u, pep, pep_resids[0])
        p2 = ca_pos(u, pep, pep_resids[1])
        pw = ca_pos(u, pep, pep_resids[-1])
        span = float(distance_array(p1.reshape(1, 3), pw.reshape(1, 3))[0, 0])
        depths = {}
        for name, plane_res, anchor in (("depth_P2_A", B_PLANE, p2), ("depth_Pomega_A", F_PLANE, pw)):
            c, n = plane_fit([lm[r] for r in plane_res])
            if np.dot(n, lm[FLOOR_RES] - c) > 0:
                n = -n
            depths[name] = round(float(-np.dot(anchor - c, n)), 2)
        rows.append({
            "pdb": pdb, "allele": allele,
            "locus": allele.split("*")[0].replace("HLA-", "") if allele else None,
            "pep_len": len(pep_resids),
            "width_B_A": round(width_B, 2), "width_F_A": round(width_F, 2),
            "pep_span_A": round(span, 2), **depths,
        })
    except Exception as e:
        skipped.append((pdb, f"{type(e).__name__}: {e}"))

def summ(key, sub=None):
    v = [r[key] for r in rows if sub is None or r["locus"] == sub]
    return {"n": len(v), "mean": round(stats.mean(v), 2),
            "sd": round(stats.stdev(v), 2) if len(v) > 1 else 0.0,
            "min": round(min(v), 2), "max": round(max(v), 2)} if v else {}

summary = {
    "n_structures_analyzed": len(rows), "n_skipped": len(skipped),
    "skipped": skipped,
    "width_B_A": summ("width_B_A"), "width_F_A": summ("width_F_A"),
    "pep_span_A": summ("pep_span_A"),
    "depth_P2_A": summ("depth_P2_A"), "depth_Pomega_A": summ("depth_Pomega_A"),
    "by_locus": {l: {"n": sum(1 for r in rows if r["locus"] == l),
                     "width_B_A": summ("width_B_A", l),
                     "width_F_A": summ("width_F_A", l),
                     "depth_P2_A": summ("depth_P2_A", l),
                     "depth_Pomega_A": summ("depth_Pomega_A", l)}
                 for l in sorted({r["locus"] for r in rows if r["locus"]})},
}
try:
    from scipy.stats import spearmanr
    pairs = [("pep_len", "pep_span_A"), ("width_B_A", "depth_P2_A"),
             ("width_F_A", "depth_Pomega_A")]
    summary["spearman"] = {}
    for a, b in pairs:
        x = [r[a] for r in rows]; y = [r[b] for r in rows]
        sp = spearmanr(x, y)
        summary["spearman"][f"{a}_vs_{b}"] = {"rho": round(float(sp.statistic), 3),
                                              "p": float(sp.pvalue)}
except Exception as e:
    summary["spearman_error"] = str(e)

json.dump({"tool": "MDAnalysis " + mda.__version__, "rows": rows, "summary": summary},
          open(OUT, "w"), indent=1)
print(f"analyzed {len(rows)} structures, skipped {len(skipped)} -> {OUT}")
for k in ("width_B_A", "width_F_A", "pep_span_A", "depth_P2_A", "depth_Pomega_A"):
    print(k, summary[k])
print("spearman", json.dumps(summary.get("spearman"), indent=None))
print("loci", {l: v["n"] for l, v in summary["by_locus"].items()})
print("skipped", skipped)
