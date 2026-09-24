"""ProDy ANM normal-mode analysis of the MHC binding groove (1DUZ, HLA-A*02:01).

Question: are the anchor pockets (B pocket for P2, F pocket for P-Omega)
unusually rigid compared with the rest of the alpha chain? A rigidity
contrast would give a dynamics-level explanation for the anchor-burial
regularity measured with FreeSASA (paper anchor-burial section).
"""
import json
import numpy as _np
for _gone, _sub in (("alltrue", "all"), ("sometrue", "any"), ("product", "prod"),
                    ("cumproduct", "cumprod"), ("allstring", None)):
    if not hasattr(_np, _gone) and _sub is not None:
        setattr(_np, _gone, getattr(_np, _sub))
import sys as _sys, types as _types
if "numpy.lib.arraysetops" not in _sys.modules:
    try:
        import numpy.lib.arraysetops  # noqa
    except ModuleNotFoundError:
        _m = _types.ModuleType("numpy.lib.arraysetops")
        _m.isin = _np.isin
        _sys.modules["numpy.lib.arraysetops"] = _m
import scipy.linalg as _sla
_orig_eigh = _sla.eigh
def _eigh_compat(a, b=None, *args, **kw):
    ev = kw.pop("eigvals", None); kw.pop("turbo", None)
    if ev is not None:
        kw["subset_by_index"] = [ev[0], ev[1]]
    return _orig_eigh(a, b, *args, **kw)
_sla.eigh = _eigh_compat
from prody import parsePDB, ANM, calcSqFlucts, calcCollectivity

B_POCKET = [7, 9, 24, 34, 45, 63, 66, 67, 70, 99]
F_POCKET = [74, 77, 80, 81, 84, 95, 97, 114, 116, 123, 143, 146, 147]

st = parsePDB("data/raw/pdb/1duz.pdb", subset="calpha")
ca = st.select("chain A and resnum <= 181")
anm = ANM("1DUZ chain A alpha1-alpha2")
anm.buildHessian(ca)
anm.calcModes(n_modes=20)
sqf = calcSqFlucts(anm[:10])
resnums = ca.getResnums()
idx = {int(r): i for i, r in enumerate(resnums)}
mean_all = float(sqf.mean())
def pocket_ms(poss):
    vals = [float(sqf[idx[p]]) for p in poss if p in idx]
    return sum(vals) / len(vals), len(vals)
b_ms, bn = pocket_ms(B_POCKET)
f_ms, fn = pocket_ms(F_POCKET)
out = {
    "tool": "ProDy 2.4.1 (ANM, 20 modes, CA, cutoff 15 A default)",
    "structure": "1DUZ chain A (HLA-A*02:01 alpha1-alpha2, resnum 1-181)",
    "n_ca_atoms": int(ca.numAtoms()),
    "mean_sqfluct_all": round(mean_all, 4),
    "B_pocket_mean_sqfluct": round(b_ms, 4), "B_pocket_n": bn,
    "F_pocket_mean_sqfluct": round(f_ms, 4), "F_pocket_n": fn,
    "B_pocket_rigidity_ratio": round(mean_all / b_ms, 3) if b_ms else None,
    "F_pocket_rigidity_ratio": round(mean_all / f_ms, 3) if f_ms else None,
    "mode1_collectivity": round(float(calcCollectivity(anm[0])), 3),
}
json.dump(out, open("results/structure_prody_anm_1duz.json", "w"), indent=1)
print(out)
