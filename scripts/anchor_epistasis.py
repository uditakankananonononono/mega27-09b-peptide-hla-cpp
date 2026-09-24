"""Anchor-epistasis estimator (novel tool): per-allele pairwise interaction
contrast between peptide position 2 and the C-terminal anchor (P-Omega),
with bootstrap confidence intervals.

Method: two-way ANOVA decomposition of the P2 x P-Omega grid of mean
log10(IC50) on 9-mers. Interaction terms I[a,b] are the double-mutant-cycle
energies of Proposition 2; ||I||_2^2 (weighted) is the interaction variance
that bounds additive-model gains. Reported per allele with 200x peptide
bootstrap CIs. No public tool reports this quantity.
"""
import json
import sys
from collections import defaultdict
sys.path.insert(0, "src")

import numpy as np

from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv

HYDROPHOBIC = set("LIVMFWA")
MIN_CELLS = 8          # minimum filled cells for a valid grid
MIN_PER_CELL = 5       # minimum observations per (P2, POmega) cell
MIN_9MERS = 1500       # minimum 9-mer pairs per allele
N_BOOT = 200
RNG = np.random.default_rng(23)


def grid_stats(pairs):
    """pairs: list of (p2, pomega, log10ic50). Returns interaction matrix stats."""
    cells = defaultdict(list)
    for a, b, g in pairs:
        cells[(a, b)].append(g)
    cells = {k: np.mean(v) for k, v in cells.items() if len(v) >= MIN_PER_CELL}
    if len(cells) < MIN_CELLS:
        return None
    aas_p2 = sorted({k[0] for k in cells})
    aas_po = sorted({k[1] for k in cells})
    if len(aas_p2) < 3 or len(aas_po) < 3:
        return None
    idx_a = {a: i for i, a in enumerate(aas_p2)}
    idx_b = {b: i for i, b in enumerate(aas_po)}
    G = np.full((len(aas_p2), len(aas_po)), np.nan)
    for (a, b), m in cells.items():
        G[idx_a[a], idx_b[b]] = m
    # two-way decomposition on filled cells (row/col grand means over available)
    row_mean = np.nanmean(G, axis=1)
    col_mean = np.nanmean(G, axis=0)
    grand = np.nanmean(G)
    I = G - row_mean[:, None] - col_mean[None, :] + grand
    mask = ~np.isnan(I)
    l2 = float(np.sqrt(np.nansum(I ** 2) / mask.sum()))  # RMS interaction
    # hydrophobic-pair contrast: mean interaction where BOTH anchors hydrophobic
    hh = [I[idx_a[a], idx_b[b]] for a in HYDROPHOBIC for b in HYDROPHOBIC
          if a in idx_a and b in idx_b and not np.isnan(I[idx_a[a], idx_b[b]])]
    hh_val = float(np.mean(hh)) if hh else float("nan")
    return {"rms_interaction": l2, "hh_contrast": hh_val, "n_cells": int(mask.sum())}


def main():
    df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
    examples = aggregate(df)
    by_allele = defaultdict(list)
    for e in examples:
        if len(e.sequence) == 9:
            by_allele[e.allele].append((e.sequence[1], e.sequence[8], e.log_ic50))
    out = {}
    for allele, pairs in sorted(by_allele.items()):
        if len(pairs) < MIN_9MERS:
            continue
        base = grid_stats(pairs)
        if base is None:
            continue
        boots = []
        n = len(pairs)
        for _ in range(N_BOOT):
            idx = RNG.integers(0, n, n)
            s = grid_stats([pairs[i] for i in idx])
            if s:
                boots.append(s["rms_interaction"])
        ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) if boots else (None, None)
        out[allele] = {**base, "n_9mers": n,
                       "rms_interaction_ci95": [round(c, 4) if c else None for c in ci]}
        print(f"{allele}: rms_I={base['rms_interaction']:.4f} "
              f"CI95=[{ci[0]:.4f},{ci[1]:.4f}] hh={base['hh_contrast']:+.4f} "
              f"cells={base['n_cells']} n={n}", flush=True)
    json.dump(out, open("results/anchor_epistasis.json", "w"), indent=2)
    sig = [a for a, v in out.items() if v["rms_interaction_ci95"][0] and v["rms_interaction_ci95"][0] > 0.05]
    print(f"\nalleles with significant anchor epistasis (CI95 lower > 0.05 log10): {len(sig)}")
    print(sig[:20])


if __name__ == "__main__":
    main()
