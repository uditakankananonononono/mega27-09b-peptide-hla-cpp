# PREREGISTRATION B1 - Unseen-allele benchmark (addendum 2026-09-27, B1)
Locked 2026-09-27 BEFORE any allele selection, training, or scoring. This is
the headline test of the locked G3 negative (allele-held-out AUROC 0.6687 vs
0.9285 in-distribution).

## Question
Does the model generalize to HLA alleles with NO training rows, when the
protocol, allele set, thresholds, and decision rules are fixed in advance?

## Data and allele set (selection rule, not cherry-picked)
- Source: data/processed/iedb_class1_human_nM.tsv (existing, sha256 on record).
- Held-out alleles: every allele with >= 50 usable rows (standard alphabet,
  length 8-14) that is NOT in the current training allele set, sorted by
  HLA supertype (G Sidney 2008 assignment), taking ALL such alleles - no
  performance-based exclusion. The exact list is computed once, committed as
  results/b1_allele_set.json with the selection code, and frozen.
- No re-selection after seeing results. Any exclusion discovered later
  (e.g. degenerate sequences) is disclosed with the reason, not silently
  dropped.

## Arms
1. Current production stack (the model family behind G1/G2).
2. ESM-2 embedding arm (B2): frozen ESM-2 (esm2_t12_35M_UR50D or the smallest
  available if memory-bound) mean-pooled embeddings + the same classifier
  head, same train rows. If the full ESM-2 model does not fit the sandbox,
  the smallest variant is used and disclosed.
3. PSSM additive baseline (the G1 comparator) - anchor for "when do deep
   models help".

## Protocol
- Train on ALL rows of the training allele set; score once per held-out
  allele. Single test evaluation per arm per allele.
- Primary endpoint: mean per-allele AUROC across held-out alleles
  (allele-equal weighting, locked).
- Secondary: per-allele AUROC distribution, pAUC0.1, calibration slope per
  allele.
- CIs: allele-level bootstrap (10,000 resamples over alleles), locked.

## Decision rules (predeclared thresholds, B4)
- Generalization SUCCESS: mean per-allele AUROC >= 0.80 AND bootstrap lower
  95% > 0.70 (the midpoint between chance and the in-distribution row).
- PARTIAL: mean >= 0.70 but lower95 <= 0.70.
- FAILURE: mean < 0.70.
- The G3 locked negative stands regardless; B1 measures whether the NEW
  protocol (and ESM-2 arm) repairs it. No arm result retroactively edits G3.
- Whichever arm is best is reported with the others alongside - no
  best-arm-only reporting.

## Falsifiers / honesty controls
- A label-shuffle arm on one held-out allele must score ~0.5 AUROC; if not,
  the harness is broken and all results are void until fixed.
- Near-duplicate audit: max 9-mer identity between each held-out allele's
  peptides and the training set reported per allele; alleles with >90%
  overlap are flagged, not removed.

## Reporting
results/b1_unseen_allele.json + a per-allele table in the paper. All three
arms reported as-is, including failures.
