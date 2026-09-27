# PREREG B1 - Protocol implementation details (locked 2026-09-27, BEFORE any
# B1 training or scoring)

PREREG_B1_UNSEEN_ALLELE_2026-09-27.md locked the arms, endpoints, thresholds
and decision rules. It left the unseen-allele inference mechanics of each arm
unspecified. This note locks those mechanics BEFORE any B1 training or
scoring run. It changes nothing in the primary prereg; it only pins
implementation choices that were undetermined there.

## Training rows (all arms)
- ALL aggregated examples (geometric-mean IC50 per (peptide, allele), the
  project-standard aggregation in src/peptidehlacpp/data/iedb_dataset.py)
  whose allele is one of the 52 training alleles in
  results/per_allele_analysis.json. No held-out-allele row touches any fit,
  in any arm, in any form.
- CNN arm only: a 5% peptide-level validation split (seed 11) carved from
  the TRAINING rows, used solely for best-epoch selection (same role as the
  production val split). Val peptides never re-enter training.

## Arm 1 - production stack, unseen-allele mode (disclosed)
- The production ensemble is per-allele PSSM + allele-conditioned CNN. A
  held-out allele has zero training rows BY DESIGN, so the per-allele PSSM
  component cannot be fit for it. Arm 1 on unseen alleles therefore reduces
  to the CNN component alone. This is disclosed as a structural property,
  not an outcome.
- The CNN's allele embedding for a held-out allele is the MEAN of the 52
  learned training-allele embeddings (zero-shot mean-embedding). Locked
  before scoring.
- Hyperparameters mirror production: PHLACNN defaults, epochs 8, batch 512,
  lr 1e-3, AdamW wd 1e-4, seed 3, dual loss (MSE on log10 IC50 + BCE on
  binder), best-epoch-by-val-AUC state restored.
- Score = CNN classifier logit.

## Arm 2 - ESM-2 arm (disclosed)
- esm2_t12_35M_UR50D, frozen, mean-pooled final-layer embeddings (480-d).
  If OOM on this sandbox, the smallest variant (esm2_t6_8M_UR50D) is used
  and disclosed in results.
- Head: logistic regression (lbfgs, C=1.0, max_iter 2000) on training-row
  embeddings. Allele-agnostic (peptide embedding only) - disclosed: this arm
  tests how far a generic protein-LM peptide representation carries across
  alleles with no allele information at all.

## Arm 3 - PSSM additive baseline, unseen-allele mode (disclosed)
- A per-allele PSSM cannot be fit for an unseen allele (no rows). Arm 3 is
  therefore a POOLED allele-agnostic AllelePSSM (ridge 1.0) fit on all
  training rows. Score = -predicted log10(IC50). Locked before scoring.

## Scoring and metrics
- Each held-out allele scored ONCE per arm on its full usable row set
  (results/b1_allele_set.json); label = geo-mean IC50 <= 500 nM.
- Per allele: AUROC (primary), pAUC0.1, calibration slope (logistic fit of
  label on score within the allele; reported only when both classes have
  >=5 rows, else null with reason).
- An allele with a single observed class has no defined AUROC: excluded
  from the primary mean, disclosed by name, still included in the duplicate
  audit. (Locked handling, not a post-hoc exclusion.)
- Primary endpoint: mean per-allele AUROC, allele-equal weighting; 10,000
  resamples over alleles, seed 23 (bootstrap CI).

## Falsifier (locked target)
- Held-out allele HLA-A*32:07 (first allele in the frozen b1_allele_set.json
  order). Its arm-1 CNN scores are evaluated against a seed-29 label shuffle
  of that allele. AUROC must fall in [0.40, 0.60]; otherwise the harness is
  void until fixed.

## 9-mer near-duplicate audit (locked implementation)
- Training 9-mer pool: every contiguous 9-mer in every training peptide,
  plus every training 8-mer (for 8-mer queries).
- Per held-out peptide: hit = any of its contiguous 9-mers in the pool
  (len 8 peptides: the peptide itself in the 8-mer pool).
- Per allele: overlap = fraction of peptides with a hit. >0.90 flagged,
  not removed.
