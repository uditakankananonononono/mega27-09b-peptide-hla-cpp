# Pre-registration addendum - 2026-09-26 (rules 1-8 revival), v2 CORRECTED
Locked BEFORE any new outcome is scored. Original registrations unchanged.
CORRECTION to the 5:02 PM version: its gap list was written before the full local
audit finished and mis-stated what exists. Verified against the live clone
(4a49f68 + this tree): the paired bootstrap ALREADY exists
(results/h2h_bootstrap.json: delta AUC +0.0117, 95% CI [0.0043, 0.0192],
p_win 1.0, 2000 replicates); the per-allele head-to-head ALREADY exists
(results/per_allele_headtohead.json: 30 alleles, ensemble wins 21); a
NetMHCpan-4.1-EL comparison via the official IEDB tools API ALREADY exists
(results/iedb_api_benchmark.json: IEDB AUC 0.8961 vs ensemble 0.9176 on the same
400-example subset); CellPPD-server candidate analysis exists
(results/cpp_cellppd_design.json); split is peptide-disjoint by construction
(split_by_peptide); the MHCflurry training-overlap caveat is documented in
scripts/head_to_head_mhcflurry.py (advantage to MHCflurry, reported honestly).

## TRUE remaining gaps (verified 5:05 PM)
- G1: the IEDB-API NetMHCpan subset (n=400) has no bootstrap CI - add one
  (10,000 replicates, paired) WITHOUT changing the subset or the scores.
- G2: per-allele record (21/30 wins) has no signed-rank test - add Wilcoxon over
  per-allele AUC deltas.
- G3: no held-out-ALLELE generalization experiment (train/test share alleles by
  design; split is peptide-disjoint). Locked protocol: allele-grouped split,
  retrain the FAST branches only (PSSM + k-mer logistic, locked), report per-allele
  mean AUC vs locally-run MHCflurry on identical allele-held-out test peptides.
- G4: CPP arm has no same-split benchmark vs a published CPP classifier - rebuild
  the CellPPD dataset split publicly and compare; or run a published tool's server
  where its license/endpoint allows, documented as in cpp_cellppd_design.json.
- G5: discovery naming - candidate set exists (cpp_novel_candidates.json +
  cpp_iedb_immunogenicity.json); needs a locked falsifiable discovery criterion and
  evidence cards per named candidate.
- G6: paper is a ~1.8k-word markdown draft, far below the 50+ text-page rule.
- G7: 0/10 judge rounds (rule 8: each must add a concrete novelty change).

## Beat gate (locked)
pHLA beat is considered ESTABLISHED only when G1+G2 pass (CIs excluding zero) and
G3 shows allele-held-out ensemble >= locally-run MHCflurry mean per-allele AUC on
the identical held-out alleles. CPP beat gate per G4. Any failed gate = documented
negative + rule-6 redirection.

## Judge rounds
Minimum 10, each producing a concrete novelty improvement (rule 8), verbatim in
JUDGE_ROUNDS.md.
