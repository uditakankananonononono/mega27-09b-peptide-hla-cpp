# PREREGISTRATION - G1 NEXTGEN RERUN (2026-09-28, locked before any scoring)

Parent approval: main-agent message 2026-09-28 00:09 IST ("queue the 09b G1
rerun after Comparator A lands, same discipline"). Trigger: the IEDB nextgen
tools API is reachable again (live-verified 2026-09-28), resolving the G1
403 egress blocker recorded in PREREG_ADDENDUM_2026-09-27.md.

## Scope (nothing else changes)
Rerun of two already-locked analyses, transport layer only:
1. G1 official-engine benchmark (originally scripts/iedb_api_benchmark.py,
   old tools-cluster-interface API): the SAME deterministic 400-example
   stratified subset of the held-out test split (rng seed 7, top-8 alleles,
   <=35 binders + <=15 non-binders per allele), rebuilt by the identical
   committed code path.
2. G1 paired bootstrap CI (locked in addendum v2, scripts/g1_iedb_bootstrap.py,
   never executed - blocked by the 403): 10,000 paired replicates, rng 123,
   identical resampling code.

## Locked mechanics
- Engine: IEDB nextgen API POST /api/v1/pipeline, tool_group mhci,
  predictors [{type: binding, method: netmhcpan_el}] (NetMHCpan EL
  percentile rank). Per (allele, exact length) batches; each peptide its
  own FASTA record with peptide_length_range [L,L] (no spurious k-mers).
- Raw API JSON cached and COMMITTED under results/g1_nextgen_raw/.
- Scores: netmhcpan_el_percentile per (peptide, allele); AUC computed on
  NEGATED percentile (higher = stronger binder), identical to the locked
  scripts; ensemble scores from the committed
  results/ensemble_test_predictions.npz (no refit, no new outcome data).
- Outputs: results/g1_iedb_nextgen_rerun.json (benchmark) and
  results/g1_bootstrap_nextgen.json (bootstrap CI).
- Decision rules and interpretation are unchanged from the original G1 and
  addendum v2; this rerun adds a cross-engine agreement note (old API
  0.8961 vs nextgen) purely as disclosure.

## Disclosures locked in advance
- Engine host/version changes (tools-cluster-interface -> nextgen); both
  serve NetMHCpan EL. Any numeric difference vs the original 400-example
  result is reported, not reconciled post-hoc.
- If the nextgen API drops or renames any allele/length group, the group is
  recorded as failed and excluded symmetrically from BOTH engines' AUCs
  (identical-example comparison), with the exclusion disclosed in the
  output JSON.
