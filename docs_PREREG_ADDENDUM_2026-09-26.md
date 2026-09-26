# Pre-registration addendum - 2026-09-26 (rules 1-8 revival)
Locked BEFORE any new outcome is scored. Governs only new work; original registrations unchanged.

## B1 - pHLA benchmark-beat hardening (locked)
Current verified state (this agent, live clone 11f88be0, 33/33 tests green 4:23 PM):
head_to_head_mhcflurry.json on the shared 6,000-example test: ensemble AUC 0.9281 /
AUC0.1 0.6026 / PPV 0.9440 vs MHCflurry AUC 0.9164 / 0.5762 / 0.9383 - a numerical lead
on a matched, locally-run comparator, but NOT yet an established beat:
- Gap A: split cleanliness unaudited (allele disjointness; IEDB overlap with
  MHCflurry's own training set - overlap would disadvantage MHCflurry, so the audit
  can only strengthen or contextualize the claim, never inflate it).
- Gap B: no paired significance test (bootstrap on the 6,000 paired scores).
- Gap C: no held-out-ALLELE evaluation vs published NetMHCpan-4.1/MHCflurry-2.0
  reported numbers.
Locked protocol: A = leakage audit script (allele disjointness proof + IEDB-assay
overlap census vs MHCflurry training data); B = 10,000-replicate paired bootstrap on
AUC/AUC0.1/PPV deltas; C = allele-grouped evaluation (no allele in test appears in
train) scored per allele vs published per-allele numbers where available.
BEAT gate: bootstrap 95% CI of (ensemble - MHCflurry) AUC delta excludes 0 AND
allele-held-out ensemble mean AUC >= published comparator mean on the same alleles.
Any failed gate = documented negative + rule-6 redirection.

## B2 - CPP arm (locked)
Gap: no matched comparison vs a published CPP predictor. Protocol: rebuild the
CellPPD benchmark (public), score the committed kmer_lr + CNN, compare against
published CellPPD/MLCPP-class numbers on the identical split; BEAT gate = exceed the
strongest published comparator on the primary metric with bootstrap CI excluding 0.
Discovery arm: the committed cpp_novel_candidates.json set re-screened with the
improved classifier; novelty + toxicity flags per candidate.

## Judge rounds
Minimum 10, each producing a concrete novelty improvement (rule 8), verbatim in
JUDGE_ROUNDS.md.
