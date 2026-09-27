# Owner mega-verdict and agent paper-response queue

Provenance: original inbound WhatsApp message wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=, 2026-09-27 12:02:56 IST. The quoted blocks below are extracted without editing from the original channel body. The queue is agent-authored and is not a quotation or owner instruction.

Full original body SHA-256: d700f12c2a01d6f21b7392305aaa6b25d7a9efb29062ce5ba95c14e22f11bc33
Lane section SHA-256: 27cedc9916d3fa92f6e2b82bca44982bfe84fe7ce3086cc50d0db2a17b65a6df

## Agent-authored response queue, locked before manuscript edits

- P0: One question: when do deep models help versus additive baselines? G3 locked allele-held-out failure is the decisive negative and overrides any global leader claim.
- P1: Distinguish same-row G2 result from unseen-allele transfer, peptide-disjoint PSSM comparison, and CPP 3-mer LR loss. No external or prospective win inferred.
- P2: Explain plausible mechanisms of G3 failure as hypotheses only: allele representation coverage, shift in binding motif, assay composition, training diversity. Do not call any mechanism demonstrated without committed audits.
- P3: Add the decision-use and abstention path, clean count semantics, and pending predeclared unseen-allele/ESM/ligand-eluted/identifiability tests without numeric outcomes. Increase substantive pages.

## Owner header, verbatim

IGNORE ABOUT ISEF DELIVERABLES, IMPROVE PAGE COUNT

## Owner lane section, verbatim

8. manuscript.pdf: Peptide-HLA + CPP (09b)
Weaknesses (20) — computational only:

G3 allele-held-out challenge fails badly (0.6687 vs 0.9285).

Mean per-allele AUROC 0.5798 vs 0.9056.

CNN does not beat PSSM on peptide-disjoint split (0.9044 vs 0.9156).

GNN trails even with extended training (0.8480).

CPP CNN loses to 3-mer LR (0.8956 vs 0.9265).

G3 is a locked negative — overrides leader-beating claim.

Same-row MHCflurry margin is small (0.9281 vs 0.9164).

Per-allele G2 Wilcoxon p=0.0104 but alleles not independent.

Assay-method slices small (competitive fluorescence n=757).

Binder prevalence inflates PPV.

Calibration under new alleles not established.

pAUC0.1 = 0.52 — poor top-of-ranking.

CPP candidates Arg/Trp-rich — likely near-duplicates.

Max Jaccard 0.176-0.444 — not globally novel.

Identifiability propositions are conditional, not empirically audited.

40-tool gate open.

1,328 records are nested, not independent.

Paper thin (18 pages).

Allele-held-out failure not mechanistically explained.

No second allele-held-out benchmark.

Additions (computational):

Run a pre-registered new unseen-allele benchmark.

Test ESM-2 embeddings for allele generalization.

Add ligand-eluted benchmark (NetMHCpan task) for fair comparison.

Predeclare decision thresholds before new cohorts.

Add a proper identifiability rank audit on the observed design.

Simplify to one question: when do deep models help vs additive baselines.

Render 12-slide storyboard.

## Owner cross-cutting section, verbatim

Cross-Cutting Computational Themes
Recurring weaknesses:

Dataset/tool count inflation (nested records counted as independent).

Long papers (47-58 pages) — not ISEF-ready.

Negative-heavy narratives that obscure positive contributions.

Single-seed headline numbers.

Ad hoc gates/thresholds rather than theoretically derived.

Homology leakage in random-split benchmarks.

No leave-family-out CV in most projects.

CIs often overlapping — point-estimate wins only.

Universal computational additions:

One primary question per paper.

One locked primary endpoint.

Cluster-level bootstrap CIs everywhere.

Leave-family-out CV as the primary protocol.

12-slide storyboard as the ISEF deliverable.

One-page summary card.

Decision tree for tool use.

Pre-registered replication within the paper itself.

"What this is NOT" section in every abstract.

Reduce tool/dataset counting to study-level units.
