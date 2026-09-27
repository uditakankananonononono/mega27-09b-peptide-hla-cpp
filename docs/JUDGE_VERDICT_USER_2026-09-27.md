# USER-PROVIDED JUDGE VERDICT - 2026-09-27 (user mega-verdict section 8)

Source (archive of record): the user's own WhatsApp message 2026-09-27 12:02:56
IST, wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=,
verified author=user (16,557 body bytes, 10-project mega-verdict). This repo's
section and the Cross-Cutting Computational Themes were extracted DIRECTLY from
that authenticated message. Note: a courier-compiled file received the same
minute did NOT byte-match her message; it was rejected as archive source.
Header directive (verbatim first words of her message): "IGNORE ABOUT ISEF
DELIVERABLES, IMPROVE PAGE COUNT" - page-length weaknesses and storyboard
items are superseded; papers grow with substantive content. Under her rule
("EACH PROJECTS NEED ONE FROM ME TO PASS", WhatsApp 10:01:47 same day) this is
this project's ONE counted judge round. Factual claims are verified
independently before adoption.

---

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

---

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
