# Evidence-locked companion analysis: binding and CPP screening

This companion expands the paper with committed outputs only. It does not rerun the science builder or turn retrospective comparisons into new experiments. Each section names its JSON source. The earlier within-allele and peptide-held-out results must be read alongside the independently failed G3 allele-held-out test.

## Matched-row binding benchmark: what actually happened

The committed `results/head_to_head_mhcflurry.json` evaluates 6,000 examples, 5,057 called binders under its defined threshold. High binder prevalence makes positive predictive value less revealing than AUROC or partial AUROC. The comparison uses predictions on the same recorded rows, not proof that each model had identical training access.

| Model | AUROC | AUROC at FPR ≤ 0.1 | PPV |
|:--|--:|--:|--:|
| MHCflurry | 0.9164 | 0.5762 | 0.9383 |
| Per-allele PSSM | 0.9202 | 0.5676 | 0.9425 |
| CNN | 0.9026 | 0.5127 | 0.9343 |
| PSSM/CNN ensemble | 0.9281 | 0.6026 | 0.9440 |

On these selected same rows, ensemble AUROC exceeds MHCflurry by 0.0118 in the stored point summaries. The committed `results/h2h_bootstrap.json` reports 2,000 resamples, a mean delta 0.01171 and a percentile 95% interval [0.00430, 0.01918]. This conditions on the sample and the evaluation protocol. It is not an independent held-out-allele replication, as G3 demonstrates. The bootstrap reports p_win=1 across the recorded resamples, which is not a literal probability of superiority in all populations.

The `results/g2_wilcoxon.json` per-allele comparison has 30 evaluable alleles and 21 ensemble wins; mean AUROC delta 0.0130, median 0.0110 and one-sided Wilcoxon p=0.0104. Alleles are not independent biological cohorts, so this test does not license a clinical-generalization claim. The sign and size of a small within-panel margin should be judged against the much larger held-out-allele failure.

## Locked G3 transfer failure, allele by allele

The separate `results/g3_allele_holdout.json` has 83,333 training examples, 6,000 evaluations and 12 held-out test alleles. These test alleles are not the same unit of generalization as a peptide-disjoint split within familiar alleles. G3 overall AUROC is 0.6687 versus 0.9285 for MHCflurry; mean per-allele AUROC 0.5798 versus 0.9056. Its recorded two-sided paired Wilcoxon p is 0.00048828125, favoring the comparator, and the gate is false.

| Held-out allele | n | Our AUROC | MHCflurry AUROC |
|:--|--:|--:|--:|
| HLA-A*02:01 | 2290 | 0.6794 | 0.9304 |
| HLA-A*02:03 | 933 | 0.6534 | 0.9142 |
| HLA-A*26:01 | 512 | 0.5553 | 0.9365 |
| HLA-A*32:01 | 154 | 0.6634 | 0.8301 |
| HLA-A*33:01 | 510 | 0.6557 | 0.8967 |
| HLA-A*69:01 | 365 | 0.6030 | 0.9496 |
| HLA-B*08:03 | 28 | 0.4470 | 0.8864 |
| HLA-B*15:03 | 113 | 0.5498 | 0.8879 |
| HLA-B*27:05 | 359 | 0.6007 | 0.9337 |
| HLA-B*39:01 | 234 | 0.5238 | 0.9330 |
| HLA-B*58:01 | 462 | 0.6937 | 0.9444 |
| HLA-C*06:02 | 40 | 0.3324 | 0.8242 |

All twelve recorded per-allele comparisons favor MHCflurry numerically. The 28-row HLA-B*08:03 and 40-row HLA-C*06:02 strata are noisy on their own, but the pooled and mean-per-allele gaps remain large. No data-dependent exclusion of those small alleles is justified after seeing the result. [PENDING: an architecture and training protocol locked before a new unseen-allele benchmark.]

## Assay-method slices and confounding limits

The committed `results/assay_method_headtohead.json` reports three assay-method strata. Their listed row counts do not sum to all 6,000 rows, so these are a subset view rather than an exhaustive decomposition; the file also records exact replication of its subset calculation (maximum difference zero).

| Method | n | Binders | MHCflurry | Ensemble | Delta |
|:--|--:|--:|--:|--:|--:|
| Competitive radioactivity | 4056 | 3677 | 0.8856 | 0.9002 | 0.0147 |
| Direct fluorescence | 1166 | 823 | 0.9539 | 0.9599 | 0.0060 |
| Competitive fluorescence | 757 | 548 | 0.9208 | 0.9230 | 0.0022 |

The three listed differences are positive but small. Multiple assay systems differ in measurement scale and chemistry. Any claim about a mutation mechanism requires assay-aware evaluation, not a post-hoc positive slice. `results/binder_confounds.json` fits a model to the 16,870-row held-out peptide panel controlling for length, assay method and locus; the ensemble coefficient is 1.3411 per standardized score with recorded odds ratio 3.823 and standard error 0.0214. This is association conditional on the recorded covariates, not causal proof or transfer validation. Its McFadden pseudo-R² changes from 0.0852 in the covariate-only model to 0.5008 when the score enters. Outcome-derived evaluation and feature construction must remain separated in any follow-up.

## Calibration and decision thresholds

The committed `results/calibration_decomposition.json` records n=16,870, base rate 0.2998 and binned Brier score 0.09839; a constant-base-rate predictor has 0.20991. Reliability term 0.000597 and resolution 0.11211 are binned decomposition quantities, not a guarantee of calibration under a new allele or clinic. Expected calibration error in `results/calibration.json` is 0.0164 for its bins. A threshold selected after observing a held-out cohort would invalidate a prospective PPV claim. [PENDING: predeclared decision threshold and calibration on new assay/allele domains.]

## Computational CPP nominations and safety boundary

The committed `results/cpp_novel_candidates.json` records 50 screened sequences, 18 meeting its implemented novelty filter, and ten sequence families. This is a filter result, not proof of biological novelty, delivery or safety. The highest-score sequences are often rich in Arg/Trp/Lys; without measured uptake and toxicity, classifier score alone is not an actionable treatment ranking.

| ID | Score | Charge | Moment | Max Jaccard |
|:--|--:|--:|--:|--:|
| 9B-CPP-1 | 0.9979 | 10.0 | 0.635 | 0.176 |
| 9B-CPP-2 | 0.9965 | 11.0 | 0.811 | 0.429 |
| 9B-CPP-3 | 0.9964 | 8.0 | 0.420 | 0.333 |
| 9B-CPP-4 | 0.9963 | 7.0 | 1.197 | 0.222 |
| 9B-CPP-5 | 0.9962 | 9.0 | 0.482 | 0.250 |
| 9B-CPP-6 | 0.9961 | 7.0 | 0.669 | 0.444 |
| 9B-CPP-7 | 0.9957 | 8.1 | 0.241 | 0.267 |
| 9B-CPP-8 | 0.9955 | 5.0 | 0.836 | 0.286 |
| 9B-CPP-9 | 0.9955 | 11.0 | 0.859 | 0.286 |
| 9B-CPP-10 | 0.9955 | 7.1 | 0.250 | 0.286 |
| 9B-CPP-11 | 0.9954 | 11.0 | 0.266 | 0.273 |
| 9B-CPP-12 | 0.9954 | 9.0 | 0.493 | 0.214 |
| 9B-CPP-13 | 0.9953 | 7.0 | 0.399 | 0.308 |
| 9B-CPP-14 | 0.9953 | 10.1 | 1.158 | 0.273 |
| 9B-CPP-15 | 0.9953 | 12.0 | 0.947 | 0.316 |
| 9B-CPP-16 | 0.9952 | 10.1 | 1.064 | 0.400 |
| 9B-CPP-17 | 0.9952 | 11.0 | 0.332 | 0.300 |
| 9B-CPP-18 | 0.9951 | 11.2 | 1.511 | 0.222 |

Sequences, in the same order as the table:

- 9B-CPP-1: ``LRILRWWRWRWRNRWKRRR``
- 9B-CPP-2: ``RLRRRLRLRLRLRLRLRRL``
- 9B-CPP-3: ``WRWRKKKRKK``
- 9B-CPP-4: ``GLWRWRWRWRKSLKK``
- 9B-CPP-5: ``RRINRRWRWKRWKKC``
- 9B-CPP-6: ``WRWRKKRRR``
- 9B-CPP-7: ``RGRWRWKKRHRPK``
- 9B-CPP-8: ``GLARWRWRWRRQ``
- 9B-CPP-9: ``RLRRRRRRNRRRNWRWL``
- 9B-CPP-10: ``WRWKWKKRWRKH``
- 9B-CPP-11: ``RLFRRRRLRRRRRRNN``
- 9B-CPP-12: ``RLRLRLRLRLRIRIRR``
- 9B-CPP-13: ``AWRWRWRCKAKKR``
- 9B-CPP-14: ``RRRRWRRLWRHLRRR``
- 9B-CPP-15: ``YLLRLRRRLRRRSRARRRR``
- 9B-CPP-16: ``RGRLRRHLRRRLRRR``
- 9B-CPP-17: ``RLRRRLRLRRFLRRRR``
- 9B-CPP-18: ``LGRLWRRLRRRRRHLRRAARSH``

Several candidates share basic-rich motifs and may be near-duplicates under a biological similarity metric even when a small reference search passes. The published `max_db_identity` entries for these candidates can be high; a pass against a defined reference database is not a global novelty claim. [PENDING: measured cellular uptake with endosomal-escape and membrane-damage controls, hemolysis testing, a broad homology audit and an assay-aware independent test.]

## Evidence map and limits

- `results/phla_benchmark.json`: 114,688 measured (peptide, allele) pairs and 60 alleles in the draft protocol. Its peptide-disjoint split supports within-panel generalization, not unseen-allele performance.
- `results/head_to_head_mhcflurry.json` and `results/h2h_bootstrap.json`: same-row comparison and paired empirical uncertainty under the selected test panel.
- `results/g2_wilcoxon.json`: per-allele within-panel comparisons. This is descriptive against the independent transfer failure.
- `results/g3_allele_holdout.json`: locked negative transfer test; it overrides any blanket leader-beating assertion.
- `results/cpp_benchmark.json` and `results/cpp_novel_candidates.json`: classifier and screen outputs, not wet-lab evidence.

These artifacts support a scoped scientific contribution: a transparent baseline-plus-ensemble panel, a failure mode exposed by allele-held-out evaluation, and testable computational CPP nominations. They do not establish a best-in-world immune prediction tool, peptide therapy, or clinical performance.
