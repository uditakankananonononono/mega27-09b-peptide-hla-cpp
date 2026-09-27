# Detailed diagnostics for binding and CPP baselines

The following are transcriptions of committed result tables, not new model runs. A selected 6,000-row head-to-head panel and the 16,870-row peptide-held-out benchmark are different datasets; the 12-allele G3 transfer challenge is different again. Each denominator stays attached to its own metric.

## Within-panel comparison by allele

`results/per_allele_headtohead.json` contains 30 allele-specific comparisons on the selected same-row head-to-head panel. Overall 21 of 30 favor the ensemble numerically, but both allele prevalence and row counts vary. A per-allele numerical gain does not ensure a held-out-allele gain, as the G3 failure demonstrates.

| Allele | Rows | Binders | MHCflurry AUC | Ensemble AUC | Difference |
|:--|--:|--:|--:|--:|--:|
| HLA-A*01:01 | 136 | 102 | 0.9478 | 0.9230 | -0.0248 |
| HLA-A*02:01 | 736 | 660 | 0.9418 | 0.9420 | +0.0002 |
| HLA-A*02:02 | 229 | 218 | 0.9133 | 0.9333 | +0.0200 |
| HLA-A*02:03 | 302 | 276 | 0.9242 | 0.9331 | +0.0089 |
| HLA-A*02:06 | 262 | 241 | 0.8718 | 0.9176 | +0.0458 |
| HLA-A*03:01 | 372 | 314 | 0.8827 | 0.8731 | -0.0096 |
| HLA-A*11:01 | 341 | 296 | 0.9273 | 0.9101 | -0.0172 |
| HLA-A*23:01 | 74 | 63 | 0.9365 | 0.9380 | +0.0014 |
| HLA-A*24:02 | 128 | 104 | 0.8678 | 0.8702 | +0.0024 |
| HLA-A*24:03 | 54 | 45 | 1.0000 | 0.9753 | -0.0247 |
| HLA-A*26:01 | 89 | 59 | 0.8712 | 0.8842 | +0.0130 |
| HLA-A*29:02 | 96 | 81 | 0.9012 | 0.9267 | +0.0255 |
| HLA-A*30:01 | 99 | 88 | 0.8740 | 0.9349 | +0.0610 |
| HLA-A*31:01 | 279 | 236 | 0.9031 | 0.9225 | +0.0194 |
| HLA-A*33:01 | 142 | 114 | 0.9583 | 0.9558 | -0.0025 |
| HLA-A*68:01 | 254 | 229 | 0.8376 | 0.9027 | +0.0652 |
| HLA-A*68:02 | 226 | 186 | 0.9148 | 0.9113 | -0.0035 |
| HLA-A*69:01 | 63 | 39 | 0.9818 | 0.9562 | -0.0256 |
| HLA-B*07:02 | 196 | 171 | 0.9415 | 0.9331 | -0.0084 |
| HLA-B*08:01 | 126 | 103 | 0.8611 | 0.9164 | +0.0553 |
| HLA-B*15:01 | 195 | 166 | 0.8880 | 0.9094 | +0.0214 |
| HLA-B*15:17 | 69 | 55 | 0.9818 | 0.9909 | +0.0091 |
| HLA-B*18:01 | 53 | 27 | 0.9145 | 0.9359 | +0.0214 |
| HLA-B*27:05 | 85 | 58 | 0.9055 | 0.9623 | +0.0568 |
| HLA-B*39:01 | 44 | 33 | 0.8650 | 0.8788 | +0.0138 |
| HLA-B*40:01 | 92 | 74 | 0.9872 | 0.9910 | +0.0038 |
| HLA-B*44:02 | 68 | 56 | 0.9420 | 0.9896 | +0.0476 |
| HLA-B*51:01 | 70 | 47 | 0.8483 | 0.8085 | -0.0398 |
| HLA-B*57:01 | 76 | 58 | 0.9665 | 0.9933 | +0.0268 |
| HLA-B*58:01 | 122 | 99 | 0.9095 | 0.9368 | +0.0272 |

Comparisons on small or binder-heavy alleles can have few negatives, so the sampling variability of AUROC is uneven. The Wilcoxon result in the companion analysis summarizes the signs but is not a cluster-robust inference over independent cohorts. [PENDING: replicate an allele-held-out benchmark without selecting alleles or thresholds after these results.]

## Peptide-length and assay-composition sensitivity

`results/length_analysis.json` records the underlying peptide-length composition, heavily concentrated at length nine. `results/per_length_headtohead.json` narrows to the selected 6,000-row comparison and reports AUROC differences with bootstrap intervals. The two denominators must not be mixed.

| Length | Head-to-head rows | Binders | Ensemble AUC | MHCflurry AUC | Delta | 95% interval |
|--:|--:|--:|--:|--:|--:|:--|
| 10 | 1600 | 1410 | 0.8928 | 0.8703 | +0.0224 | [+0.0027, +0.0434] |
| 11 | 36 | 30 | 0.8500 | 0.9056 | -0.0580 | [-0.2857, +0.1010] |
| 8 | 28 | 23 | 0.8783 | 0.5826 | +0.2977 | [+0.0267, +0.5925] |
| 9 | 4330 | 3588 | 0.9392 | 0.9345 | +0.0046 | [-0.0029, +0.0120] |

The 8-mer stratum has only 28 rows and 11-mer has 36; no mechanistic length preference should be asserted from either. Length nine has 4,330 rows and a delta interval crossing zero. Length ten has 1,600 rows and a positive within-panel delta, but it is a post-hoc stratum and does not reverse the G3 transfer failure.

## CPP learning curve and compositional feature audit

The committed `results/cpp_learning_curve.json` fits a three-mer logistic model at several training fractions and scores the same 183-row fixed test. Four recorded seeds are used below; at the full fraction they produce identical AUROC. This is a data-efficiency diagnostic for this one split, not independent test-set replication.

| Training fraction | Training sequences | Mean AUROC | Across-seed SD |
|--:|--:|--:|--:|
| 0.10 | 92 | 0.8092 | 0.0289 |
| 0.25 | 230 | 0.8655 | 0.0070 |
| 0.50 | 461 | 0.9030 | 0.0059 |
| 0.75 | 692 | 0.9177 | 0.0040 |
| 1.00 | 922 | 0.9265 | 0.0000 |

The 0.10-to-1.00 training-size increase moves mean AUROC from 0.8092 to 0.9265; it cannot by itself prove that the CNN would catch up with enough data, because the learning curve here belongs to the logistic baseline only. `results/cpp_kmer_feature_weights.json` confirms exact refit-test AUROC 0.926523 and records 8,000 vocabulary features. Positive weights are strongest for short basic motifs; these are correlational model coefficients, not experimental membrane-transport mechanisms.

| Positive 3-mer | Logistic coefficient |
|:--|--:|
| `RRR` | +3.568 |
| `KKK` | +1.650 |
| `KKR` | +1.422 |
| `RKK` | +1.267 |
| `KRR` | +0.970 |
| `KLA` | +0.914 |
| `QRR` | +0.857 |
| `ARR` | +0.833 |
| `KRK` | +0.824 |
| `LAL` | +0.821 |
| `RRQ` | +0.794 |
| `RQR` | +0.789 |
| `LAK` | +0.765 |
| `LLK` | +0.726 |
| `RRA` | +0.681 |

Of the top 50 positive terms, 46 contain R or K, four contain W/Y/F, and 18 contain L/I/V according to the committed composition summary. Feature frequency and coefficient signs depend on the training set and background negatives; they may reflect source and amino-acid composition rather than true uptake. The separate CPP screen still needs measured intracellular delivery and hemolysis.

## Interpretation

These detailed tables explain two apparent contradictions without resolving them by rhetoric: the ensemble can win narrowly on a selected familiar-allele panel while failing badly on novel alleles; and a transparent k-mer baseline can beat a CNN on a small filtered CPP set while still being a poor proxy for measured cell entry. The next contribution must come from a new locked evaluation or wet-lab data, not extra pages of favorable post-hoc slices.
