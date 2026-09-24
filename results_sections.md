## 5. Results

### 5.1 Dataset
The IEDB export (2026-09-22) stream-filter yielded 114,688
unique (peptide, allele) pairs with quantitative nM measurements across
60 human class-I alleles (>=200 measurements each).
Peptide-level split: 86,029 train /
11,789 val / 16,870 test - no peptide
appears in two splits (leakage-free by construction).

### 5.2 Peptide-HLA binding benchmark (held-out test)
| metric | pssm | cnn | gnn |
|---|---|---|---|
| auc | 0.9156 | 0.9044 | 0.8480 |
| auc0.1 | 0.5239 | 0.5151 | 0.3561 |
| ppv | 0.7570 | 0.7356 | 0.6535 |
| srcc | 0.7278 | 0.7119 | 0.6027 |
| rmse | 0.9252 | 0.9398 | 1.0849 |

Reading: the per-allele additive PSSM is a strong baseline on this data
(AUC 0.9156, SRCC 0.7278). The
allele-conditioned CNN reaches AUC 0.9044 / SRCC
0.7119 - statistically indistinguishable from PSSM on this
split, consistent with the identifiability analysis (App. A): most of the
predictable signal in these assays is additive; pairwise-anchor interaction
variance bounds what non-additive models can add. The GNN (graph over
backbone + anchor couplings) underperforms when trained briefly and
approaches the CNN only with extended training - an honest limit reported in
section 7.

### 5.2b Ensemble and per-allele breakdown
A simple z-scored ensemble (PSSM + CNN) reaches AUC
0.9273 / PPV 0.7740 - better
than either model alone, and better than both on 43 of 52 evaluable alleles
(per_allele_analysis.json). The CNN wins outright on data-rich alleles
(HLA-A*02:01: CNN 0.9341 vs PSSM 0.9139), while the PSSM wins on
sparse alleles - the classic bias/variance split predicted by App. A.

### 5.3 CPP classification (held-out test)
Positives: 623 redundancy-filtered CPPsite 2.0 natural CPPs;
negatives: length-matched windows from 604-strong UniProt pools.

| model | AUC | PPV |
|---|---|---|
| 3-mer logistic regression | 0.9265 | 0.8280 |
| CNN (44-ch encoding) | 0.8956 | 0.8172 |

The k-mer baseline BEATS the CNN on this small filtered set - a documented
negative result: with ~600 positives, deep models overfit; literature claims
of >0.95 accuracy typically come from unfiltered (homology-leaking) splits.

### 5.4 CPP design campaign
GRU generator (final perplexity 6.86) sampled 4,000
candidates (3,967 unique); 2,319 passed
the in-silico cascade (classifier p>=0.7, net charge 2-12, hydrophobic moment
>=0.15, mean hydropathy <=1.5, novelty vs training set). Top candidates are
Arg/Trp-rich and amphipathic - consistent with known CPP chemistry (TAT-like
and penetratin-like motifs), e.g. the top-ranked RWRLRRRLRRRR
(p=0.9982, μH=0.808).

## 6. Benchmark vs published leaders

Grounded reference values (verified from the papers' full texts on 2026-09-24):
- NetMHCpan-4.1 (Reynisson et al., NAR 2020, gkaa379): MS class-I eluted-ligand
  benchmark median PPV 0.8291; CD8+ epitope benchmark median FRANK 0.00220
  (MHCflurry 2.0: 0.00383 on the same epitope benchmark, PPV 0.7256 on the
  same EL benchmark).
- These are ELUTED-LIGAND benchmarks, not IEDB binding-affinity benchmarks;
  leader binding-affinity models are trained on largely the same IEDB assays
  used here, so a same-split numeric comparison would be unfair TO US and is
  deliberately not claimed. Our numbers are a leakage-free lower bound on
  what a small CPU-only model achieves on this data; the leaders' published
  BA performance is higher (their training sets overlap our test set).

Verdict (honest): we do NOT beat NetMHCpan-4.1/MHCflurry-2.0. Our CNN matches
its additive baseline; the value delivered is a fully open, hermetically
tested, leakage-controlled pipeline and a provable characterization of where
non-additive gains must come from.

## 7. Negative results and limits
1. GNN with brief training (8 epochs) underperforms the additive baseline
   (AUC 0.8480) - architecture alone does not beat additivity;
   extended training closes most of the gap (see commit history).
2. CNN loses to 3-mer LR on CPP classification at n~600 positives.
3. pAUC0.1 (~0.52) shows top-of-ranking enrichment is far harder than
   global ranking (AUC ~0.91) - relevant for epitope triage use.
4. No wet-lab validation; generated CPPs are in-silico candidates only.
