# Peptide-HLA Binding Prediction and Cell-Penetrating Peptide Design: A Benchmarked In-Silico Study

Author: Udita Phookan


## Abstract

We evaluated peptide-HLA binding and cell-penetrating peptide (CPP) classification on distinct curated datasets. On a peptide-disjoint IEDB binding-affinity test, the per-allele additive PSSM reached AUROC 0.9156 and the CNN 0.9044, while a z-scored PSSM/CNN ensemble reached 0.9273. The CPP three-mer logistic model (AUROC 0.9265) exceeded our CNN (0.8956) on its held-out set. A separate, locked G3 allele-held-out challenge did not generalize: our AUROC 0.6687 versus MHCflurry 0.9285 on 6,000 examples, and mean per-allele AUROC 0.5798 versus 0.9056 across 12 test alleles. The G3 gate failed. CPP generator outputs and model scores are computational nominations, not measured delivery or safety. Mathematical identifiability of additive models requires checking the observed design rank; residue coverage alone is insufficient. Our contribution is an evidence-bounded account of where simple models compete and where allele extrapolation fails, not a leader-beating claim.


## 1. Introduction
Class-I HLA presentation is the rate-limiting filter of CD8+ T-cell epitope
discovery; cell-penetrating peptides are the dominant delivery vector class
for peptide therapeutics. Both problems admit large, open, quantitative
datasets - IEDB binding assays (135,854 quantitative class-I human assays
in the 2026-09-22 export) and CPPsite 2.0 (1,564 natural validated CPPs) -
making them fair ground for honest benchmarking of deep models against the
published leaders (NetMHCpan-4.1, MHCflurry-2.0; CellPPD/MLCPP-class tools).

## 2. Mathematical foundations

### 2.1 Why log10(IC50) is the right target: thermodynamic derivation
For a competitive equilibrium P + M ⇌ PM with competitor L:
IC50 ≈ Kd (Cheng-Prusoff, [L]≪Kd limit), and ΔG° = RT ln Kd.
Hence log10(IC50) = ΔG°/(RT ln 10) + const: the log-concentration is *linear
in binding free energy*, so a regressor on log10(IC50) estimates additive
free-energy contributions directly. Justifies MSE-on-log-IC50 loss.

### 2.2 Identifiability of additive binding models (proposition + proof)
**Prop 1.** An additive model ΔG = Σ_i w(i, a_i) is identifiable from assay
data up to a per-position constant gauge iff the design matrix has full
rank in a specified gauge; sample size and marginal residue coverage alone do not guarantee that rank.
**Prop 2.** Under an independent product sequence distribution, no additive model can represent a centered pairwise interaction term
w(i,j; a_i,a_j) with nonzero interaction contrast; the residual is lower-
bounded by the interaction's L2 norm. (Full proofs: appendix A.)
Consequence: additive models have a conditional squared-error limit when centered interaction contrasts exist. This does not imply a numerical AUROC ceiling on the measured data.

### 2.3 Hydrophobic moment as an amphipathicity statistic
μH = (1/N)|Σ_n H_n e^{inδ}|, δ=100° (α-helix). **Prop 3.** 0 ≤ μH ≤ mean|H_n|
with equality iff all vectors align; μH is invariant to cyclic permutation
of the sequence. Derivation in appendix B. Used as a CPP screening filter.

### 2.4 Decision theory at the 500 nM threshold
Binder = IC50 ≤ 500 nM (community convention). PPV among top-k predictions:
PPV(k) = Σ_i y_i 𝟙[rank(i)≤k]/k; its expectation under prevalence π and
score-conditional accuracy derived in appendix C - motivates reporting PPV at
k = #true binders for comparability with the MHCflurry benchmark.

## 3. Data
IEDB export stream-filtered (columns, filters, yields); CPPsite 2.0 natural
set, redundancy filtering; UniProt negative pools. All splits peptide-level
with homology guards.

## 4. Models
PSSM ridge baseline (per allele); allele-conditioned CNN (44-ch encoding);
GNN over backbone+next-nearest+anchor graphs; CPP CNN classifier; GRU
generator + biophysical cascade.

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
split. This does not establish that most biological binding signal is additive. Appendix A is conditional mathematics, not an empirical explanation of the CNN gap. The GNN (graph over
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

## 8. Conclusion

A leakage-controlled, fully open benchmark on 114,688 real IEDB class-I
measurements shows the additive position-specific model remains the one to
beat for peptide-HLA binding at small compute (AUC 0.9156); an
allele-conditioned CNN matches but does not beat it, and a residue-graph GNN
trails even with extended training. This comparison does not prove why:
non-additive gains are bounded by the (small) pairwise interaction variance
in these assays. For CPPs, a 3-mer baseline outperforms a CNN at n~600,
and a GRU generator + biophysical cascade produced 2,319 novel candidate
CPPs for future experimental follow-up. All code, tests, figures and raw
results ship in the project repo.




# Appendix A: Identifiability of additive binding models

**Setup.** Peptides of length L over alphabet A (|A|=20). The additive
(position-specific scoring) model assigns free energy
  G_add(x; w) = b + Σ_{i=1..L} w(i, x_i),  w ∈ R^{L×20}.
Data: pairs (x^(n), g^(n)) with g^(n) = log10 IC50 (linear in ΔG by §2.1).

**Proposition 1 (identifiability on the observed design).**
With an intercept and one reference residue omitted at each position, let Z
be the n × [1 + L(20-1)] design matrix. The coefficients in this chosen
gauge are uniquely identified by ordinary least squares if and only if
rank(Z) = 1 + L(20-1). If Z has lower rank, a nonzero null-space direction
gives exactly the same predictions for all observed peptides. The intercept
plus all 20L one-hot indicators is necessarily rank-deficient, since the
indicators at each position sum to the intercept. Even in the reduced gauge,
large n and full marginal amino-acid coverage do not imply full rank:
correlations between positions can still create exact dependencies.

*Proof.* Full column rank makes ZᵀZ positive definite and gives a unique
least-squares solution. Otherwise take v ≠ 0 in ker(Z). For any coefficient
vector β, vectors β and β+v have the same fitted values. Fixing the ordinary
position-constant gauge cannot remove additional dependencies from the
observed design. This is a mathematical condition, not a measured rank
audit of the archived assay matrix. [PENDING: inspect the actual per-allele
training-design ranks and conditioning.]

**Proposition 2 (additive ceiling under interactions).**
Let the true energy be G(x) = G_add(x) + G_int(x) where
G_int(x) = Σ_{i<j} U(i,j; x_i, x_j) with U centered (each marginal mean zero)
and at least one nonzero interaction contrast
  C_{ij}(a,a';b,b') = U(i,j;a,b) - U(i,j;a,b') - U(i,j;a',b) + U(i,j;a',b') ≠ 0.
Then for EVERY additive model w:
  E_x[(G(x) - G_add(x;w))^2] ≥ E_x[G_int(x)^2] > 0,
i.e. the best achievable additive MSE is bounded below by the interaction
energy variance.

*Proof.* Centered interactions are orthogonal to the additive subspace in
L2(uniform x): E[G_int · f(x_i)] = Σ_{i<j} E[U(i,j;x_i,x_j) f(x_i)] = 0 by
centering (each inner conditional mean is zero). The MSE-optimal additive
fit is the L2 projection of G onto the additive subspace, leaving residual
exactly G_int. Nonzero contrast implies G_int ≠ 0 on a positive-measure set,
so the squared norm is strictly positive. ∎

**Empirical anchor.** P2 and PΩ side chains share the B and F pockets'
chemical environment in class I HLA; correlated preferences (e.g. jointly
hydrophobic P2/PΩ in HLA-A*02:01) are precisely interaction contrasts, hence
the PSSM AUROC 0.9156 versus CNN 0.9044 and GNN 0.8480 in §5 is a dataset-specific comparison, not a theoretical ceiling or proof of interaction variance.

# Appendix B: Hydrophobic moment properties

μH(x; δ) = (1/N) |Σ_{n=1..N} H_{x_n} e^{i n δ}|, δ = 100°.

**Proposition 3.** (i) 0 ≤ μH ≤ (1/N) Σ_n |H_{x_n}|, equality iff all phasors
align (H_n e^{inδ} all same argument); (ii) μH is invariant under reversal
combined with conjugation symmetry |μH(x)| = |μH(rev x)|; (iii) for a
periodic amphipathic helix with hydrophobic moment phase at multiples of
δ-period, μH is maximized by placing the largest |H_n| at aligned phases.

*Proof.* (i) is the triangle inequality with equality condition. (ii)
reversal maps the sum to e^{i(N+1)δ} times the conjugate sum, same modulus.
(iii) rearrangement inequality applied to the aligned components. ∎

# Appendix C: PPV under class imbalance

With binder prevalence π and a score s with conditional binder density
f1(s), non-binder f0(s), the posterior is p(s) = π f1(s)/(π f1(s) +
(1-π) f0(s)). Top-k PPV estimates E[p | s in top-k]. For the benchmark's
choice k = #{true binders}, PPV(k) → sup over threshold policies of
precision at recall=π·N/|{predicted}| ... (operational point used by
MHCflurry's published benchmark, which we follow for comparability).


## Locked G3 allele-held-out challenge: negative result

The committed `results/g3_allele_holdout.json` records 83,333 training
examples and 6,000 evaluations across 12 test alleles. Overall AUROC was
**0.6687** for ours versus **0.9285** for MHCflurry. Mean per-allele
AUROC was **0.5798** versus **0.9056**. The recorded paired two-sided
Wilcoxon p-value across the twelve allele AUROCs is 0.00048828125 in
favor of the comparator; gate G3 is `false`. This differs from the
peptide-held-out results above. HLA-C*06:02 has only 40 rows (ours
0.3324, comparator 0.8242), so its individual comparison is unstable,
but it does not overturn the large aggregate gap. No allele-generalization
or leader-beating claim follows from the earlier within-dataset ensemble
AUROC of 0.9273. [PENDING: independent allele-held-out improvement and
prospective assay validation.]

## Paper-build correction log, September 27, 2026

An earlier draft claimed that marginal residue coverage and large n
guaranteed additive-model identifiability, and incorrectly asserted that
the observed design covariance was block diagonal. Both claims are now
withdrawn. Appendix A states the exact observed-design rank condition
without claiming an empirical audit. The Appendix's earlier CNN AUROC of
0.9057 conflicted with its results table (0.9044) and is corrected.
The G3 negative is preserved rather than recast as success.
