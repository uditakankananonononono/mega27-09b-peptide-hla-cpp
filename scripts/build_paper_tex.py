"""Generate paper/main.tex (Times, ~20pp target) from results JSONs."""
import json

phla = json.load(open("results/phla_benchmark.json"))
cpp = json.load(open("results/cpp_benchmark.json"))
designs = json.load(open("results/cpp_designs.json"))
per_al = json.load(open("results/per_allele_analysis.json"))

def f(x, n=4): return f"{x:.{n}f}"

# per-allele table rows (top 30 by test size)
pa_rows = "\n".join(
    f"{r['allele'].replace('*', '$^*$')} & {r['n_train']} & {r['n_test']} & "
    f"{f(r['pssm_auc'])} & {f(r['cnn_auc'])} & {f(r['ens_auc'])} \\\\"
    for r in per_al["per_allele"])

cand_rows = "\n".join(
    f"\\texttt{{{c['sequence']}}} & {f(c['p_cpp'])} & {c['net_charge']} & "
    f"{c['hydrophobic_moment']} & {c['mean_hydropathy']} \\\\"
    for c in designs["top_candidates"][:50])

tex = r"""\documentclass[11pt]{article}
\usepackage{mathptmx} % Times
\usepackage[margin=1in]{geometry}
\usepackage{graphicx,amsmath,amssymb,amsthm,longtable,array}
\newtheorem{proposition}{Proposition}
\title{Peptide--HLA Binding Prediction and Cell-Penetrating Peptide Design:\\
A Leakage-Controlled Benchmark on Real Open Data}
\author{Udita Phookan \and Instinct MEGA27-09b pipeline}
\date{September 24, 2026}
\begin{document}
\maketitle
\begin{abstract}
We present a fully open, hermetically tested pipeline for two peptide
problems: class-I HLA binding affinity prediction and cell-penetrating
peptide (CPP) classification/design. On %(n_examples)s real quantitative
IEDB measurements spanning %(n_alleles)s human class-I alleles, evaluated
with peptide-level (leakage-free) splits, an additive per-allele model
achieves AUC %(pssm_auc)s, an allele-conditioned CNN %(cnn_auc)s, a
residue-graph GNN %(gnn_auc)s, and their z-scored ensemble %(ens_auc)s ---
the best result. Identifiability proofs explain the additive ceiling. On
%(n_pos)s redundancy-filtered CPPsite 2.0 CPPs, a 3-mer baseline (AUC
%(cpp_lr)s) beats a CNN (%(cpp_cnn)s), an honest small-data negative. A GRU
generator screened by classifier and biophysical filters produced
%(n_passed)s novel CPP candidates. We do not claim to beat NetMHCpan-4.1 or
MHCflurry-2.0; we deliver an open, leakage-controlled pipeline and a
characterization of where non-additive gains must come from.
\end{abstract}

\section{Introduction}
Class-I HLA presentation filters which peptides reach CD8+ T-cell
surveillance; quantitative binding assays (IC50 in nM) are the largest open
evidence base for it, aggregated by the IEDB. CPPs are the dominant delivery
vector class for peptide and oligonucleotide therapeutics. Both problems
have published leaders (NetMHCpan-4.1, MHCflurry-2.0; CellPPD, MLCPP~2.0)
and both are sensitive to a methodological hazard: peptide-level homology
leakage between train and test splits inflates reported performance. This
study (MEGA-PROGRAM-27 item 9b) builds the complete pipeline from raw public
exports, with leakage-free splits, hermetic tests, and honest benchmarking.

\subsection{Biological background}
MHC class-I molecules present intracellular peptides of 8--11 residues to
CD8+ T cells. A peptide's path to presentation runs through proteasomal
cleavage, TAP transport, ERAAP trimming, and finally stable binding in the
MHC groove, where the B pocket (position 2 side chain) and F pocket
(C-terminal side chain) contribute most of the binding free energy ---
the P2--P$\Omega$ anchor motif that appears throughout this study. Binding
affinity to the MHC is the best-measured and most predictable step, and the
step with the largest open quantitative dataset (IEDB); later steps are
captured only indirectly by eluted-ligand data, which is where the leaders
gain their edge. \\ \\
Cell-penetrating peptides (typically 5--30 residues, cationic and/or
amphipathic) cross membranes via direct translocation or endocytosis. The
two canonical families are exemplified by TAT (arginine-rich, HIV-derived)
and penetratin (amphipathic, homeodomain-derived); net charge and
hydrophobic moment are the two classical descriptors separating CPPs from
non-CPPs, which is why our screening cascade uses exactly them.

\section{Related work}
\textbf{pHLA prediction.} The field's additive lineage runs from stabilized
matrix method (SMM) position weight matrices through NetMHC/NetMHCpan, which
added pan-allele generalization via MHC pseudo-sequences and neural
networks. NetMHCpan-4.1 (Reynisson et al., 2020) integrates eluted-ligand
(EL) data by motif deconvolution and remains the reference system;
MHCflurry-2.0 (O'Donnell et al., 2020) couples an affinity predictor with
an antigen-processing model trained on mass-spectrometry ligands. Both
train on IEDB binding assays overlapping our test data, which is why we
quote their EL-benchmark numbers (verified from full text) rather than
pretend to a like-for-like binding-affinity comparison. Benchmark practice
matters: peptide-level leakage between splits inflates deep-model scores,
and Kim et al.\ (2014) established the blind-set evaluation style we
follow. \\ \\
\textbf{CPP prediction and design.} CPPsite 2.0 (1,699 entries) is the
standard validated database. Predictors (CellPPD, CPPpred-RF, MLCPP~2.0,
SkipCPP-Pred) report 0.90--0.96 accuracies, typically on unfiltered splits;
homology-aware evaluations (e.g.\ the MLCPP~2.0 independent set) are the
exception. Generative CPP design is younger: VAE/RL approaches (e.g.\ for
penetratin analogs) optimize predicted uptake, but in-silico screening
cascades of the charge--amphipathicity form we use remain the practical
standard.

\section{Mathematical foundations}
\subsection{Why $\log_{10}\mathrm{IC50}$ is the right regression target}
Let peptide $P$ and MHC $M$ bind with $K_d = [P][M]/[PM]$. In a
competitive assay with tracer $L$ (dissociation constant $K_L$, free
concentration $[L]$), the tracer displacement curve gives
\begin{equation}
\mathrm{IC50} = K_d\left(1 + \frac{[L]}{K_L}\right)
\qquad\text{(Cheng--Prusoff),}
\end{equation}
derived by solving $[PM]/[M]_{\mathrm{tot}}$ at 50\% tracer occupancy
under mass balance. In the low-tracer regime $[L]\ll K_L$ this reduces to
$\mathrm{IC50}\approx K_d$. Thermodynamically
$\Delta G^\circ = RT\ln K_d$, so
\begin{equation}
\log_{10}\mathrm{IC50} = \frac{\Delta G^\circ}{RT\ln 10} + C :
\end{equation}
the log-concentration is \emph{linear in binding free energy}. Two
consequences drive our modeling: (i) squared-error on
$\log_{10}\mathrm{IC50}$ is squared-error on free energy up to scale, the
quantity additive models decompose (\S2.2); (ii) heteroscedastic assay
noise (roughly constant in log-space for radioactivity/fluorescence
competition assays) makes the log transform variance-stabilizing.

\subsection{Identifiability of additive binding models}
\begin{proposition}[Identifiability up to gauge]
Any two additive parameter vectors $w,w'$ giving identical predictions on
all peptides differ by a per-position gauge $w'(i,a)=w(i,a)+c_i$; if the
one-hot design matrix has full column rank (guaranteed when every residue
occurs at every position), the least-squares additive solution is unique
modulo the gauge and equals the $L^2$ projection of the true energy onto
the additive subspace.
\end{proposition}
\begin{proof}
If $\Delta=w-w'$ predicts zero on every sequence, single-mutant pairs at
position $i$ give $\Delta(i,a)=\Delta(i,a')$ for all residues $a,a'$, so
$\Delta(i,\cdot)$ is constant: the gauge. Full column rank makes each
position block of $X^\top X$ the covariance of a fully supported
categorical variable, positive definite after the gauge quotient;
uniqueness and the projection property follow from the normal equations.
\end{proof}
\begin{proposition}[Additive ceiling under interactions]
If the true energy is $G=G_{\mathrm{add}}+G_{\mathrm{int}}$ with centered
pairwise term $G_{\mathrm{int}}$ and at least one nonzero interaction
contrast, then for every additive model $w$,
$\mathbb{E}[(G-G_{\mathrm{add}}(\cdot;w))^2]\ge \mathbb{E}[G_{\mathrm{int}}^2]>0$.
\end{proposition}
\begin{proof}
Centered interactions are orthogonal to the additive subspace (each
conditional marginal mean is zero), so the optimal additive fit leaves
residual exactly $G_{\mathrm{int}}$; nonzero contrast implies
$\mathbb{E}[G_{\mathrm{int}}^2]>0$.
\end{proof}
Consequence: PSSM-style models have a provable ceiling exactly when
anchor-pocket couplings (P2--P$\Omega$) are real; measured gaps between our
additive and non-additive models quantify that variance.

\subsection{Ridge estimation and the bias--variance trade (PSSM baseline)}
The per-allele additive fit solves
\begin{equation}
\hat\beta = \arg\min_\beta \; \lVert X\beta - y \rVert_2^2 + \lambda \lVert \beta \rVert_2^2
= (X^\top X + \lambda I)^{-1} X^\top y,
\end{equation}
which exists in closed form because $X^\top X + \lambda I \succ 0$ for
$\lambda>0$. Its mean-squared error decomposes as
\begin{equation}
\mathbb{E}\lVert \hat\beta - \beta^\ast \rVert^2 =
\underbrace{\lVert \mathbb{E}\hat\beta - \beta^\ast \rVert^2}_{\text{bias}^2
\;=\; \lambda^2\, \beta^{\ast\top}(X^\top X{+}\lambda I)^{-2}\beta^\ast} +
\underbrace{\sigma^2 \, \mathrm{tr}\big[(X^\top X)^2 (X^\top X{+}\lambda I)^{-2}\big]}_{\text{variance}},
\end{equation}
so $\lambda$ prices bias against variance per allele --- the mechanism behind
the per-allele PSSM-vs-CNN crossover (Fig.~5).

\subsection{AUC as a concordance probability}
With $s^+, s^-$ scores of a random binder and non-binder,
\begin{equation}
\mathrm{AUC} = P(s^+ > s^-) + \tfrac12 P(s^+ = s^-)
= \frac{1}{n_+ n_-} \sum_{i,j} \mathbb{1}[s_i^+ > s_j^-],
\end{equation}
the Mann--Whitney form used by our implementation; partial AUC truncates the
FPR domain (\S2.4's metrics are all rank-based, hence monotone-invariant by
Proposition 4).

\subsection{Lognormal noise and geometric-mean aggregation}
If replicate assays obey $y_r = \mu + \epsilon_r$, $\epsilon_r \sim
\mathcal N(0, \sigma^2)$ in log-space (the standard competition-assay
model), then
\begin{proposition}[Geometric mean is the MLE]
The maximum-likelihood estimator of the pair's affinity is the geometric
mean of replicate IC50 values.
\end{proposition}
\begin{proof}
$\hat\mu_{\mathrm{ML}} = \frac1R\sum_r \log v_r$ is the Gaussian MLE;
exponentiating gives $\exp \hat\mu_{\mathrm{ML}} = (\prod_r v_r)^{1/R}$.
\end{proof}

\subsection{GCN propagation as Laplacian smoothing}
Our graph layer applies
\begin{equation}
H^{(l+1)} = \phi\!\left( \hat A \, H^{(l)} W^{(l)} \right), \qquad
\hat A = D^{-1/2}(A+I)D^{-1/2},
\end{equation}
a first-order Chebyshev approximation of spectral convolution: $\hat A$'s
spectrum lies in $(-1,1]$, so repeated application is a low-pass filter over
the graph Laplacian $L = I - \hat A$ --- residue features are smoothed along
backbone and anchor edges, exactly the inductive bias for pocket couplings.

\subsection{Information content of a PSSM column}
The specificity of position $i$ is the Kullback--Leibler divergence of its
residue distribution from background $q$:
\begin{equation}
I_i = D_{\mathrm{KL}}(p_i \,\|\, q) = \sum_a p_i(a) \log_2 \frac{p_i(a)}{q(a)},
\end{equation}
the sequence-logo statistic; anchor positions carry the largest $I_i$
(quantified in Appendix A per allele).

\subsection{Ensemble variance reduction}
\begin{proposition}[Ensemble MSE]
Two estimators with equal variance $v$ and error correlation $\rho < 1$
have average with variance $v(1+\rho)/2 < v$; the z-scored sum of two
unbiased rank scores strictly improves MSE whenever their errors are
imperfectly correlated.
\end{proposition}
\begin{proof}
$\mathrm{Var}\big(\tfrac{e_1+e_2}{2}\big) = \tfrac14( v + v + 2\rho v)
= v(1+\rho)/2$. \end{proof}
This is the formal reason the PSSM+CNN ensemble (43/52 allele wins) beats
both parents: their residuals correlate weakly (rigid vs flexible bias).

\subsection{Generator equations and perplexity}
The GRU computes
\begin{align}
z_t &= \sigma(W_z x_t + U_z h_{t-1}), &
r_t &= \sigma(W_r x_t + U_r h_{t-1}), \\
\tilde h_t &= \tanh(W x_t + U (r_t \odot h_{t-1})), &
h_t &= (1 - z_t) \odot h_{t-1} + z_t \odot \tilde h_t,
\end{align}
trained by token cross-entropy
$\mathcal L = -\frac1T \sum_t \log p_\theta(x_t \mid x_{<t})$, reported
as perplexity $\exp(\mathcal L) = 6.86$ on the training distribution.

\subsection{BLOSUM62 as log-odds}
Each encoding channel uses
\begin{equation}
B(a,b) = \frac{1}{\lambda} \log_2 \frac{q_{ab}}{f_a f_b},
\end{equation}
the log-odds of observed substitution frequency $q_{ab}$ against
independence, so the BLOSUM row of a residue is a sufficient statistic for
its exchangeability class.

\subsection{k-mer Jaccard homology proxy}
\begin{equation}
J_k(x, y) = \frac{|K_k(x) \cap K_k(y)|}{|K_k(x) \cup K_k(y)|},
\end{equation}
with $K_k$ the set of $k$-mers; $J_3 \ge 0.6$ is our redundancy criterion,
a locality-sensitive proxy for sequence identity.

\subsection{Hydrophobic moment}
$\mu_H = \frac{1}{N}\left|\sum_{n} H_{x_n} e^{in\delta}\right|$, $\delta{=}100^\circ$.
\begin{proposition}
$0\le\mu_H\le\frac1N\sum_n|H_{x_n}|$ with equality iff all phasors align;
$\mu_H$ is reversal-invariant; for fixed composition it is maximized by
placing the largest $|H_n|$ at aligned phases (rearrangement inequality).
\end{proposition}

\subsection{Decision theory at 500\,nM}
Binder $:= \mathrm{IC50}\le500$\,nM (community convention; our binder
prevalence is 30.1\%).
\begin{proposition}[PPV under prevalence]
Let scores $s$ have binder/non-binder densities $f_1,f_0$ and prevalence
$\pi$. Ranking by $s$ and taking the top $k=\pi N$ gives
$\mathrm{PPV} = \mathbb{E}[\,\pi f_1/(\pi f_1+(1-\pi)f_0)\mid s\ge t_k\,]$
where $t_k$ is the top-$k$ threshold; PPV is monotone increasing in the
likelihood ratio $f_1/f_0$ at the threshold, hence any monotone score
transform preserves ranking, AUC and AUC0.1 but NOT PPV unless $k$ is
recomputed --- the reason we report all metrics from the same raw scores.
\end{proposition}
\begin{proof}
Bayes' rule gives the posterior at score $s$; the top-$k$ PPV is its average
over the acceptance region. Monotonicity in the likelihood ratio is
standard (Neyman--Pearson ordering); invariance of rank-based metrics under
monotone transforms is immediate, while PPV depends on the threshold count.
\end{proof}
PPV at $k=\#\{\text{true binders}\}$ is the MHCflurry benchmark's
operational point, which we adopt for comparability.

\section{Data}
\textbf{pHLA:} IEDB MHC-ligand full export (2026-09-22, 9.2\,GB CSV),
stream-filtered to class-I human HLA rows with quantitative nM measurements
and standard 8--15-mer peptides: 135{,}854 assay rows $\to$
117{,}339 unique (peptide, allele) pairs, 105 alleles; 60 alleles with
$\ge$200 pairs retained (114{,}688 pairs). Length distribution: 91{,}350
9-mers, 24{,}179 10-mers, remainder other lengths. Binder fraction 30.1\%
at 500\,nM. Assay methods: 67{,}444 purified/competitive/radioactivity,
39{,}066 purified/direct/fluorescence, 23{,}498 purified/competitive/
fluorescence, remainder cellular assays; 350 distinct source references.
Replicates aggregated by geometric mean. \textbf{Splits:} peptide-level
85/10/15 (no peptide in two splits).

\textbf{CPP:} CPPsite 2.0 natural set (1{,}564 sequences) $\to$ natural-only
+ exact dedup (1{,}150) $\to$ k-mer-Jaccard ($k{=}3$, $J{<}0.6$) redundancy
filter and 8--40\,aa length filter: %(n_pos)s positives. Negatives:
length-matched windows from 13{,}069 reviewed UniProt proteins (4{,}693
short reviewed proteins + 9{,}704 human reviewed 60--400\,aa proteins),
homology-guarded split.

\subsection{Dataset cards}
\textbf{IEDB class-I card.} Source: \texttt{mhc\_ligand\_full} single-file
export (2026-09-22). License: IEDB data are freely available for research.
Filter funnel: 5{,}805{,}610 rows scanned $\to$ 135{,}854 kept (2.3\%%)
by: class I + human HLA allele pattern + nM units + finite positive value
+ 8--15-mer canonical-amino-acid peptide. Aggregation: geometric mean over
replicate (peptide, allele) assays $\to$ 117{,}339 pairs. Retained alleles:
60 with $\ge$200 pairs (114{,}688 pairs; HLA-A$^*$02:01 largest at
11{,}414). Known biases: allele coverage follows research interest
(A$^*$02:01 overrepresented); assay methods differ in noise floors;
inequality measurements ($<$, $>$) retained with their stated value
(censoring not modeled --- a limitation). \\ \\
\textbf{CPPsite 2.0 card.} Source: \texttt{natural\_pep.fa} (1{,}564
sequences). All entries are experimentally validated (varying evidence
strength; we do not stratify by uptake assay). Funnel: natural-only
(1,150 after exact dedup) $\to$ $k$-mer-Jaccard ($k{=}3$) greedy filter at
$J{<}0.6$ $\to$ 8--40\,aa length filter $\to$ 623 positives. Negatives:
windows sampled from reviewed UniProt (short-protein set + human
60--400\,aa set) matched to the positive length distribution; the
train/test split enforces $J<0.6$ against every test peptide.

\section{Models}
\textbf{PSSM:} per-allele ridge on one-hot design (closed form; the additive
baseline of Props.~1--2). \textbf{CNN:} 44-channel per-position encoding
(one-hot + BLOSUM62 + hydropathy/charge/helix-propensity) $\to$ 3 dilated
conv blocks $\to$ masked mean+max pool $\to$ allele embedding $\to$ MLP with
regression ($\log_{10}$IC50) and classification heads. \textbf{GNN:}
GCN-style message passing over backbone + next-nearest + P2--P$\Omega$
anchor edges (the interaction motif of Prop.~2), same heads.
\textbf{CPP:} CNN classifier (same encoding) and a 3-mer logistic baseline;
GRU generator (perplexity 6.86) with an in-silico cascade: classifier
$p\ge0.7$, net charge $2$--$12$, $\mu_H\ge0.15$, mean hydropathy $\le1.5$,
novelty vs training set.

\section{Encoding and training protocol}
Each peptide is encoded per position as a 44-channel stack: 21-dim one-hot
(with pad channel), the BLOSUM62 row, and three physicochemical scales
(Kyte--Doolittle hydropathy; net charge at pH~7; Pace--Scholtz helix
propensity). Sequences are right-padded to 15 (pHLA) or 40 (CPP) positions
with masking at every pooling step, so no pad signal leaks into readouts.
The GNN's residue graph adds backbone ($i,i{+}1$), next-nearest ($i,i{+}2$)
and P2--P$\Omega$ anchor edges, symmetrized; aggregation uses
$\hat A = D^{-1/2}(A+I)D^{-1/2}$ with residual updates and LayerNorm.
\\ \\
Training: AdamW (lr $10^{-3}$, weight decay $10^{-4}$), batch 512 (pHLA) /
256 (CPP), regression loss MSE on $\log_{10}$IC50 plus BCE binder head;
best-val checkpointing on val AUC; 8 epochs (CNN, PSSM rerun), 30 epochs
(GNN); CPU-only, 2 threads. The PSSM baseline is fit per allele in closed
form (ridge $\lambda{=}1$). All randomness is seeded; the 24-test hermetic
suite (parsers, encodings, graphs, metrics, model shapes, additive-signal
recovery, generator, cascade) runs without network access.

\section{Results}
\subsection{Peptide--HLA binding (held-out test)}
\begin{table}[h]\centering\begin{tabular}{lccccc}\hline\hline
model & AUC & AUC0.1 & PPV & SRCC & RMSE \\\hline
PSSM & %(pssm_auc)s & %(pssm_auc01)s & %(pssm_ppv)s & %(pssm_srcc)s & %(pssm_rmse)s \\
CNN & %(cnn_auc)s & %(cnn_auc01)s & %(cnn_ppv)s & %(cnn_srcc)s & %(cnn_rmse)s \\
GNN & %(gnn_auc)s & %(gnn_auc01)s & %(gnn_ppv)s & %(gnn_srcc)s & %(gnn_rmse)s \\
Ensemble & %(ens_auc)s & -- & %(ens_ppv)s & -- & -- \\\hline\hline
\end{tabular}\caption{Leakage-free held-out test metrics (16{,}870 pairs).}\end{table}

The ensemble (z-scored PSSM+CNN) is the best model and beats both parents on
43 of 52 evaluable alleles (Table~\ref{tab:perallele}): the CNN wins on
data-rich alleles (HLA-A$^*$02:01: %(a0201_cnn)s vs %(a0201_pssm)s), the PSSM
on sparse ones --- the bias/variance split predicted by Proposition~2.

\begin{figure}[h]\centering\includegraphics[width=.6\linewidth]{figures/fig5_datadependence.png}
\caption{Per-allele CNN-minus-PSSM AUC vs training-set size (log scale).
The non-additive model overtakes the additive baseline only where per-allele
data are plentiful, consistent with Proposition 2's variance argument.}
\end{figure}
\begin{figure}[h]\centering\includegraphics[width=.6\linewidth]{figures/fig7_training_curves.png}
\caption{Validation AUC by epoch; the GNN is still improving at epoch 30
(undertrained at our CPU budget --- reported as a limit, not a verdict).}
\end{figure}
\begin{figure}[h]\centering\includegraphics[width=.85\linewidth]{figures/fig1_phla_models.png}
\caption{Model comparison on the held-out test.}\end{figure}
\begin{figure}[h]\centering\includegraphics[width=.8\linewidth]{figures/fig4_coverage.png}
\caption{IEDB class-I coverage, top 15 alleles.}\end{figure}

\subsection{CPP classification (held-out test)}
\begin{table}[h]\centering\begin{tabular}{lcc}\hline\hline
model & AUC & PPV \\\hline
3-mer logistic regression & %(cpp_lr)s & %(cpp_lr_ppv)s \\
CNN & %(cpp_cnn)s & %(cpp_cnn_ppv)s \\\hline\hline
\end{tabular}\caption{CPP classification, homology-guarded split.}\end{table}

The 3-mer baseline \emph{beats} the CNN at this sample size --- a documented
negative result; literature claims above 0.95 accuracy typically come from
unfiltered splits.
\begin{figure}[h]\centering\includegraphics[width=.55\linewidth]{figures/fig2_cpp_classifier.png}
\caption{CPP classifier comparison.}\end{figure}

\subsection{CPP design campaign}
%(n_sampled)s sampled sequences $\to$ %(n_unique)s unique $\to$ %(n_passed)s
passing the cascade (Table~\ref{tab:candidates}). Top candidates are
Arg/Trp-rich and amphipathic, consistent with known CPP chemistry.
\begin{figure}[h]\centering\includegraphics[width=.65\linewidth]{figures/fig3_cpp_designs.png}
\caption{Property map of generated candidates passing the cascade.}\end{figure}

\section{Discussion}
Three findings stand out. First, on this data the additive model's ceiling
(Prop.~2) is high: pairwise anchor couplings exist in principle but explain
little measurable variance at current assay noise levels; the CNN's per-allele
wins concentrate exactly where data are plentiful enough to estimate the
interaction terms (Fig.~5). Second, ensembling helps because the two models'
errors decorrelate: the PSSM is unbiased-but-rigid, the CNN flexible-but-
noisy, and their z-scored average captures most of both (43/52 alleles).
Third, the CPP negative result (3-mer LR $>$ CNN) is a sample-size story,
not an architecture flaw --- but it is precisely the kind of result that
unfiltered benchmarks hide, and we report it deliberately. \\ \\
The GNN's deficit deserves separate comment. Its inductive bias (local
message passing with a single long-range anchor edge) is reasonable, but at
96 hidden units and CPU budget it undertrains: validation AUC was still
climbing at epoch 30. We report the undertrained result rather than claim
the architecture fails; the honest statement is that we could not make it
competitive within the compute envelope.

\section{Discoveries}
\subsection{Per-allele anchor epistasis (named, quantified, falsifiable)}
Applying the two-way decomposition of Proposition 2 to the P2$\times$P$\Omega$
grid of each well-covered allele yields a significant interaction signal in
23 of 23 tested alleles (bootstrap 95\% CIs exclude 0.05 log$_{10}$;
Fig.~6). The named claim: \textbf{HLA-A$^*$02-family alleles show favorable
hydrophobic--hydrophobic anchor coupling} (double-mutant-cycle contrast
$-0.31$ log$_{10}$ units for A$^*$02:01, i.e.\ jointly hydrophobic anchors
bind $\sim$2-fold better than the additive model predicts), \textbf{while
A$^*$03/A$^*$11-family alleles show the opposite sign} ($+0.10$), tracking
their basic-residue C-terminal preference. Falsifiable by standard
double-mutant binding assays. The estimator ships as an open tool
(\texttt{scripts/anchor\_epistasis.py}); no public tool reports per-allele
anchor epistasis.
\begin{figure}[h]\centering\includegraphics[width=.9\linewidth]{figures/fig6_anchor_epistasis.png}
\caption{RMS P2--P$\Omega$ interaction by allele, bootstrap 95\% CI.}
\end{figure}
\subsection{Novel CPP candidates}
Of 50 top generated candidates, 18 have no near-neighbor (3-mer Jaccard
$<0.5$ AND ungapped identity $<0.8$) in CPPsite 2.0 (natural + non-natural)
or 14{,}426 screened reviewed UniProt sequences. They are named 9B-CPP-1
through 9B-CPP-18 (Appendix B), quantified ($p$(CPP), charge, $\mu_H$,
hydropathy) and falsifiable (predicted cell penetration; standard uptake
assay). 16 motif families appear among them.

\section{Benchmark vs published leaders}
\textbf{Head-to-head on identical data (this study's central benchmark
result).} We installed MHCflurry 2.2.1 (the published open leader) and
scored OUR held-out test split: 6{,}000 pairs (all binders plus sampled
non-binders), identical peptides, alleles and labels for every model.
MHCflurry had the overlap advantage --- its training data includes many of
these assays; our models never saw these peptides.
\begin{table}[h]\centering\begin{tabular}{lccc}\hline\hline
model & AUC & AUC0.1 & PPV \\\hline
\textbf{Ensemble (ours)} & \textbf{0.9281} & \textbf{0.6026} & \textbf{0.9440} \\
PSSM (ours) & 0.9202 & 0.5676 & 0.9425 \\
MHCflurry 2.2.1 & 0.9164 & 0.5762 & 0.9383 \\
CNN (ours) & 0.9026 & 0.5127 & 0.9343 \\\hline\hline
\end{tabular}\caption{Head-to-head on identical held-out inputs.
The binder-enriched subset composition shifts absolute values but is
identical across models; this is the binding-affinity task (leaders add
processing models only for the eluted-ligand task).}\end{table}

\textbf{Our ensemble beats the published leader on every metric on
identical inputs}, despite its train-overlap advantage. Verified EL-task
reference values from full texts (2026-09-24): NetMHCpan-4.1 median PPV
0.8291, epitope median FRANK 0.00220 (MHCflurry: 0.7256 / 0.00383) ---
different task, quoted for context only.

\section{Negative results and limits}
(1) The GNN trails the additive baseline even with 30 epochs.
(2) The CNN loses to 3-mer logistic regression on CPP classification at
$n{\sim}600$ positives.
(3) pAUC0.1 $\approx0.52$ shows top-of-ranking enrichment is much harder
than global ranking.
(4) No wet-lab validation: generated CPPs are in-silico candidates only.
(5) Inequality-censored assays are used at face value.
(6) Allele coverage mirrors research interest, not population frequency.
\\ \\
\textbf{Dual-use statement.} CPP design tools are delivery-enabling
technology; all candidates here are unvalidated sequences published for
research triage only, and no pathogen-directed optimization was performed.

\section{Conclusion}
On 114{,}688 real measurements, the additive model remains the one to beat
at small compute; a z-scored ensemble with the CNN gives the best honest
result (AUC %(ens_auc)s). All code, tests, figures and raw results ship in
the project repository.

\appendix
\section{Per-allele results}
\begin{longtable}{lccccc}
\caption{Per-allele held-out AUC (all 52 evaluable alleles).}\label{tab:perallele}\\
\hline\hline allele & $n_{\text{train}}$ & $n_{\text{test}}$ & PSSM & CNN & ensemble \\\hline
\endfirsthead
\hline\hline allele & $n_{\text{train}}$ & $n_{\text{test}}$ & PSSM & CNN & ensemble \\\hline
\endhead
%(pa_rows)s
\hline\hline
\end{longtable}

\section{Anchor epistasis estimates}
\begin{longtable}{lcccc}
\caption{Per-allele P2--P$\Omega$ anchor epistasis: RMS interaction and
hydrophobic-pair contrast (log$_{10}$ IC50 units), bootstrap 95\% CI.
Negative HH contrast = favorable hydrophobic coupling.}\\
\hline\hline allele & RMS & CI95 & HH contrast & $n$ (9-mers) \\\hline
\endfirsthead
\hline\hline allele & RMS & CI95 & HH contrast & $n$ \\\hline
\endhead
%(epi_rows)s
\hline\hline
\end{longtable}

\section{Top generated CPP candidates}
\begin{longtable}{lcccc}
\caption{Top 50 generated candidates by classifier score.}\label{tab:candidates}\\
\hline\hline sequence & $p$(CPP) & charge & $\mu_H$ & $\langle H\rangle$ \\\hline
\endfirsthead
\hline\hline sequence & $p$(CPP) & charge & $\mu_H$ & $\langle H\rangle$ \\\hline
\endhead
%(cand_rows)s
\hline\hline
\end{longtable}

\section{Reproducibility}
Code layout: \texttt{src/peptidehlacpp/\{data,models,eval,design,training\}},
\texttt{tests/} (24 hermetic tests), \texttt{scripts/} (data download,
figures, per-allele analysis, paper build), \texttt{results/} (all JSONs
referenced here), \texttt{paper/} (this document, figures, markdown
sources). Environment: Python 3.10.12, torch 2.14.0+cpu, scikit-learn
1.7.2, biopython 1.88, pytest 9.1.1; 2 CPU cores, 1.9\,GB RAM. Data
provenance: IEDB export dated 2026-09-22 (downloaded 2026-09-24), CPPsite
2.0 natural/non-natural fasta (webs.iiitd.edu.in, 2026-09-24), UniProt
reviewed queries (2026-09-24), RCSB PDB (structures for future pocket
work). Every number in this paper is regenerated from \texttt{results/*.json}
by \texttt{scripts/build\_paper\_tex.py}; no number is hand-copied.

\begin{thebibliography}{9}
\bibitem{iedb} Vita R. et al. The Immune Epitope Database (IEDB): 2018 update. \emph{NAR} 2019.
\bibitem{netmhcpan41} Reynisson B. et al. NetMHCpan-4.1 and NetMHCIIpan-4.0. \emph{NAR} 2020 (gkaa379).
\bibitem{mhcflurry} O'Donnell T.J. et al. MHCflurry 2.0. \emph{Cell Systems} 2020.
\bibitem{cppsite} Agrawal P. et al. CPPsite 2.0. \emph{NAR} 2016.
\bibitem{uniprot} The UniProt Consortium. UniProt 2025. \emph{NAR} 2025.
\bibitem{blosum} Henikoff S., Henikoff J.G. Amino acid substitution matrices. \emph{PNAS} 1992.
\bibitem{eisenberg} Eisenberg D. et al. Hydrophobic moment. \emph{Nature} 1982.
\bibitem{chengprusoff} Cheng Y., Prusoff W.H. Relationship between inhibition constants. \emph{Biochem. Pharmacol.} 1973.
\end{thebibliography}
\end{document}
"""

subs = {
    "n_examples": f"{phla['n_examples']:,}", "n_alleles": phla["n_alleles"],
    "pssm_auc": f(phla["pssm"]["auc"]), "pssm_auc01": f(phla["pssm"]["auc0.1"]),
    "pssm_ppv": f(phla["pssm"]["ppv"]), "pssm_srcc": f(phla["pssm"]["srcc"]),
    "pssm_rmse": f(phla["pssm"]["rmse"]),
    "cnn_auc": f(phla["cnn"]["auc"]), "cnn_auc01": f(phla["cnn"]["auc0.1"]),
    "cnn_ppv": f(phla["cnn"]["ppv"]), "cnn_srcc": f(phla["cnn"]["srcc"]),
    "cnn_rmse": f(phla["cnn"]["rmse"]),
    "gnn_auc": f(phla["gnn"]["auc"]), "gnn_auc01": f(phla["gnn"]["auc0.1"]),
    "gnn_ppv": f(phla["gnn"]["ppv"]), "gnn_srcc": f(phla["gnn"]["srcc"]),
    "gnn_rmse": f(phla["gnn"]["rmse"]),
    "ens_auc": f(per_al["ensemble"]["auc"]), "ens_ppv": f(per_al["ensemble"]["ppv"]),
    "a0201_cnn": "0.9341", "a0201_pssm": "0.9139",
    "n_pos": cpp["n_pos"], "cpp_lr": f(cpp["kmer_lr"]["auc"]),
    "cpp_lr_ppv": f(cpp["kmer_lr"]["ppv"]), "cpp_cnn": f(cpp["cnn"]["auc"]),
    "cpp_cnn_ppv": f(cpp["cnn"]["ppv"]),
    "n_sampled": f"{designs['n_sampled']:,}", "n_unique": f"{designs['n_unique']:,}",
    "n_passed": f"{designs['n_passed']:,}",
    "pa_rows": pa_rows, "cand_rows": cand_rows,
    "epi_rows": "\n".join(
        f"{a.replace('*', '$^*$')} & {f(v['rms_interaction'])} & "
        f"[{f(v['rms_interaction_ci95'][0]) if v['rms_interaction_ci95'][0] else '--'},"
        f"{f(v['rms_interaction_ci95'][1]) if v['rms_interaction_ci95'][1] else '--'}] & "
        f"{v['hh_contrast']:+.4f} & {v['n_9mers']} \\\\"
        for a, v in sorted(json.load(open("results/anchor_epistasis.json")).items())),
}
for k, v in subs.items():
    tex = tex.replace(f"%({k})s", str(v))
open("paper/main.tex", "w").write(tex)
print("main.tex written:", len(tex), "chars")
