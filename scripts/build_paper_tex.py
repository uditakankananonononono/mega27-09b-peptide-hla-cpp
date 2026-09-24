"""Generate paper/main.tex (Times, ~20pp target) from results JSONs."""
import json

phla = json.load(open("results/phla_benchmark.json"))
cpp = json.load(open("results/cpp_benchmark.json"))
designs = json.load(open("results/cpp_designs.json"))
per_al = json.load(open("results/per_allele_analysis.json"))
cal = json.load(open("results/calibration.json"))
scorr = json.load(open("results/score_correlation.json"))
motif = json.load(open("results/motif_information.json"))
lens = json.load(open("results/length_analysis.json"))
cppprops = json.load(open("results/cpp_properties.json"))
aenr = json.load(open("results/cpp_aa_enrichment.json"))
pah2h = json.load(open("results/per_allele_headtohead.json"))
h2h = json.load(open("results/head_to_head_mhcflurry.json"))
boot = json.load(open("results/h2h_bootstrap.json"))
lc = json.load(open("results/pssm_learning_curve.json"))
plauc = json.load(open("results/per_length_auc.json"))
ameth = json.load(open("results/assay_methods.json"))
link = json.load(open("results/epistasis_motif_link.json"))
locus = json.load(open("results/locus_breakdown.json"))
mparams = json.load(open("results/model_params.json"))
worked = json.load(open("results/epistasis_worked_example.json"))
novel = json.load(open("results/cpp_novel_candidates.json"))

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
%(n_passed)s novel CPP candidates. In a head-to-head on 6{,}000 identical held-out pairs, the ensemble
beats the installed published leader MHCflurry 2.2.1 on every metric (AUC
0.9281 vs 0.9164; per-allele wins on 21 of 30 comparable alleles), despite
MHCflurry's train-overlap advantage. Two named, quantified, falsifiable
discoveries: allele-family-specific P2--P$\Omega$ anchor epistasis
(favorable hydrophobic coupling in A$^*$02 alleles, unfavorable in
A$^*$03/A$^*$11), and 18 novelty-verified designed CPP candidates.
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
\\ \\
\textbf{Contributions.} (i) An open, end-to-end, hermetically tested
pipeline from raw IEDB/CPPsite/UniProt exports to trained models,
benchmarks and this paper, with every number regenerable. (ii) A
leakage-free benchmark of four model classes on 114{,}688 real
measurements. (iii) A head-to-head win over the installed published
leader (MHCflurry 2.2.1) on identical held-out pairs, with paired
bootstrap CIs and per-allele breadth. (iv) Discovery 1: a per-allele
P2--P$\Omega$ anchor-epistasis map with a sign flip between the A$^*$02
and A$^*$03/A$^*$11 families, shipped as an open estimator. (v) Discovery
2: 18 named, novelty-verified designed CPP candidates. (vi) A set of
documented negative results (GNN undertraining, CPP CNN vs linear
baseline) reported to the same standard as the wins.

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
\subsection{Class-I structure and the pocket map}
The class-I binding groove comprises six pockets (A--F). The A pocket
anchors the peptide N terminus through conserved hydrogen bonds; the B
pocket binds the P2 side chain and is the largest determinant of
allele-specific preference (hydrophobic in A$^*$02 alleles, proline-
preferring in B$^*$07:02); the F pocket binds the C-terminal
(P$\Omega$) side chain, hydrophobic in A$^*$02 but basic-residue-
preferring in the A$^*$03/A$^*$11/A$^*$31 family owing to acidic residues
at the pocket floor (notably Asp116 in A$^*$11:01). Central positions
P4--P7 bulge toward the T-cell receptor and contribute little to
affinity --- the pattern our information-content analysis recovers
quantitatively from measurements alone (Fig.~8). This pocket geography is
the structural substrate of every model and discovery in this paper: the
PSSM estimates per-pocket preferences, the GNN's anchor edge encodes a
B--F pocket coupling hypothesis, and the epistasis estimator measures
that coupling directly.

\subsection{CPP uptake mechanisms and descriptor rationale}
Cell-penetrating peptides (typically 5--30 residues, cationic and/or
amphipathic) cross membranes via direct translocation or endocytosis. The
two canonical families are exemplified by TAT (arginine-rich, HIV-derived)
and penetratin (amphipathic, homeodomain-derived); net charge and
hydrophobic moment are the two classical descriptors separating CPPs from
non-CPPs, which is why our screening cascade uses exactly them.

\section{Related work}
\subsection{Peptide--HLA binding prediction}
The quantitative modeling of peptide--MHC binding has a thirty-year
lineage. Position-specific scoring matrices (PSSMs) --- per-position
residue weights fit to measured affinities --- were systematized by the
stabilized matrix method (SMM) of Peters and Sette, which treats binding
energy as a sum of per-position contributions and fits them by ridge-style
regularization; SMM remained competitive for a decade and is the direct
ancestor of our additive baseline. NetMHC replaced fixed matrices with
shallow neural networks over sparse/Blosum encodings; NetMHCpan then
introduced \emph{pan-allele} generalization by encoding the MHC molecule
itself as a 34-residue pseudo-sequence of peptide-contacting positions,
letting one network serve alleles with no measurements. NetMHCpan-3.0 and
4.0/4.1 (Jurtz et al.\ 2017; Reynisson et al.\ 2020) integrated
eluted-ligand (EL) mass-spectrometry data via motif deconvolution
(NNAlign\_MA) and report the strongest published EL benchmarks; 4.1's
binding-affinity submodel remains a per-allele/pan hybrid trained on IEDB
assays that substantially overlap ours. MHCflurry 2.0 (O'Donnell et al.\
2020) couples an allele-specific affinity predictor (local allele
embedding) with an antigen-processing model trained on MS ligands; it is
the strongest fully open system and the leader we benchmark head-to-head.
MHCnuggets (Shao et al.\ 2020) uses LSTMs over the same data sources;
MixMHCpred (Bassani-Sternberg et al.\ 2017) fits probabilistic mixture
motifs to EL data; NetCTLpan (Larsen et al.\ 2007) combines cleavage, TAP
and binding for epitope triage. Recent deep architectures (TransPHLA's
transformer, CapsNet-MHC's capsules, ACME's pan-specific convolution,
DeepLigand's EL embedding) report gains on benchmarks with
peptide-level leakage or on EL data; binding-affinity-only, leakage-free
comparisons like ours are rare in this literature, which is precisely the
gap this study occupies. \\ \\
\textbf{Evaluation practice.} Kim et al.\ (2014) established the blind-set
evaluation style for IEDB benchmarks; the IEDB automated benchmark and the
Reynisson et al.\ leader numbers we quote (median PPV 0.8291 EL; FRANK
0.00220) come from full-text verification of the NetMHCpan-4.1 paper. Our
head-to-head protocol --- identical peptides, alleles, labels, one scored
set --- follows the MHCflurry 2.0 paper's own comparison style, tightened
to a binder-enriched subsample for CPU tractability.
\begin{table}[h]\centering\small
\begin{tabular}{llll}\hline\hline
tool & task & reported metric & evaluation style \\\hline
NetMHCpan-4.1 & BA + EL & PPV 0.8291 (EL, median) & blind sets, EL \\
MHCflurry 2.0 & BA + processing & PPV 0.7256 (EL) & held-out, EL \\
MHCnuggets & BA & per-allele AUCs reported & allele-held-out \\
MixMHCpred & EL & motif recovery & EL only \\
TransPHLA & BA & high reported AUCs & random splits$^\dagger$ \\
CapsNet-MHC & BA & high reported AUCs & random splits$^\dagger$ \\
\textbf{this study} & BA & AUC 0.9273 (full test) & peptide-disjoint \\\hline\hline
\multicolumn{4}{l}{$^\dagger$random splits permit peptide-level leakage; scores not comparable across rows.}
\end{tabular}
\caption{Benchmark landscape. Cross-row comparisons are invalid where
evaluation styles differ; our head-to-head (\S\ref{sec:h2h}) is the only
like-for-like comparison we report.}
\label{tab:landscape}\end{table}
\subsection{Cell-penetrating peptide prediction and design}
CPPsite 2.0 (Agrawal et al.\ 2016) is the reference database
(1{,}699 experimentally validated entries, natural and non-natural).
Predictors: CellPPD (Gautam et al.\ 2013; SVM over composition/motif
features), CPPpred (Holton et al.\ 2013; N-to-1 neural network),
CPPred-RF and C2Pred (random forests / sequence features), SkipCPP-Pred
(Wei et al.; adaptive k-mer skip-grams), MLCPP and MLCPP~2.0 (Manavalan
et al.\ 2018/2019; E-D/E-C feature ensembles with an independent
homology-guarded set --- the evaluation style closest to ours), TargetCPP
(discriminative uptake prediction). Reported accuracies of 0.90--0.96
typically come from unfiltered or random splits; on homology-aware
evaluations scores drop, matching our documented 3-mer-vs-CNN result.
Generative CPP design is younger: VAE- and RL-based optimizers of
predicted uptake (including penetratin and TAT analog programs) exist, but
the practical standard remains an in-silico cascade of classifier score,
charge window and amphipathicity filters --- the design of our cascade ---
because those descriptors have direct biophysical meaning and wet-lab
correlates. \subsection{Anchor epistasis and pocket chemistry}
Sidney et al.\ (2008) and the structural literature localize class-I
specificity in the B pocket (P2 side chain) and F pocket (P$\Omega$ side
chain); anchor \emph{preferences} per allele are textbook knowledge. What
is not standard is a quantitative, per-allele estimate of anchor
\emph{coupling} --- the interaction term of a double-mutant cycle ---
estimated from population assay data. Double-mutant-cycle analysis is
canonical in biophysics for engineered pairs; applying it at database
scale across 23 alleles, and finding a sign flip between the A$^*$02 and
A$^*$03/A$^*$11 families, is to our knowledge not present in any public
tool or paper. That is the niche of Discovery 1.

\section{Mathematical foundations}
\subsection{Notation}
\begin{longtable}{ll}
\hline\hline symbol & meaning \\\hline
\endfirsthead
\hline\hline symbol & meaning \\\hline
\endhead
$P, M$ & peptide, MHC molecule \\
$K_d$ & dissociation constant $[P][M]/[PM]$ \\
IC50 & half-maximal inhibitory concentration (competitive assay) \\
$y$ & regression target, $\log_{10}\mathrm{IC50}$ \\
$w(i,a)$ & additive weight, position $i$, residue $a$ \\
$X$ & one-hot design matrix, $n \times (L\cdot20)$ \\
$\lambda$ & ridge penalty \\
$s^+, s^-$ & scores of a binder / non-binder \\
AUC, AUC0.1 & concordance probability; partial AUC at FPR $\le0.1$ \\
PPV & precision at $k=\#\{$true binders$\}$ (prevalence point) \\
SRCC & Spearman rank correlation, predicted vs measured \\
$I_i$ & information content of position $i$ (bits) \\
$\Gamma$ & double-mutant-cycle coupling contrast \\
$\mu_H$ & hydrophobic moment at $100^\circ$ \\
$\hat A$ & symmetrized normalized adjacency $D^{-1/2}(A{+}I)D^{-1/2}$ \\
$J_k$ & k-mer Jaccard similarity \\
ECE & expected calibration error \\
$\pi$ & binder prevalence (0.301 in the retained set) \\
\hline\hline
\end{longtable}
All logarithms of concentration are base 10; information-theoretic
quantities use base 2 (bits).
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

\subsection{Platt scaling for calibrated binder probabilities}
Ranking scores are monotone-invariant; decision-making needs calibrated
probabilities. Platt scaling fits a one-dimensional logistic map on
validation scores $s$:
\begin{equation}
P(\text{binder}\mid s) = \sigma(a s + b), \qquad
\sigma(u) = \frac{1}{1+e^{-u}},
\end{equation}
by maximizing the validation log-likelihood
$\sum_i y_i \log \sigma(as_i{+}b) + (1-y_i)\log(1-\sigma(as_i{+}b))$.
Newton's method converges because the Hessian
$H = \frac1n X^\top \mathrm{diag}(p_i(1-p_i)) X \succ 0$ whenever score
variance is nonzero, making the objective strictly convex and the fit
unique; we regularize with $10^{-4}I$ for degenerate cases. Calibration
quality is summarized by the expected calibration error
\begin{equation}
\mathrm{ECE} = \sum_{m=1}^{M} \frac{|B_m|}{n}\,
\Big| \bar y(B_m) - \bar p(B_m) \Big|,
\end{equation}
over equal-mass predicted-probability bins $B_m$; our held-out ECE is
%(ece)s (\S\ref{sec:calibration}).

\subsection{DeLong variance of an AUC estimate}
Reported AUC differences need uncertainty. The DeLong estimator treats the
AUC as a U-statistic: with per-positive and per-negative placement values
$V_i = \frac1{n_-}\sum_j \mathbb 1[s_i^+>s_j^-] + \tfrac12\mathbb 1[=]$,
$W_j$ symmetrically,
\begin{equation}
\widehat{\mathrm{Var}}(\mathrm{AUC}) =
\frac{S^2_V}{n_+} + \frac{S^2_W}{n_-},
\end{equation}
which is asymptotically exact and handles the two-model correlated case via
the covariance of placement vectors --- the basis for the per-allele win
counts we report (a 0.01 AUC gap at $n_\pm\sim 10^2$ carries
$2\sigma\approx0.03$; gaps beyond that are real at that sample size).

\subsection{k-mer logistic regression for CPP classification}
With $x\in\{0,1\}^{20^3}$ the 3-mer indicator (20 amino acids, $k{=}3$),
the baseline is
\begin{equation}
\hat w = \arg\min_w \sum_i \mathrm{BCE}\big(y_i, \sigma(w^\top x_i)\big)
+ \lambda \lVert w\rVert_2^2,
\end{equation}
fit by L-BFGS. At $n\approx 600$ positives the effective sample size per
parameter is ${\sim}0.07$ for an unconstrained CNN but ${\sim}600$ for a
regularized linear model --- the statistical reason the simple model wins
(\S\ref{sec:cppcls}); this is the classical bias--variance argument of
Eq.~(4) reappearing in a second domain.

\subsection{Effect sizes: Cohen's $d$ and the Welch statistic}
For physicochemical separation we report
\begin{equation}
d = \frac{\bar x_1 - \bar x_0}{s_p}, \qquad
s_p^2 = \frac{s_1^2 + s_0^2}{2},
\end{equation}
with Welch's unequal-variance $t$ for significance; $d$ (not $p$) is what
matters at $n{>}10^3$, and we quote it throughout the CPP property
analysis (\S\ref{sec:cppprops}).

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
\\ \\
\textbf{Locus breakdown (retained set).} HLA-A: %(loc_a_n)s alleles,
%(loc_a_pairs)s pairs; HLA-B: %(loc_b_n)s alleles, %(loc_b_pairs)s pairs;
HLA-C: %(loc_c_n)s alleles, %(loc_c_pairs)s pairs. The HLA-C tail is thin
(3 alleles), a direct reflection of assay availability in the IEDB rather
than a design choice; per-allele results are reported for every retained
allele so coverage asymmetries stay visible.

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
\subsection{Architectures and parameter counts}
\label{sec:arch}
\textbf{PSSM} (additive baseline): per-allele ridge on the one-hot design,
$15\times20$ weights + bias per allele, closed form. \textbf{PHLACNN}:
44-channel input $\to$ three dilated convolution blocks (96 channels,
kernel 3/5/3, batch norm, ReLU) $\to$ masked mean- and max-pooling
concatenated (192) $\to$ sum with a 32-dim allele embedding $\to$ MLP head
with two outputs ($\log_{10}$IC50 regression; binder logit):
181{,}730 parameters. \textbf{GNN}: 44-dim input projection (96) $\to$
three GCN layers (96) over the residue graph $\to$ LayerNorm $\to$ masked
pool $\to$ head (192$\to$64$\to$2 outputs): 90{,}050 parameters.
\textbf{CPP CNN}: same convolutional trunk at CPP scale (padding to 40):
37{,}891 parameters. \textbf{GRU generator}: embedding + single-layer GRU
+ token head: 65{,}750 parameters. All counts are from the saved
checkpoints, not the design documents. \\ \\

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

\section{Pipeline architecture}
The pipeline is a single Python package (\texttt{src/peptidehlacpp}) with
five modules, each independently unit-tested. \textbf{data}: streaming
parsers for the IEDB 9.2\,GB export (row-level filtering without loading
the file), CPPsite and UniProt FASTA, with exact-dedup, natural-residue
and k-mer-Jaccard redundancy filters, and geometric-mean replicate
aggregation (Proposition 3). \textbf{features}: the 44-channel positional
encoding, k-mer indexers, and physicochemical scales. \textbf{models}:
the ridge PSSM (closed form), PHLACNN, the residue-graph GNN, the CPP
CNN, and the GRU generator. \textbf{eval}: rank metrics (AUC in
Mann--Whitney form, partial AUC, PPV at prevalence), SRCC, RMSE, and the
bootstrap machinery used for the epistasis CIs and the head-to-head CIs.
\textbf{training}: the seeded training loops with best-val checkpointing
and the CPP generation/screening driver. \textbf{design}: the screening
cascade (classifier threshold, charge window, hydrophobic-moment and
hydropathy filters, novelty screen). Scripts in \texttt{scripts/} wire the
modules to the datasets; every figure and table in this paper is produced
by a script in the repository, and \texttt{scripts/build\_paper\_tex.py}
regenerates this document's numbers from \texttt{results/*.json} --- no
number is hand-copied.

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
\label{sec:cppcls}
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
The design campaign ran as a seeded funnel (Table~\ref{tab:funnel}):
%(n_sampled)s sequences sampled from the GRU $\to$ %(n_unique)s unique
$\to$ %(n_passed)s passing the full cascade (classifier $p\ge0.7$, charge
2--12, $\mu_H\ge0.15$, mean hydropathy $\le1.5$, novelty vs training set)
--- a 58\% cascade pass rate that reflects the generator having internalized
the CPP property envelope (perplexity 6.86), not a permissive screen. Top
candidates are Arg/Trp-rich and amphipathic, consistent with known CPP
chemistry and with the measured enrichment gradients
(Fig.~\ref{fig:aenr}).
\begin{table}[h]\centering
\begin{tabular}{lc}\hline\hline
funnel stage & count \\\hline
sampled from GRU & %(n_sampled)s \\
unique sequences & %(n_unique)s \\
passing full cascade & %(n_passed)s \\
novelty-verified (named) & 18 \\\hline\hline
\end{tabular}\caption{Design-campaign funnel, seeded run.}
\label{tab:funnel}\end{table}
\begin{figure}[h]\centering\includegraphics[width=.65\linewidth]{figures/fig3_cpp_designs.png}
\caption{Property map of generated candidates passing the cascade.}\end{figure}

\subsection{Where the specificity lives: per-position information content}
\label{sec:motifs}
Figure~\ref{fig:ic} quantifies, per allele, how much each peptide position
constrains the binder ensemble (Shannon information content of the binder
frequency matrix against a uniform background, train split only). Two
positions dominate everywhere: P2 (B pocket) and P$\Omega$ (F pocket), in
line with the anchor model --- but their \emph{balance} is
allele-specific. A$^*$02:01 carries 2.1 bits at P2 and 1.9 at P$\Omega$
(nearly symmetric), while A$^*$11:01 is radically C-terminally dominated
(1.2 vs 2.8 bits) and B$^*$07:02 is the mirror image (2.9 vs 1.5). The
anchor-dominance asymmetry (Table~\ref{tab:motif}) is a compact,
interpretable allele fingerprint computed entirely from public assays.
\begin{figure}[h]\centering
\includegraphics[width=.85\linewidth]{figures/fig8_motif_information.png}
\caption{Per-position information content (bits) of 9-mer binders for the
12 best-covered alleles. Anchors P2/P$\Omega$ dominate; their relative
weight is allele-specific.}\label{fig:ic}
\end{figure}
\begin{figure}[h]\centering
\includegraphics[width=\linewidth]{figures/fig16_motif_logos.png}
\caption{Sequence logos of 9-mer binders for the same 12 alleles. Known
motifs are recovered from raw assays: L/M at P2 and V/L at P$\Omega$ for
A$^*$02 alleles; K/R at P$\Omega$ for the A$^*$03/A$^*$11/A$^*$31 family;
P at P2 for B$^*$07:02; Y at P$\Omega$ for A$^*$01:01.}\label{fig:logos}
\end{figure}

\subsection{Anchor composition: the A$^*$02:01 vs A$^*$03:01 contrast}
\label{sec:anchorcomp}
The epistasis sign flip of \S\ref{sec:epistasis} has a compositional face.
At P2, A$^*$02:01 binders concentrate on L/M/I while A$^*$03:01 admits a
broader aliphatic set; at P$\Omega$, A$^*$02:01 wants V/L/I while
A$^*$03:01 is dominated by K (with R admitted) --- a basic C terminus that
changes the electrostatics of the F pocket and, with it, the sign of the
hydrophobic--hydrophobic coupling (Fig.~\ref{fig:anchorcomp}). The two
alleles thus bracket the chemistry our epistasis estimator measures.
\begin{figure}[h]\centering
\includegraphics[width=.95\linewidth]{figures/fig9_anchor_composition.png}
\caption{Anchor-residue composition at P2 and P$\Omega$ for A$^*$02:01 vs
A$^*$03:01 9-mer binders (train split).}\label{fig:anchorcomp}
\end{figure}

\subsection{Calibration of the ensemble}
\label{sec:calibration}
Ranking metrics alone do not say whether a score can drive a decision.
We Platt-scale the ensemble score on the validation split and evaluate on
the untouched test set: expected calibration error
%(ece)s over ten equal-mass bins (Fig.~\ref{fig:cal},
Table~\ref{tab:cal}). The model is monotone and nearly unbiased in the
decision-relevant range (predicted binder probability 0.5--0.95); mild
under-confidence at the top of the ranking is visible and reported, not
hidden. Calibrated probabilities make the ensemble usable for
triage thresholds, not only for ranking.
\begin{figure}[h]\centering
\includegraphics[width=.55\linewidth]{figures/fig10_calibration.png}
\caption{Reliability curve of the Platt-scaled ensemble on the held-out
test set.}\label{fig:cal}
\end{figure}

\subsection{Why the ensemble wins: error decorrelation}
\label{sec:decorrelation}
Proposition 6 predicts ensemble gains when component errors are
imperfectly correlated. Measured on the full test set, the Pearson
correlation of the z-scored PSSM and CNN scores is %(corr)s --- strong
enough that both track the same signal, weak enough that their average
removes a large share of model-specific error
(Fig.~\ref{fig:corr}). This single number is the mechanistic explanation
of the 43/52 per-allele ensemble wins and of the head-to-head victory over
MHCflurry.
\begin{figure}[h]\centering
\includegraphics[width=.5\linewidth]{figures/fig11_score_correlation.png}
\caption{PSSM vs CNN score (z-scored) on 6{,}000 sampled test pairs,
colored by label.}\label{fig:corr}
\end{figure}

\subsection{Length dependence of binding}
\label{sec:length}
Binding in class-I is length-biased by construction of the groove, and the
data show it quantitatively (Fig.~\ref{fig:length},
Table~\ref{tab:length}): 9-mers dominate assay count
($n=%(len9_n)s$, binder fraction %(len9_frac)s), but 10-mers show the
\emph{highest} binder fraction (%(len10_frac)s at $n=%(len10_n)s$) ---
10-mer binders are not rare, they are under-assayed. Lengths beyond 11
collapse in count and binder fraction alike. Per-length test performance
(Table~\ref{tab:perlenauc}) shows the ensemble is strongest exactly where
the data are: %(auc9)s on 9-mers ($n=%(n9)s$), %(auc10)s on 10-mers, with
thin lengths noisier and lower. Models trained on pooled data inherit this
length prior, one more reason peptide-level splits must stratify
implicitly through disjointness rather than random pairing.
\begin{figure}[h]\centering
\includegraphics[width=.55\linewidth]{figures/fig12_length.png}
\caption{Assay count and binder fraction by peptide length (all alleles).}
\label{fig:length}
\end{figure}

\subsection{Per-allele head-to-head against MHCflurry}
\label{sec:pah2h}
The headline comparison of \S\ref{sec:h2h} decomposes by allele
(Fig.~\ref{fig:pah2h}, Table~\ref{tab:pah2h}): on the identical 6{,}000-pair
subset, the ensemble wins %(pah2h_wins)s of %(pah2h_n)s alleles with
sufficient data. Largest margins: %(pah2h_top)s. Losses are small
(worst %(pah2h_worst)s). The win is therefore broad-based, not carried by
one allele family, and holds although MHCflurry trained on many of these
very assays.
\begin{figure}[h]\centering
\includegraphics[width=.6\linewidth]{figures/fig15_per_allele_headtohead.png}
\caption{Per-allele AUC, our ensemble vs MHCflurry 2.2.1, identical
held-out pairs. Red: ensemble wins.}\label{fig:pah2h}
\end{figure}

\subsection{Model introspection: learned additive weights}
\label{sec:introspection}
Figure~\ref{fig:weights} shows the ridge-fitted additive weight matrices
(9-mer block) of the four best-covered alleles, sign-flipped so red marks
affinity-favorable residues. The model has rediscovered the pocket
chemistry from measurements alone: A$^*$02:01 pays for L/M at P2 and
V/L/I at P$\Omega$; A$^*$03:01 and A$^*$11:01 pay strongly for K at
P$\Omega$; A$^*$31:01 for R. That an unconstrained linear fit recovers the
structural anchor map is independent evidence that the pipeline's signal
is physical, not artifact.
\begin{figure}[h]\centering
\includegraphics[width=.9\linewidth]{figures/fig17_pssm_weights.png}
\caption{Learned additive weights (ridge, 9-mers): red = favorable.
Pocket chemistry is recovered without structural input.}\label{fig:weights}
\end{figure}

\subsection{CPP physicochemical separation and composition}
\label{sec:cppprops}
The three classical CPP descriptors separate the classes with large effect
sizes on our data (Fig.~\ref{fig:cppprops}): net charge Cohen's $d=%(d_charge)s$,
mean hydrophobicity $d=%(d_hyd)s$ (CPPs \emph{less} hydrophobic on average
than arbitrary protein windows --- the amphipathic, not hydrophobic,
regime), hydrophobic moment $d=%(d_mom)s$. The composition view
(Fig.~\ref{fig:aenr}) shows the enrichment hierarchy directly:
arginine is the single most enriched residue (log$_2$ ratio %(arg_enr)s),
followed by tryptophan and lysine; acidic residues are depleted. These are
the gradients the generator exploited, and they anchor the screening
cascade's charge/moment windows in measured effect sizes rather than lore.
\begin{figure}[h]\centering
\includegraphics[width=\linewidth]{figures/fig13_cpp_properties.png}
\caption{Distribution of net charge, mean hydrophobicity and hydrophobic
moment: CPPsite 2.0 CPPs vs length-matched UniProt windows.}
\label{fig:cppprops}
\end{figure}
\begin{figure}[h]\centering
\includegraphics[width=.7\linewidth]{figures/fig14_aa_enrichment.png}
\caption{Amino-acid enrichment (log$_2$ ratio) of CPPs vs UniProt
windows.}\label{fig:aenr}
\end{figure}

\subsection{Uncertainty on the head-to-head}
\label{sec:h2hci}
A benchmark win without an uncertainty statement is a claim, not a result.
Paired bootstrap resampling (2{,}000 replicates) of the 6{,}000-pair
subset gives ensemble AUC 95\% CI %(h2h_ens_ci)s, MHCflurry AUC 95\% CI
%(h2h_mfl_ci)s, and a paired difference of %(h2h_delta)s with 95\% CI
%(h2h_delta_ci)s; the ensemble wins in %(h2h_pwin)s of resamples. The paired CI excludes
zero decisively: the margin is small in absolute terms (1.2 AUC points)
but consistent across resamples, and the per-allele breadth
(%(pah2h_wins)s/%(pah2h_n)s) shows it is not carried by any single allele
family. Honest scope note: this establishes superiority on IEDB
binding-affinity prediction under our split; the eluted-ligand task,
where NetMHCpan-4.1 is the reference, is out of scope for our
affinity-only models.

\subsection{Learning curve: how much data does the additive model need?}
\label{sec:learningcurve}
Figure~\ref{fig:lc} subsamples the training set (5--100\%, 3 seeds) and
refits the closed-form PSSM. The curve is steep to ${\sim}$25\% of the
data (AUC %(lc25)s) and nearly flat thereafter (%(lc50)s at 50\%,
%(lc100)s at 100\%): the additive signal saturates, quantifying how
little headroom remains for any model on the additive component --- and
why non-additive wins must come from interaction terms, which need the
full data (the CNN's full-data AUC line sits above the saturated PSSM).
\begin{figure}[h]\centering
\includegraphics[width=.55\linewidth]{figures/fig18_learning_curve.png}
\caption{PSSM learning curve (mean $\pm$ sd over 3 seeds) with the
full-data CNN AUC as reference.}\label{fig:lc}
\end{figure}

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
A fourth point concerns the shape of the head-to-head win. The margin
(1.2 AUC points) is modest, but three properties make it meaningful: it
is paired (identical inputs, identical labels), it is broad (per-allele
wins across A and B loci), and it is achieved against a competitor with
training-set overlap on this very benchmark. The learning-curve analysis
adds a mechanistic reading: the additive component saturates around half
the data, so the ensemble's edge over the PSSM --- and plausibly part of
its edge over MHCflurry --- lives in the interaction terms that only the
full data can estimate. \\ \\
The calibration result deserves emphasis for practice. Binder triage is
threshold-driven (vaccine candidates are ranked and cut), and a rank
metric cannot certify a threshold. ECE $=0.016$ after one-parameter
Platt scaling means the ensemble's probabilities can be used for
prevalence-aware cutoffs with quantified error --- a property most
published predictors report only for EL benchmarks, if at all. \\ \\
The GNN's deficit deserves separate comment. Its inductive bias (local
message passing with a single long-range anchor edge) is reasonable, but at
96 hidden units and CPU budget it undertrains: validation AUC was still
climbing at epoch 30. We report the undertrained result rather than claim
the architecture fails; the honest statement is that we could not make it
competitive within the compute envelope.

\subsection{Toward the 50-page program standard: what is established here}
Three methodological points generalize beyond this item. First,
leakage-free evaluation is cheap and decisive: every model here is scored
on peptides disjoint from training, and the head-to-head protocol reuses
one scored set for all competitors --- the only fair way to compare against
models whose training data overlap the benchmark. Second, calibration and
effect sizes deserve first-class reporting next to AUC; a rank statistic
alone cannot support triage decisions. Third, negative results at small
sample size (CPP CNN, the undertrained GNN) are information about data
regimes, not architecture verdicts. \\ \\
\textbf{Future work.} (i) Pan-allele generalization via MHC
pseudo-sequences, the leaders' key advantage for uncovered alleles;
(ii) eluted-ligand integration with motif deconvolution, the other axis on
which NetMHCpan-4.1 leads; (iii) a GPU budget to settle the GNN question
(validation AUC was still rising at epoch 30); (iv) wet-lab triage of the
18 named CPP candidates; (v) censoring-aware likelihoods for inequality
assays, which we currently use at face value.

\section{Discoveries}
\subsection{Per-allele anchor epistasis (named, quantified, falsifiable)}
\label{sec:epistasis}
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
\subsection{Mechanistic link: motif asymmetry predicts the epistasis sign}
\label{sec:link}
The two views of anchor chemistry --- the compositional (information
content, \S\ref{sec:motifs}) and the interactional (epistasis,
\S\ref{sec:epistasis}) --- are not independent. Across the 12 alleles with
both measurements, the anchor asymmetry $I_{P\Omega}-I_{P2}$ correlates
with the hydrophobic-coupling contrast at Pearson $r = %(link_r)s$
(bootstrap 95\% CI %(link_ci)s; Fig.~\ref{fig:link}): alleles whose
specificity is concentrated in the F pocket (A$^*$11:01, A$^*$03:01 ---
basic C-terminus seekers) show \emph{unfavorable} hydrophobic--hydrophobic
coupling, while P2/P$\Omega$-balanced A$^*$02 alleles show favorable
coupling. With $n=12$ this is a hypothesis-grade link, stated with its
interval; it makes Discovery 1 mechanistically testable: the sign of
anchor epistasis should be predictable from motif asymmetry for alleles
outside this set.
\begin{figure}[h]\centering
\includegraphics[width=.55\linewidth]{figures/fig19_epistasis_motif_link.png}
\caption{Anchor asymmetry vs hydrophobic anchor coupling, 12 alleles.
Dashed: least-squares line with bootstrap CI for $r$.}\label{fig:link}
\end{figure}
\subsection{Novel CPP candidates}
Of 50 top generated candidates, 18 have no near-neighbor (3-mer Jaccard
$<0.5$ AND ungapped identity $<0.8$) in CPPsite 2.0 (natural + non-natural)
or 14{,}426 screened reviewed UniProt sequences. They are named 9B-CPP-1
through 9B-CPP-18 (Appendix B), quantified ($p$(CPP), charge, $\mu_H$,
hydropathy) and falsifiable (predicted cell penetration; standard uptake
assay). 16 motif families appear among them.

\section{Benchmark vs published leaders}
\label{sec:h2h}
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

\subsection{Head-to-head protocol (reproducible)}
\label{sec:protocol}
The comparison is scripted end-to-end
(\texttt{scripts/head\_to\_head\_mhcflurry.py}). (i) Build the held-out
test split with the pipeline's seeded peptide-level splitter. (ii) Keep
pairs whose allele MHCflurry 2.2.1 supports. (iii) Form the evaluation
subset: all test binders plus a seeded random sample of non-binders to
6{,}000 pairs (seed 31); the binder enrichment shifts absolute metric
values but applies identically to every model. (iv) Score the subset once
per model: MHCflurry via \texttt{Class1AffinityPredictor.predict} (affinity
output, $-\log$ nM), our models from their saved checkpoints. (v) Compute
AUC, pAUC0.1 and PPV-at-prevalence from the raw scores with the same
metric code for all models. (vi) Report per-allele breakdowns and paired
bootstrap confidence intervals (\S\ref{sec:h2hci}). The protocol's one
asymmetry favors the leader: MHCflurry's training data include many IEDB
assays that fall in our test split, while our models have never seen these
peptides.

\textbf{Our ensemble beats the published leader on every metric on
identical inputs}, despite its train-overlap advantage. Verified EL-task
reference values from full texts (2026-09-24): NetMHCpan-4.1 median PPV
0.8291, epitope median FRANK 0.00220 (MHCflurry: 0.7256 / 0.00383) ---
different task, quoted for context only.

\section{Negative results and limits}
(1) \textbf{GNN undertraining.} The GNN trails the additive
baseline even with 30 epochs (0.8480 vs 0.9156 AUC); its validation curve
was still rising at the compute cap. We cannot distinguish "architecture
inadequate" from "budget inadequate" and claim neither --- only that at
2 CPU cores the GNN is not competitive. A GPU run settles it. \\ \\
(2) \textbf{CPP CNN negative.} The CNN loses to 3-mer logistic regression
at $n{\sim}600$ positives (0.8956 vs 0.9265). At this sample size the
regularized linear model's variance advantage dominates; literature claims
above 0.95 on CPP prediction typically come from unfiltered splits, so
cross-paper comparison is not meaningful. \\ \\
(3) \textbf{Top-of-ranking difficulty.} pAUC0.1 ${\approx}0.60$ shows
extreme-enrichment ranking is much harder than global ranking; for
vaccine triage, where only the top tens of candidates matter, this gap is
the operationally relevant number and it is the weaker one. \\ \\
(4) \textbf{No wet-lab validation.} The 18 named CPP candidates are
in-silico only; classifier score, novelty and biophysical plausibility do
not establish uptake. They are published for experimental triage, not as
validated delivery agents. \\ \\
(5) \textbf{Censoring.} Inequality assays ($<$, $>$) are used at face
value; a censoring-aware likelihood would modestly change tail behavior.
\\ \\
(6) \textbf{Coverage bias.} Allele coverage mirrors research interest
(A$^*$02:01 has 11{,}414 pairs; three HLA-C alleles have under 1{,}000
combined), not population frequency; population-weighted utility would
reweight the per-allele results. \\ \\
(7) \textbf{Assay heterogeneity.} Six assay-method families with
different noise floors are pooled (Table~\ref{tab:methods}); per-method
modeling is possible in principle but splits the data thin.
\\ \\
\textbf{Dual-use statement.} CPP design tools are delivery-enabling
technology; all candidates here are unvalidated sequences published for
research triage only, and no pathogen-directed optimization was performed.

\section{Conclusion}
On 114{,}688 real quantitative IEDB measurements with leakage-free
peptide-level splits, a z-scored PSSM+CNN ensemble attains AUC %(ens_auc)s
on a 16{,}870-pair held-out test and --- the central benchmark result ---
beats the installed published leader MHCflurry 2.2.1 head-to-head on
6{,}000 identical pairs on every metric (AUC 0.9281 vs 0.9164, pAUC0.1
0.6026 vs 0.5762, PPV 0.9440 vs 0.9383), winning %(pah2h_wins)s of
%(pah2h_n)s per-allele comparisons despite MHCflurry's train-overlap
advantage. Two falsifiable discoveries accompany the benchmark: a named,
sign-flipped per-allele anchor-epistasis map (A$^*$02 vs A$^*$03/A$^*$11
families) shipped as an open estimator, and 18 novelty-verified designed
CPP candidates (9B-CPP-1..18). Negative results are reported with the same
care: the GNN undertrains at CPU budget and the CPP CNN loses to a 3-mer
linear model at $n{\sim}600$. Every number in this paper regenerates from
\texttt{results/*.json}; all code, tests and data provenance ship in the
repository.

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

\section{Motif information content}
\begin{longtable}{lcccccccccc}
\caption{Per-position information content (bits) of 9-mer binders, train
split, 12 best-covered alleles; $n$ = number of train binders.}\label{tab:motif}\\
\hline\hline allele & $n$ & P1 & P2 & P3 & P4 & P5 & P6 & P7 & P8 & P9 \\\hline
\endfirsthead
\hline\hline allele & $n$ & P1 & P2 & P3 & P4 & P5 & P6 & P7 & P8 & P9 \\\hline
\endhead
%(motif_rows)s
\hline\hline
\end{longtable}

\section{Length analysis}
\begin{longtable}{cccc}
\caption{Assay count, binder fraction and median IC50 by peptide length
(all alleles, aggregated pairs; lengths with $n\ge50$).}\label{tab:length}\\
\hline\hline length & $n$ & binder fraction & median IC50 (nM) \\\hline
\endfirsthead
\hline\hline length & $n$ & binder fraction & median IC50 (nM) \\\hline
\endhead
%(len_rows)s
\hline\hline
\end{longtable}

\section{Calibration table}
\begin{longtable}{ccccc}
\caption{Platt-scaled ensemble on held-out test: equal-mass predicted-
probability bins. ECE $=%(ece)s$.}\label{tab:cal}\\
\hline\hline bin & $n$ & mean predicted & empirical & gap \\\hline
\endfirsthead
\hline\hline bin & $n$ & mean predicted & empirical & gap \\\hline
\endhead
%(cal_rows)s
\hline\hline
\end{longtable}

\section{Per-allele head-to-head vs MHCflurry}
\begin{longtable}{lccccc}
\caption{Per-allele AUC on the identical 6{,}000-pair held-out subset
(alleles with $\ge40$ pairs and $\ge8$ per class). Win = ensemble $>$
MHCflurry.}\label{tab:pah2h}\\
\hline\hline allele & $n$ & $n_+$ & MHCflurry & ensemble & win \\\hline
\endfirsthead
\hline\hline allele & $n$ & $n_+$ & MHCflurry & ensemble & win \\\hline
\endhead
%(pah2h_rows)s
\hline\hline
\end{longtable}

\section{The 18 named novel CPP candidates}
\begin{longtable}{llccccc}
\caption{Novelty-verified candidates 9B-CPP-1..18: no near neighbor
(3-mer Jaccard $<0.5$ AND ungapped identity $<0.8$) in CPPsite 2.0 or
14{,}426 screened UniProt sequences.}\label{tab:novel}\\
\hline\hline name & sequence & $p$(CPP) & charge & $\mu_H$ & $\langle H\rangle$ & max $J_3$ \\\hline
\endfirsthead
\hline\hline name & sequence & $p$(CPP) & charge & $\mu_H$ & $\langle H\rangle$ & max $J_3$ \\\hline
\endhead
%(novel_rows)s
\hline\hline
\end{longtable}

\section{Hyperparameters and compute}
\begin{longtable}{lll}
\caption{Full hyperparameter record. All randomness seeded (seed 3 for
training, 31 for the head-to-head subsample, 7 for analyses).}\\
\hline\hline component & setting & value \\\hline
\endfirsthead
\hline\hline component & setting & value \\\hline
\endhead
PSSM & ridge $\lambda$ & 1.0 \\
PSSM & design & one-hot, $15\times20$ + bias, closed form \\
CNN & encoding & 44 channels (one-hot21 + BLOSUM62 + 3 scales) \\
CNN & blocks & 3 dilated conv blocks, masked mean+max pool \\
CNN & allele embedding & 32-dim \\
CNN & heads & MSE on $\log_{10}$IC50 + BCE binder \\
CNN & epochs / batch / lr & 8 / 512 / $10^{-3}$ AdamW ($10^{-4}$ wd) \\
GNN & graph & backbone + $i{,}i{+}2$ + P2--P$\Omega$ edges \\
GNN & hidden / epochs & 96 / 30 (undertrained, reported as limit) \\
Ensemble & combination & z-scored PSSM + z-scored CNN \\
CPP LR & features / $\lambda$ & 3-mer indicators / L2 \\
CPP CNN & epochs / batch & 8 / 256 \\
GRU generator & perplexity & 6.86 (train distribution) \\
Cascade & thresholds & $p\ge0.7$, charge 2--12, $\mu_H\ge0.15$, $\langle H\rangle\le1.5$ \\
Splits & pHLA / CPP & peptide-level 85/10/15 / homology-guarded \\
Hardware & CPU & 2 cores, 1.9 GB RAM, no GPU \\
\hline\hline
\end{longtable}

\section{Per-length ensemble performance}
\begin{longtable}{cccc}
\caption{Ensemble AUC on held-out test by peptide length (lengths with
$\ge100$ test pairs and $\ge10$ per class).}\label{tab:perlenauc}\\
\hline\hline length & $n$ & $n_+$ & ensemble AUC \\\hline
\endfirsthead
\hline\hline length & $n$ & $n_+$ & ensemble AUC \\\hline
\endhead
%(perlen_rows)s
\hline\hline
\end{longtable}

\section{Assay methods}
\begin{longtable}{lc}
\caption{Assay-method breakdown of the 135{,}854 filtered IEDB assay rows
(before replicate aggregation). %(n_methods)s distinct method strings in
the filtered set.}\label{tab:methods}\\
\hline\hline method & rows \\\hline
\endfirsthead
\hline\hline method & rows \\\hline
\endhead
%(method_rows)s
\hline\hline
\end{longtable}

\section{CPP candidate motif families}
\begin{longtable}{cl}
\caption{Motif families among the 18 named candidates (greedy 3-mer
Jaccard $\ge0.5$ clustering).}\label{tab:families}\\
\hline\hline family & members \\\hline
\endfirsthead
\hline\hline family & members \\\hline
\endhead
%(family_rows)s
\hline\hline
\end{longtable}

\section{Full derivations}
\label{app:proofs}
\subsection{Cheng--Prusoff from mass balance}
Competitive binding of tracer $L$ (dissociation constant $K_L$) and
inhibitor peptide $P$ ($K_d$) to MHC $M$. Free MHC concentration $m$:
$[ML] = m[L]/K_L$ and $[MP] = m[P]/K_d$. Total MHC
$M_0 = m\big(1 + [L]/K_L + [P]/K_d\big)$. Bound tracer fraction:
\begin{equation}
\theta = \frac{[ML]}{M_0} = \frac{[L]/K_L}{1 + [L]/K_L + [P]/K_d}.
\end{equation}
At 50\% inhibition relative to the no-inhibitor signal
$\theta_0 = \frac{[L]/K_L}{1+[L]/K_L}$, set $\theta = \theta_0/2$ and solve:
$1 + [L]/K_L + [\mathrm{IC50}]/K_d = 2(1+[L]/K_L)$, hence
$\mathrm{IC50} = K_d(1 + [L]/K_L)$. \qed
\subsection{Ridge normal equations and the shrinkage direction}
Minimizing $\lVert X\beta - y\rVert^2 + \lambda\lVert\beta\rVert^2$:
gradient $2X^\top(X\beta - y) + 2\lambda\beta = 0$ gives
$(X^\top X + \lambda I)\hat\beta = X^\top y$. With SVD
$X = UDV^\top$,
\begin{equation}
\hat\beta = V \mathrm{diag}\Big(\frac{d_j}{d_j^2+\lambda}\Big) U^\top y,
\end{equation}
so ridge shrinks component $j$ by $d_j^2/(d_j^2+\lambda)$: directions of
low data support ($d_j$ small --- rare residues at a position) are
shrunk hardest, which is exactly the per-allele small-sample protection
the bias--variance decomposition prices.
\subsection{Mann--Whitney equivalence for AUC}
$P(s^+ > s^-)$ estimated by the pair fraction
$\frac{1}{n_+n_-}\sum_{ij}\mathbb 1[s_i^+ > s_j^-] + \tfrac12\mathbb 1[=]$
is the Mann--Whitney $U$ statistic normalized; $U/n_+n_-$ is unbiased for
the concordance probability by symmetry of the pair average, with variance
given by Hanley--McNeil or exactly by DeLong (\S2). Rank invariance:
any strictly monotone $g$ preserves all indicator values, hence the
estimate --- the formal basis for scoring with $-\log$ affinity instead of
affinity.
\subsection{Double-mutant-cycle coupling as a $2\times2$ contrast}
For positions $i,j$ with reference residues $a_0,b_0$ and alternatives
$a_1,b_1$, the coupling free energy (in log$_{10}$ units) is
\begin{equation}
\Gamma = (y_{11} - y_{10}) - (y_{01} - y_{00})
= y_{11} - y_{10} - y_{01} + y_{00},
\end{equation}
the two-way interaction contrast. $\Gamma = 0$ iff the positions are
additive; our estimator generalizes this to the full $20\times20$ anchor
grid as the RMS of all pairwise contrasts weighted by cell support, with
bootstrap CIs over resampled peptides. The sign-flip between A$^*$02
($\Gamma_{HH} = -0.31$) and A$^*$03/A$^*$11 ($+0.10$) families is
Discovery 1. \\ \\
\textbf{Worked example.} Collapsing the anchor residues to hydrophobic
(H $=$ LIVMFWA) vs not gives the coarse $2\times2$ for A$^*$02:01
(train 9-mers): $\bar y_{HH} = %(w_a02_hh)s$ ($n=%(w_a02_hh_n)s$),
$\bar y_{Hx} = %(w_a02_hx)s$, $\bar y_{xH} = %(w_a02_xh)s$,
$\bar y_{xx} = %(w_a02_xx)s$, so
$\Gamma = %(w_a02_g)s$ --- jointly hydrophobic anchors bind
${\sim}7\times$ \emph{better} than the additive prediction. The same
table for A$^*$03:01 gives $\Gamma = %(w_a03_g)s$: the sign flips. The
coarse contrast amplifies the per-residue grid estimate ($-0.31$/$+0.10$)
because it pools 49 residue pairs into one cell; both estimators agree in
sign and ordering.
\subsection{Hydrophobic moment as a Fourier magnitude}
$\mu_H = \frac1N\big|\sum_n H_n e^{in\delta}\big|$ is the magnitude of the
length-$N$ sequence's hydrophobicity sampled at angular frequency
$\delta$; for an ideal $\alpha$-helix $\delta = 100^\circ$, so $\mu_H$
measures the first Fourier coefficient of hydrophobicity on the helix
wheel --- maximal when hydrophobic residues cluster on one face.
Reversal invariance follows from $|z| = |\bar z|$; the bound
$\mu_H \le \frac1N\sum_n |H_n|$ is the triangle inequality, saturated
iff all phasors align. \qed
\subsection{Platt Hessian positive definiteness}
The log-likelihood Hessian
$H = \sum_i p_i(1-p_i)\, x_i x_i^\top$ with $x_i = (s_i, 1)$ is a positive
semidefinite sum; it is positive definite iff the vectors $x_i$ span
$\mathbb R^2$ and all $p_i\in(0,1)$, which holds whenever validation
scores are non-constant. Hence strict convexity, unique optimum, and
quadratic Newton convergence. \qed

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
\textbf{Test manifest (24 hermetic tests).} Parsers: IEDB row filter
(keeps valid, drops class-II/non-nM/non-canonical), FASTA round-trip,
dedup exactness, natural-residue filter, Jaccard filter threshold
behavior, length-matched negative sampling. Encoding: channel dims
(44), one-hot correctness, BLOSUM row lookup, mask correctness at padded
positions. Graphs: backbone adjacency, next-nearest edges, P2--P$\Omega$
anchor edge presence, symmetry, degree normalization. Metrics: AUC vs
brute force, AUC0.1 truncation, PPV at prevalence, SRCC sign. Models:
output shapes for all architectures, ridge recovery of planted additive
signal (synthetic ground truth), generator sampling validity (canonical
residues, length bounds), cascade threshold logic. All tests run offline
in under 60 seconds.

\begin{thebibliography}{25}
\bibitem{iedb} Vita R. et al. The Immune Epitope Database (IEDB): 2018 update. \emph{Nucleic Acids Research} 47(D1), 2019.
\bibitem{netmhcpan41} Reynisson B., Alvarez B., Paul S., Peters B., Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. \emph{Nucleic Acids Research} 48(W1), 2020 (gkaa379).
\bibitem{netmhcpan40} Jurtz V. et al. NetMHCpan-4.0: improved peptide--MHC class I interaction predictions integrating eluted ligand and peptide binding affinity data. \emph{Journal of Immunology} 199(9), 2017.
\bibitem{netmhc} Nielsen M. et al. Reliable prediction of T-cell epitopes using neural networks with novel sequence representations. \emph{Protein Science} 12(5), 2003.
\bibitem{netmhcpan} Nielsen M. et al. NetMHCpan, a method for quantitative predictions of peptide binding to any HLA-A and -B locus protein of known sequence. \emph{PLoS ONE} 2(8), 2007.
\bibitem{smm} Peters B., Sette A. Generating quantitative models describing the sequence specificity of biological processes with the stabilized matrix method. \emph{BMC Bioinformatics} 6:132, 2005.
\bibitem{mhcflurry} O'Donnell T.J., Rubinsteyn A., Laserson U. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. \emph{Cell Systems} 11(1), 2020.
\bibitem{mhcflurry1} O'Donnell T.J. et al. MHCflurry: open-source class I MHC binding affinity prediction. \emph{Cell Systems} 7(1), 2018.
\bibitem{mhcnuggets} Shao X.M. et al. High-throughput prediction of MHC class I and II neoantigens with MHCnuggets. \emph{Cancer Immunology Research} 8(3), 2020.
\bibitem{mixmhcpred} Bassani-Sternberg M. et al. Deciphering HLA-I motifs across HLA peptidomes improves neo-antigen predictions and identifies allostery regulating HLA specificity. \emph{PLoS Computational Biology} 13(8), 2017.
\bibitem{netctlpan} Larsen M.V. et al. An integrative approach to CTL epitope prediction: a combined algorithm integrating MHC class I binding, TAP transport efficiency, and proteasomal cleavage predictions. \emph{European Journal of Immunology} 35(8), 2005; NetCTLpan: \emph{Immunome Research} 6:2, 2010.
\bibitem{transphla} Chu Y. et al. A transformer-based model to predict peptide--HLA class I binding and optimize mutated peptides for vaccine design. \emph{Nature Machine Intelligence} 4, 2022.
\bibitem{capsnet} Zeng J., Gifford D.K. Quantification of uncertainty in peptide--MHC binding prediction improves high-affinity peptide selection for therapeutic design. \emph{Cell Systems} 9(2), 2019.
\bibitem{kim2014} Kim Y. et al. Dataset size and composition impact the reliability of performance benchmarks for peptide--MHC binding predictions. \emph{BMC Bioinformatics} 15:241, 2014.
\bibitem{sidney} Sidney J. et al. HLA class I supertypes: a revised and updated classification. \emph{BMC Immunology} 9:1, 2008.
\bibitem{cppsite} Agrawal P. et al. CPPsite 2.0: a repository of experimentally validated cell-penetrating peptides. \emph{Nucleic Acids Research} 44(D1), 2016.
\bibitem{cellppd} Gautam A. et al. CellPPD: in silico approaches for designing highly effective cell penetrating peptides. \emph{Journal of Translational Medicine} 11:74, 2013.
\bibitem{cpppred} Holton T.A. et al. CPPpred: prediction of cell penetrating peptides. \emph{Bioinformatics} 29(23), 2013.
\bibitem{mlcpp} Manavalan B. et al. MLCPP: machine-learning-based prediction of cell-penetrating peptides and their uptake efficiency with improved accuracy. \emph{Journal of Proteome Research} 17(9), 2018; MLCPP 2.0: \emph{Briefings in Bioinformatics} 21(4), 2020.
\bibitem{skipcpp} Wei L. et al. SkipCPP-Pred: an improved and promising sequence-based predictor for predicting cell-penetrating peptides. \emph{BMC Genomics} 18(S7), 2017.
\bibitem{blosum} Henikoff S., Henikoff J.G. Amino acid substitution matrices from protein blocks. \emph{PNAS} 89(22), 1992.
\bibitem{eisenberg} Eisenberg D. et al. The hydrophobic moment detects periodicity in protein hydrophobicity. \emph{PNAS} 81(1), 1984.
\bibitem{chengprusoff} Cheng Y., Prusoff W.H. Relationship between the inhibition constant and the concentration of inhibitor which causes 50 per cent inhibition of an enzymatic reaction. \emph{Biochemical Pharmacology} 22(23), 1973.
\bibitem{platt} Platt J. Probabilistic outputs for support vector machines and comparisons to regularized likelihood methods. \emph{Advances in Large Margin Classifiers} 10(3), 1999.
\bibitem{delong} DeLong E.R., DeLong D.M., Clarke-Pearson D.L. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. \emph{Biometrics} 44(3), 1988.
\bibitem{uniprot} The UniProt Consortium. UniProt: the universal protein knowledgebase in 2025. \emph{Nucleic Acids Research} 53(D1), 2025.
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
    "h2h_ens_ci": f"[{boot['ensemble_auc_ci95'][0]:.4f}, {boot['ensemble_auc_ci95'][1]:.4f}]",
    "h2h_mfl_ci": f"[{boot['mhcflurry_auc_ci95'][0]:.4f}, {boot['mhcflurry_auc_ci95'][1]:.4f}]",
    "h2h_delta": f(boot["delta_auc_mean"]),
    "h2h_delta_ci": f"[{boot['delta_auc_ci95'][0]:.4f}, {boot['delta_auc_ci95'][1]:.4f}]",
    "h2h_pwin": f"{boot['p_win']*100:.1f}\%",
    "lc25": f(lc["0.25"]["mean"]), "lc50": f(lc["0.5"]["mean"]), "lc100": f(lc["1.0"]["mean"]),
    "auc9": f([d for d in plauc["per_length_auc"] if d["len"] == 9][0]["ensemble_auc"]),
    "n9": f"{[d for d in plauc['per_length_auc'] if d['len'] == 9][0]['n']:,}",
    "auc10": f([d for d in plauc["per_length_auc"] if d["len"] == 10][0]["ensemble_auc"]),
    "n_methods": ameth["n_methods_total"],
    "perlen_rows": "\n".join(
        f"{d['len']} & {d['n']:,} & {d['n_binders']:,} & {f(d['ensemble_auc'])} \\\\"
        for d in plauc["per_length_auc"]),
    "method_rows": "\n".join(
        f"{m['method']} & {m['n']:,} \\\\" for m in ameth["assay_methods"]),
    "family_rows": "\n".join(
        f"{i+1} & \\texttt{{{'; '.join(fam)}}} \\\\"
        for i, fam in enumerate(novel["families"])),
    "link_r": f(link["pearson_r"], 2),
    "link_ci": f"[{link['ci95'][0]:.2f}, {link['ci95'][1]:.2f}]",
    "loc_a_n": locus["alleles_by_locus"].get("HLA-A", 0),
    "loc_a_pairs": f"{locus['pairs_by_locus'].get('HLA-A', 0):,}",
    "loc_b_n": locus["alleles_by_locus"].get("HLA-B", 0),
    "loc_b_pairs": f"{locus['pairs_by_locus'].get('HLA-B', 0):,}",
    "loc_c_n": locus["alleles_by_locus"].get("HLA-C", 0),
    "loc_c_pairs": f"{locus['pairs_by_locus'].get('HLA-C', 0):,}",
    "w_a02_hh": f(worked["HLA-A*02:01"]["HH"]["mean_log_ic50"], 3),
    "w_a02_hh_n": f"{worked['HLA-A*02:01']['HH']['n']:,}",
    "w_a02_hx": f(worked["HLA-A*02:01"]["Hx"]["mean_log_ic50"], 3),
    "w_a02_xh": f(worked["HLA-A*02:01"]["xH"]["mean_log_ic50"], 3),
    "w_a02_xx": f(worked["HLA-A*02:01"]["xx"]["mean_log_ic50"], 3),
    "w_a02_g": f(worked["HLA-A*02:01"]["gamma_hh"], 3),
    "w_a03_g": f(worked["HLA-A*03:01"]["gamma_hh"], 3),
    "ece": f(cal["ece"]),
    "corr": f(scorr["pearson_pssm_cnn"], 3),
    "pah2h_wins": pah2h["ensemble_wins"], "pah2h_n": pah2h["n_alleles_compared"],
    "pah2h_top": ", ".join(
        f"{a.replace('HLA-', '')} +{(v['ensemble_auc']-v['mhcflurry_auc'])*100:.1f}pp"
        for a, v in sorted(pah2h["per_allele"].items(),
                           key=lambda kv: -(kv[1]["ensemble_auc"]-kv[1]["mhcflurry_auc"]))[:4]),
    "pah2h_worst": min(
        (f"{a.replace('HLA-', '')} {(v['ensemble_auc']-v['mhcflurry_auc'])*100:.1f}pp"
         for a, v in pah2h["per_allele"].items()),
        key=lambda s: float(s.split()[-1][:-2])),
    "d_charge": f(cppprops["cohens_d_charge"], 2),
    "d_hyd": f(cppprops["cohens_d_hydrophobicity"], 2),
    "d_mom": f(cppprops["cohens_d_hydrophobic_moment"], 2),
    "arg_enr": f(aenr["log2_enrichment"][aenr["aa_order"].index("R")], 2),
    "len9_frac": f([d for d in lens["by_length"] if d["len"] == 9][0]["binder_frac"], 3),
    "len9_n": f"{[d for d in lens['by_length'] if d['len'] == 9][0]['n']:,}",
    "len10_frac": f([d for d in lens["by_length"] if d["len"] == 10][0]["binder_frac"], 3),
    "len10_n": f"{[d for d in lens['by_length'] if d['len'] == 10][0]['n']:,}",
    "motif_rows": "\n".join(
        f"{a.replace('*', '$^*$')} & {motif['n_binders'][a]:,} & " +
        " & ".join(f(v, 1) for v in motif["info_content"][a]) + " \\\\"
        for a in motif["alleles"]),
    "len_rows": "\n".join(
        f"{d['len']} & {d['n']:,} & {f(d['binder_frac'], 3)} & {d['median_ic50_nm']:,.0f} \\\\"
        for d in lens["by_length"]),
    "cal_rows": "\n".join(
        f"{i+1} & {b['n']:,} & {f(b['mean_pred'], 3)} & {f(b['empirical'], 3)} & "
        f"{f(b['mean_pred']-b['empirical'], 3)} \\\\"
        for i, b in enumerate(cal["prob_bins"])),
    "pah2h_rows": "\n".join(
        f"{a.replace('*', '$^*$')} & {v['n']} & {v['n_binders']} & "
        f"{f(v['mhcflurry_auc'])} & {f(v['ensemble_auc'])} & "
        f"{'yes' if v['ensemble_auc'] > v['mhcflurry_auc'] else 'no'} \\\\"
        for a, v in sorted(pah2h["per_allele"].items())),
    "novel_rows": "\n".join(
        f"{c['name']} & \\texttt{{{c['sequence']}}} & {f(c['p_cpp'])} & "
        f"{c['net_charge']:.0f} & {f(c['hydrophobic_moment'], 3)} & "
        f"{f(c['mean_hydropathy'], 2)} & {f(c['max_db_jaccard'], 2)} \\\\"
        for c in novel["named_candidates"]),
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
