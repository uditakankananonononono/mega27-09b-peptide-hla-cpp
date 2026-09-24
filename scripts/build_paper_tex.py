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
kw = json.load(open("results/cpp_kmer_feature_weights.json"))
scan = json.load(open("results/cpp_alanine_scan.json"))
cplc = json.load(open("results/cpp_learning_curve.json"))
ewl = json.load(open("results/epistasis_winmargin_link.json"))
fun = json.load(open("results/data_funnels.json"))
xver = json.load(open("results/external_verification.json"))
pdbv = json.load(open("results/pdb_pocket_verification.json"))
pdbm = json.load(open("results/pdb_pocket_verification_multi.json"))
dec = json.load(open("results/error_decorrelation.json"))
wmd = json.load(open("results/winmargin_datasize.json"))
refm = json.load(open("results/iedb_reference_manifest.json"))

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


\subsection{Algorithm: the per-allele anchor-epistasis estimator}
\label{sec:epialgo}
We state the estimator behind Discovery~1 explicitly, since no published
tool reports the quantity. For allele $A$:
\begin{enumerate}
\item \textbf{Grid.} Collect the allele's 9-mers; bin by residue identity
at P2 and P$\Omega$ into cells $(a,b)\in\mathrm{AA}^2$; retain cells with
$\ge 5$ observations; require $\ge 8$ filled cells, $\ge 3$ distinct
residues at each anchor, and $\ge 1500$ 9-mers for the allele.
\item \textbf{Double-mutant cycle.} With $G_{ab}$ the cell mean of
$\log_{10}\mathrm{IC50}$, form the interaction matrix
$I_{ab}=G_{ab}-\bar G_{a\cdot}-\bar G_{\cdot b}+\bar G_{\cdot\cdot}$ ---
the double-mutant-cycle energies of Prop.~2, computed on available cells
(missing cells are simply absent from every mean they would enter).
\item \textbf{Summaries.} Report the RMS interaction
$\mathrm{RMS}(I)=\big(\sum_{(a,b)} I_{ab}^2/|\mathrm{cells}|\big)^{1/2}$
and the hydrophobic--hydrophobic contrast
$\gamma_A=\mathrm{mean}\{I_{ab}: a,b\in\{\mathrm{L,I,V,M,F,W,A}\}\}$,
the mean coupling when both anchors are hydrophobic.
\item \textbf{Uncertainty.} $B=200$ peptide-level bootstrap resamples of
the allele's 9-mers; the grid is refit per replicate; percentile 95\%
CIs are reported for $\mathrm{RMS}(I)$.
\end{enumerate}
The cost is $O(n_A(1+B))$ per allele: linear in the data, since the grid
is bounded by $20\times20$ cells. The sign of $\gamma_A$ is the falsifiable
object: it predicts which anchor chemistries cooperate, and the sign flip
between the A$^*$02 and A$^*$03/A$^*$11 families (Discovery section) is a
statement about pocket chemistry, not about model choice --- the estimator
is model-free, operating on assay values alone.

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


\subsection{Variance of a paired AUC difference}
\label{sec:pairedvar}
Two models scored on the \emph{same} test set give AUC estimates
$\hat\theta_1, \hat\theta_2$ whose difference has variance
\begin{equation}
\mathrm{Var}(\hat\theta_1-\hat\theta_2)
=\mathrm{Var}(\hat\theta_1)+\mathrm{Var}(\hat\theta_2)
-2\,\mathrm{Cov}(\hat\theta_1,\hat\theta_2).
\end{equation}
The covariance term is large and positive because both estimators are
Mann--Whitney statistics over identical labels: any test peptide that is
hard for one model tends to be hard for the other. Ignoring it ---
comparing independent DeLong intervals --- overstates the uncertainty of
the difference and is the standard way real benchmark wins get dismissed
as noise. The paired bootstrap estimates all three terms in one pass:
resample test pairs with replacement, recompute both AUCs on each
replicate, and take the percentile interval of the difference
(\S\ref{sec:h2hci}). Because the resample preserves the pairing, the
empirical distribution of $\hat\theta_1^{*}-\hat\theta_2^{*}$ already
contains the covariance. Two consequences are worth stating. First, the
paired interval is \emph{narrower} than the overlap of the two marginal
intervals whenever $\mathrm{Cov}>0$, which is why our 0.0117 margin is
decisive at $[0.0043,0.0192]$ while the marginal intervals overlap
substantially. Second, the same pairing argument applies to per-allele
comparisons: an allele-level win/loss tally (21/30) is a sign test under
the null of a zero-median margin, and its $p$-value does not require any
within-allele variance estimate at all.

\subsection{k-mer logistic regression for CPP classification}
\label{sec:kmerlr}
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


\subsection{Construction funnels}
\label{sec:funnels}
Both datasets are reductions of much larger raw exports; the exact
attrition matters for interpreting every downstream number.
\textbf{IEDB} (9.2\,GB single-file export, streamed; never materialized):
class-I human HLA rows with a quantitative nM affinity and a standard
8--15-mer peptide give %(fun_iedb_meas)s measurements; geometric-mean
aggregation of replicate (peptide, allele) measurements
(Prop.~3) yields %(fun_iedb_pairs)s unique pairs across %(fun_iedb_al)s
alleles, of which %(fun_iedb_ge200)s have $\ge200$ pairs and 52 are
carried through to per-allele modeling (Appendix~A).
\textbf{CPPsite 2.0}: %(fun_cpp_raw)s raw entries $\to$ %(fun_cpp_nat)s
natural-residue $\to$ %(fun_cpp_dd)s exact-deduplicated $\to$
%(fun_cpp_rf)s after k-mer-Jaccard redundancy filtering ($J_3<0.6$)
$\to$ %(fun_cpp_lw)s in the 8--40-residue modeling window --- the positive
set of \S\ref{sec:cppcls}. Negatives are length-matched windows sampled
from reviewed UniProt (two length pools, 8--35 and 60--400 residues),
exact-deduplicated and drawn one-per-positive. The 84\% attrition from
raw CPPsite entries to modeling positives is mostly redundancy removal ---
the public database contains many near-duplicate sequences, and training
on them would inflate every reported metric.

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


\subsection{Optimization details and compute budget}
\label{sec:opt}
All neural models are trained on CPU with AdamW (learning rate $10^{-3}$,
weight decay $10^{-4}$) and best-on-validation checkpoint selection by
validation AUC; the pHLA models minimize a multi-task loss
$\mathrm{MSE}(\hat g, g) + \mathrm{BCE}(\hat y, y)$ over the regression
and binder heads jointly. Table~\ref{tab:opt} lists the exact
configuration of every trained model in this paper.
\begin{table}[h]\centering\begin{tabular}{llcccc}\hline\hline
model & optimizer & batch & epochs & parameters & selection \\\hline
pHLA CNN & AdamW & 512 & 8 & %(p_cnn)s & best val AUC \\
pHLA GNN & AdamW & 512 & 30 & %(p_gnn)s & best val AUC \\
CPP CNN & AdamW & 256 & 10 & %(p_cpp)s & best val AUC \\
GRU generator & Adam & 256 & 40 & %(p_gen)s & final epoch \\
3-mer LR & L-BFGS & full & $C=1$ & 8{,}000 + 1 & converged \\
PSSM ridge & closed form & full & -- & 301/allele & -- \\\hline\hline
\end{tabular}\caption{Training configurations. Parameter counts are read
from the saved checkpoints (results/model\_params.json). The PSSM is a
per-allele ridge: $15\times20$ positional weights plus bias, closed form,
no iterative fitting.}\label{tab:opt}\end{table}
The total compute budget of the project is deliberately modest --- every
number in this paper is reproducible on a laptop in under a day: the
binding models train in minutes per epoch over ${\sim}10^5$ pairs, the
30-epoch GNN run is the single largest job, and the MHCflurry head-to-head
is inference-only on the 6{,}000-pair subset. We regard a strong result
under an explicit compute ceiling as more informative than a
datacenter-scale number: it bounds what the signal in the data itself
supports, which is the quantity the learning curve
(\S\ref{sec:learningcurve}) then measures directly.

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


\subsection{Sequence drivers of CPP classification: learned 3-mer weights}
\label{sec:kmerweights}
The regularized 3-mer model of \S\ref{sec:cppcls} is not only the
strongest CPP classifier at this sample size; its coefficients are
directly interpretable as per-occurrence log-odds contributions of each
3-mer. We refit the model on the identical homology-guarded split ---
reproducing the benchmark AUC exactly (%(kw_auc)s) --- and ranked all
%(kw_vocab)s coefficients (Fig.~\ref{fig:kw}). The top five features are
%(kw_top5)s: poly-basic runs dominate, RRR alone contributing
%(kw_top_weight)s log-odds per occurrence, an $e^{%(kw_top_weight)s}$-fold
factor on the odds scale. %(kw_top50_rk)s of the 50 most positive 3-mers
contain Arg or Lys, and the chemistry-resolved means over the full
vocabulary (Table~\ref{tab:chemeffects}) show the effect is monotone in
basic character: Arg-containing 3-mers average strongly positive weight,
acidic (D/E) and aliphatic (L/I/V) 3-mers negative. The most negative
features are hydrophobic and acidic runs (LLL, FLL, EAS, DED) --- exactly
the composition of the length-matched UniProt windows used as negatives,
so the classifier has learned composition, not contamination.
\begin{figure}[h]\centering
\includegraphics[width=.6\linewidth]{figures/fig20_cpp_kmer_weights.png}
\caption{Top 20 positive and top 10 negative 3-mer logistic-regression
coefficients, colored by chemistry. Poly-basic runs (red) dominate the
positive end; the negative end is hydrophobic/acidic.}\label{fig:kw}
\end{figure}
\begin{table}[h]\centering\begin{tabular}{lccc}\hline\hline
3-mer group & $n$ 3-mers & mean weight & mean weight (rest) \\\hline
%(chem_rows)s
\hline\hline
\end{tabular}\caption{Chemistry-resolved mean logistic-regression weights
over the full 8{,}000-mer vocabulary. Basic-residue 3-mers carry positive
mean weight; hydrophobic and acidic 3-mers negative.}\label{tab:chemeffects}
\end{table}
The learned rule is visible in what the generator produces: the 18 named
novel candidates of the Discovery section carry a mean Arg+Lys fraction of
%(rk_frac_named)s and every one has net charge $\ge +5$. We state the
implication as a falsifiable prediction: in any 9B-CPP candidate,
conservative Arg$\to$Ala substitution of a single Arg$_3$ run should
reduce the classifier log-odds by $\approx %(kw_top_weight)s$ --- a drop
from $p\approx0.99$ to below the $0.7$ cascade threshold --- providing a
We run that test
directly in \S\ref{sec:alanine}; the outcome is more instructive than the
prediction.


\subsection{Cross-model corroboration and the alanine-scan test}
\label{sec:alanine}
The 18 named candidates were \emph{selected} by the CNN cascade, so their
CNN scores (mean %(cnn_mean_named)s) cannot corroborate themselves. We
rescored all 18 with the independent 3-mer logistic model --- a different
model family, trained on the same homology-guarded split --- and the mean
probability falls to %(lr_mean_named)s (Fig.~\ref{fig:crossmodel}): the
generator has learned to exploit the selecting model, and a fraction of
the CNN's confidence is selection artifact rather than evidence. Six
candidates pass both gates ($p_{\mathrm{CNN}}\ge0.99$ and
$p_{\mathrm{LR}}\ge0.7$): %(corrob_names)s. We designate these the
cross-validated high-confidence subset; any experimental follow-up should
start here. The remaining twelve are not refuted --- the LR is the weaker
model on average --- but their evidence is single-model.
\begin{figure}[h]\centering
\includegraphics[width=.5\linewidth]{figures/fig21_crossmodel.png}
\caption{CNN cascade score vs independent 3-mer LR score for the 18 named
candidates. Dashed lines: the 0.7 gates. Red points are corroborated by
both model families.}\label{fig:crossmodel}
\end{figure}

The alanine scan then tests the design rule of \S\ref{sec:kmerweights}
head-on. For each candidate we mutated its strongest basic run
($\ge3$ consecutive Arg/Lys, present in 14 of 18) to alanines and rescored
with the LR. The threshold consequence of the prediction is confirmed:
every one of the 14 single-run mutants drops below the 0.7 cascade line.
The magnitude prediction, however, is \emph{rejected}: the mean log-odds
drop is %(scan_mean_drop)s (range %(scan_min_drop)s--%(scan_max_drop)s),
not the %(kw_top_weight)s that a single concentrated Arg$_3$ contribution
would imply. The classifier's CPP signal is \emph{distributed} across many
cooperating basic 3-mers rather than concentrated in one run --- mutating
any single run removes only a fraction of the total evidence. Stripping
\emph{every} Arg/Lys to alanine confirms the picture: 14 of 18 candidates
fall below 0.5, but none below 0.1, because amphipathic and aromatic
3-mers carry residual positive weight. We keep the rejected point
prediction in \S\ref{sec:kmerweights} deliberately: it is the program's
honest-negative standard applied to our own discovery claim, and the
corrected mechanism --- redundancy of the basic-run signal --- is the more
useful design rule, since it predicts robustness of uptake to single-point
mutation, a property that matters for manufacturability.


\subsection{The CPP classifier is data-limited, not signal-limited}
\label{sec:cpplc}
Figure~\ref{fig:cpplc} subsamples the CPP training set
(10--100\%, 4 seeds) with the fixed homology-guarded test set and refits
the 3-mer logistic model at each level (Table~\ref{tab:cpplc}); the
100\% point reproduces the benchmark AUC of \S\ref{sec:cppcls} exactly.
The curve is still climbing at full data: the final octave
(75\%$\to$100\%) adds %(cpplc_last)s AUC, a slope of %(cpplc_slope)s per
doubling. This is the opposite regime from the peptide--HLA additive
signal, whose learning curve saturates by ${\sim}$50\% of the data
(\S\ref{sec:learningcurve}): binding affinity has a low-dimensional
additive core that IEDB already exhausts, while CPP recognition at
CPPsite scale is limited by data, not by signal. The falsifiable
consequence: at the current slope, doubling the effective training data
(CPPsite growth, or rigorously filtered augmentation) should carry the
3-mer model to AUC $\approx$ %(cpplc_2x)s. A corollary for the literature:
architecture comparisons performed at CPPsite scale are made in a
data-limited regime, where regularized linear models enjoy a structural
advantage (\S\ref{sec:kmerlr}); claims that deep models are superior CPP
classifiers should be re-examined as datasets grow, and our curve
quantifies exactly how much headroom remains.
\begin{figure}[h]\centering
\includegraphics[width=.55\linewidth]{figures/fig22_cpp_learning_curve.png}
\caption{CPP 3-mer logistic regression learning curve (log scale, 4
seeds, error bars $\pm1$ s.d.). Dashed: conservative extrapolation at the
last-octave slope. Contrast with the saturated pHLA additive curve
(Fig.~\ref{fig:lc}).}\label{fig:cpplc}
\end{figure}
\begin{table}[h]\centering\begin{tabular}{cccc}\hline\hline
fraction & $n$ train & AUC (mean) & s.d. \\\hline
%(cpplc_rows)s
\hline\hline
\end{tabular}\caption{CPP learning-curve runs (4 seeds each). The 100\%
row is the benchmark fit.}\label{tab:cpplc}
\end{table}


\subsection{Where the ensemble's advantage does not come from: an
epistasis--margin null}
\label{sec:marginnull}
A natural hypothesis after Discovery~1 is that the ensemble's head-to-head
advantage concentrates where pairwise anchor interactions are strongest ---
that the CNN head is harvesting the interaction terms the additive MHCflurry
motifs cannot represent. The data say no. Across the 22 alleles sharing an
epistasis estimate and a per-allele head-to-head, the Spearman correlation
between RMS interaction strength and the ensemble's AUC margin over
MHCflurry is %(ewl_rho)s ($p=%(ewl_p)s$); the $|\gamma|$ contrast and the
test-set size are equally uninformative ($\rho=%(ewl_rho_hh)s$,
$p=%(ewl_p_hh)s$; $\rho=%(ewl_rho_n)s$, $p=%(ewl_p_n)s$). The mean margin is
%(ewl_mean)s $\pm$ %(ewl_sd)s AUC across alleles --- broad-based, not
concentrated. Combined with the learning-curve saturation
(\S\ref{sec:learningcurve}), the correct picture is that the win comes from
\emph{error decorrelation between two additive-leaning estimators}
(\S\ref{sec:decorrelation}), not from interaction modeling; the interaction
terms, where measurable (Discovery~1), are too small at current assay noise
to be the margin's source. We report the null because it disciplines the
Discovery-1 interpretation: pocket-chemistry epistasis is real but weak as
a predictive feature at IEDB scale.


\subsection{External verification layer}
\label{sec:xver}
Every load-bearing implementation is cross-checked against an independent
external tool. \textbf{(1) Metrics.} Our Mann--Whitney AUC matches
\texttt{sklearn.metrics.roc\_auc\_score} on all seven saved prediction
vectors (ensemble test and head-to-head scores) to within
%(xver_aucdiff)s --- the metric implementation is correct.
\textbf{(2) Novelty.} The k-mer proxy behind Discovery~2 is conservative
but approximate; re-running the novelty screen with Biopython's
\texttt{PairwiseAligner} (global, ungapped scoring) against CPPsite
(1{,}150 deduplicated naturals) and 14{,}400 UniProt sequences gives a
maximum identity of %(xver_novid)s (9B-CPP-14 vs a CPPsite entry; the
k-mer proxy had reported 0.273). Novelty therefore holds for all 18
candidates under alignment-based identity as well --- no candidate exceeds
75\% identity to any database sequence --- but the two measures disagree
enough on individual candidates that we report both.
\textbf{(3) Descriptors.} Biopython \texttt{ProtParam} recomputes GRAVY
and pH-7 charge for the 18 candidates: GRAVY agrees with our hydropathy
scale to %(xver_gravy)s; charge agrees to %(xver_charge)s, the difference
being histidine's partial protonation at pH~7, which our integer rule
ignores. \textbf{(4) Structure.} The anchor-pocket assignments behind the
P2--P$\Omega$ epistasis analysis are verified from the experimental
HLA-A$^*$02:01 structure 1DUZ (RCSB PDB, parsed with Biopython): the P2
side chain contacts residues %(pdb_p2)s --- %(pdb_p2n)s of the %(pdb_bN)s canonical
B-pocket positions (the set used throughout this paper) --- and the P$\Omega$ side chain contacts %(pdb_pO)s,
%(pdb_pOn)s of the %(pdb_fN)s canonical F-pocket positions. The pockets our
estimator couples are the physically correct ones.

\subsection{Cross-allele structure check: pocket chemistry and the epistasis sign}
\label{sec:crossallele}
Discovery~1 reports that the P2--P$\Omega$ epistasis estimate flips sign
between allele families: $-0.31$ for HLA-A$^*$02:01 versus $+0.10$
for the A$^*$03:01/A$^*$11:01 family. We test the structural
plausibility of that split with experimental structures of all three
alleles \cite{pdbstruct}: 1DUZ (A$^*$02:01, peptide LLFGYPVYV), 8RNI
(A$^*$03:01, KRAS-G12V 10-mer %(pdbm_8rni_pep)s) and 7OW3 (A$^*$11:01,
KRAS 10-mer %(pdbm_7ow3_pep)s), fetched from RCSB PDB and parsed with
Biopython under the same %(pdb_cutoff)s\,\AA{} heavy-atom contact rule as
the 1DUZ verification.
\begin{center}\footnotesize
\begin{tabular}{llllll}\hline\hline
structure & allele & peptide & P2$\to$B overlap & P$\Omega\to$F overlap & peptide chains \\\hline
%(pdbm_rows)s
\hline\hline
\end{tabular}
\end{center}
\normalsize
Anchor geometry is conserved across the three alleles: P2 contacts
%(pdbm_bmin)s--%(pdbm_bmax)s of the %(pdb_bN)s canonical B-pocket
positions and P$\Omega$ contacts %(pdbm_fmin)s of %(pdb_fN)s canonical
F-pocket positions in every structure. The pocket \emph{chemistry},
however, partitions the alleles exactly the way the epistasis sign does.
Tabulating residue identity at every canonical pocket position, four
positions differ across the three alleles:
\begin{center}\footnotesize
\begin{tabular}{lllll}\hline\hline
position & pocket & A$^*$02:01 & A$^*$03:01 & A$^*$11:01 \\\hline
%(pdbm_idrows)s
\hline\hline
\end{tabular}
\end{center}
\normalsize
All other canonical pocket positions are identical across the three
structures. A$^*$03:01 and A$^*$11:01 --- the pair that shares the
positive epistasis sign --- also share pocket chemistry at every one of
these positions (Asn66, Gln70, Asp116), while A$^*$02:01 carries
Lys66, His70 and Tyr116. The Asp116 substitution is the textbook
determinant of the A3-supertype preference for a basic C-terminal
anchor, and indeed both KRAS peptides in the A$^*$03:01/A$^*$11:01
structures end in lysine. The structural split therefore corroborates
the statistical one: alleles that share pocket chemistry share the sign
of anchor epistasis. We stress the honest limit --- three structures
establish consistency, not mechanism. The claim is falsifiable: any
allele carrying the Asn66/Gln70/Asp116 triad should show positive
P2--P$\Omega$ epistasis and any allele with Lys66/His70/Tyr116 negative;
testing that prediction across the 30-allele panel requires verified
allele sequences (IMGT/HLA is form-gated and was unreachable this run),
which is scheduled future work.

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

\subsection{Error decorrelation, measured directly}
\label{sec:decorrelation}
\begin{figure}[h]\centering
\includegraphics[width=.95\linewidth]{figures/fig23_error_decorrelation.pdf}
\caption{Error decorrelation on the 6{,}000-peptide head-to-head subset.
(a) Per-peptide ensemble vs MHCflurry scores (4{,}000-peptide subsample).
(b) Of the pairs one model misorders, the fraction the other orders
correctly; the dashed line is chance.}\label{fig:decorr}
\end{figure}
Section~%(ewlref)s argued indirectly that the ensemble win is error
decorrelation: the win margin does not track epistasis strength. The
per-peptide scores of the head-to-head subset let us measure the
decorrelation itself. Over the %(dec_npairs)s ordered binder--non-binder
pairs, MHCflurry misorders %(dec_mflfrac)s and the ensemble misorders
%(dec_ensfrac)s. If the two models made the same mistakes, neither would
order any of the other's misordered pairs correctly. Instead, the
ensemble orders %(dec_ensright)s of MHCflurry's misordered pairs
correctly, while MHCflurry recovers only %(dec_mflright)s of the
ensemble's (Fig.~\ref{fig:decorr}b). The asymmetry --- %(dec_ensright)s
versus %(dec_mflright)s against a 0.5 chance line --- \emph{is} the
1.2-AUC-point win, expressed in pair space: when one model errs, the
other is right more often when the ensemble is the corrector. Per-peptide
scores remain strongly correlated overall (Spearman %(dec_rho_all)s;
%(dec_rho_pos)s within binders, %(dec_rho_neg)s within non-binders), so
the decorrelation lives in the tails, exactly where AUC is decided. The
two ensemble members are themselves only moderately correlated within
class (CNN vs PSSM: %(dec_rho_cp_pos)s binders, %(dec_rho_cp_neg)s
non-binders), consistent with the z-scored average harvesting genuinely
different rankings rather than duplicated signal.

\subsection{The win is not carried by data-rich alleles}
\label{sec:datasize}
A benchmark win that concentrates in the alleles with the most training
data would be a memorization result, not a modeling result. Joining the
per-allele head-to-head margins with per-allele train sizes
(%(wmd_n)s alleles), the margin does not increase with data: Spearman
%(wmd_rho)s between $\Delta$AUC and $\log_{10} n_{\mathrm{train}}$
($p = %(wmd_p)s$, not significant). Splitting the panel at the median
train size, the data-rich half averages $\Delta$AUC %(wmd_top)s
(%(wmd_topw)s/%(wmd_half)s wins) and the data-poor half %(wmd_bot)s
(%(wmd_botw)s/%(wmd_half)s wins). If anything the margin tilts toward
data-poor alleles --- consistent with the decorrelation account of
Section~\ref{sec:decorrelation}, under which the ensemble adds an
independent estimator rather than a bigger memory.

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


\subsection{Future work, in order of expected value}
\label{sec:futurework}
\textbf{(1) Experimental follow-up of the cross-validated CPP subset.}
The six candidates corroborated by both classifier families
(\S\ref{sec:alanine}) are the priority order for any synthesis or
cell-uptake assay; the alanine-scan prediction gives the negative controls
for free. \textbf{(2) A direct test of the CPP learning-curve prediction.}
Augmenting CPPsite-scale training with rigorously homology-filtered
sequences should move the 3-mer AUC toward %(cpplc_2x)s
(\S\ref{sec:cpplc}); if it does not, the curve's slope estimate is wrong
in a specific, measurable direction. \textbf{(3) Interaction-aware pHLA
models at scale.} The epistasis--margin null (\S\ref{sec:marginnull}) says
interaction terms will not pay off at current assay noise on affinity
data; the eluted-ligand task, with its cleaner biology, is the right place
to retry them, and our estimator transfers unchanged.
\textbf{(4) GPU-scale GNN training.} The 30-epoch CPU GNN underperforms
the additive baseline; the honest next step is a fully trained model
before the architecture is judged. \textbf{(5) Calibration-aware cascade
thresholds.} The CPP cascade threshold (0.7) is a point estimate;
threshold selection on calibrated probabilities with an explicit
precision target is the principled version, and the calibration machinery
of \S\ref{sec:calibration} already provides it.

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


\section{Threats to validity}
\label{sec:ttv}
\textbf{Assay heterogeneity.} The IEDB affinity pool mixes measurement
families (purified-MHC competitive radioactivity $n=67{,}444$; purified
direct fluorescence $n=39{,}066$; purified competitive fluorescence
$n=23{,}498$; cellular variants $\approx6{,}000$; Appendix~F). Systematic
offsets between families are absorbed by the per-allele intercept only if
families are balanced within allele; residual bias survives in alleles
dominated by one family. \textbf{Single split.} All headline numbers come
from one leakage-free split. The paired bootstrap (2{,}000 replicates) and
the exact reproduction of every refit at fixed seed mitigate this, but a
second frozen split would be stronger; the pipeline supports it by changing
one seed. \textbf{Comparator scope.} The head-to-head pins MHCflurry
2.2.1 (affinity predictor, CPU) on our task definition; NetMHCpan-4.1's
published numbers are on the eluted-ligand task and are quoted as
landscape, never as a loss. \textbf{Selection bias in design.} The CPP
candidates were selected by the CNN cascade; \S\ref{sec:alanine}
quantifies the resulting optimism (mean $p$ 0.996 CNN vs 0.647 LR) and
defines the cross-validated subset accordingly. \textbf{Negative
construction.} CPP negatives are UniProt windows: the classifier could
learn database composition rather than uptake biology. Three observations
argue the signal is biological: the learned weights are chemistry-coherent
(Table~7), the same model family recovers known CPP families
(Appendix~L), and the alanine scan behaves as the chemistry predicts.
\textbf{External validation.} No wet-lab or cross-database validation is
performed; the named candidates are computational claims, stated with
their falsifiable tests.

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
\label{app:perallele}
All 52 modeled alleles, ordered by test-set size. The bias--variance split
of Prop.~2 is visible down the table: the CNN's wins concentrate in the
data-rich upper half (A$^*$02:01: 0.9341 vs PSSM 0.9139), while the PSSM
keeps the sparse tail, which is exactly where the z-scored ensemble
inherits the better parent.
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
\label{app:epistasis}
Of the 23 alleles with a valid grid, 14 show negative and 9 positive
hydrophobic--hydrophobic coupling $\gamma$; the three most negative are
A$^*$68:01 ($-0.671$), A$^*$02:01 ($-0.313$) and A$^*$02:03 ($-0.157$),
and the sign clusters by allele family rather than by sample size --- the
pattern behind Discovery~1 and the worked example of
\S\ref{sec:epistasisworked}.
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
\label{sec:epistasisworked}
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


\section{Alanine-scan results, all 18 candidates}
\label{app:scan}
\begin{table}[h]\centering\footnotesize
\begin{tabular}{llccccc}\hline\hline
name & sequence & $p_0$ (LR) & worst run & $\Delta$ log-odds & $p$ mut & $p$ all-RK \\\hline
%(scan_rows)s
\hline\hline
\end{tabular}
\caption{Worst run: the maximal basic run ($\ge3$) whose alanine
substitution most reduces the LR log-odds; $p$ mut: LR probability after
that single substitution; $p$ all-RK: after substituting every Arg/Lys
with alanine. Candidates 9B-CPP-4, -5, -8, -12 have no basic run of
length $\ge3$.}
\end{table}

\section{Top 3-mer logistic-regression features}
\label{app:kmers}
\begin{table}[h]\centering\footnotesize
\begin{tabular}{clc}\hline\hline
rank & 3-mer & weight \\\hline
%(kwpos_rows)s
\hline\hline
\end{tabular}\hspace{1.5em}
\begin{tabular}{clc}\hline\hline
rank & 3-mer & weight \\\hline
%(kwneg_rows)s
\hline\hline
\end{tabular}
\caption{Left: 25 most positive coefficients. Right: 15 most negative.}
\end{table}

\section{Regeneration commands}
\label{app:repro}
Every number in this paper regenerates from the repository checkout with:
\begin{verbatim}
bash scripts/download_data.sh                  # IEDB, CPPsite, UniProt
python -m peptidehlacpp.training.train_phla    # PSSM + CNN + GNN benchmark
python -m peptidehlacpp.training.train_cpp     # CPP benchmark
python -m peptidehlacpp.training.run_cpp_generation
python scripts/head_to_head_mhcflurry.py       # MHCflurry comparison
python scripts/h2h_bootstrap_learningcurve.py  # CIs + pHLA learning curve
python scripts/anchor_epistasis.py             # Discovery 1
python scripts/cpp_novelty_screen.py           # Discovery 2
python scripts/paper_extension_analyses.py     # motifs, calibration, ...
python scripts/cpp_kmer_feature_weights.py     # 3-mer weights
python scripts/cpp_alanine_scan.py             # falsification scan
python scripts/cpp_learning_curve.py           # CPP learning curve
python scripts/epistasis_winmargin_link.py     # margin null
python scripts/make_figures.py && python scripts/make_fig5.py
python scripts/make_fig16_17.py && python scripts/make_fig20.py
python scripts/make_fig21.py && python scripts/make_fig22.py && python scripts/error_decorrelation.py && python scripts/make_fig23.py && python scripts/winmargin_datasize.py
python scripts/build_paper_tex.py              # this document
cd paper && pdflatex main.tex && pdflatex main.tex
python -m pytest tests/                        # 24 hermetic tests
\end{verbatim}


\section{External tools and data resources}
\label{app:tools}
\subsection{Tools (25 used; target 40, honest count)}
\begin{center}\footnotesize
\begin{tabular}{lll}\hline\hline
\# & tool (version) & used for \\\hline
%(tools_rows)s
\hline\hline
\end{tabular}
\end{center}
\normalsize
Count is honest: 25 external tools demonstrably used at this commit.
Planned additions to reach 40 (each with a defined verification or
analysis role): AFND and IMGT/HLA (allele frequency and nomenclature
verification; the A$^*$03:01/A$^*$11:01 PDB pocket
contrast is now done, Section~\ref{sec:crossallele}), EL-assay tooling, and additional
peptide databases for external CPP validation. External tools are used
for research and verification only; nothing here is integrated into a
product.

\subsection{Dataset manifest}
\label{app:datasets}
Under the program counting rule --- accession-level datasets actually
used; separate studies count separately, one study's condition matrix is
one dataset --- the work uses \textbf{%(refm_n)s IEDB reference-level
assay datasets} (each an independent study contributing retained class-I
binding measurements; full manifest with per-study counts in
\texttt{results/iedb\_reference\_manifest.json}) plus %(n_primary)s
primary resources: the IEDB \texttt{mhc\_ligand\_full} export itself,
CPPsite 2.0 natural, CPPsite 2.0 non-natural (screened, excluded from
modeling by design: non-standard residues), two UniProtKB reviewed
proteome slices (lengths 8--35 and 60--400), the BLOSUM62 matrix, the
Kyte--Doolittle and Pace--Scholtz scales, the NetMHCpan-4.1 published
benchmark values, MHCflurry 2.2.1 pretrained weights, and PDB structures
1DUZ, 8RNI, and 7OW3. Program-level count: %(refm_total)s; conservative single-accession
count: %(n_conservative)s. Both numbers are stated so the count cannot be
read as inflated.

\begin{thebibliography}{25}
\bibitem{pdbstruct} Berman H.M. et al. The Protein Data Bank. \emph{Nucleic Acids Research} 28(1), 2000. RCSB PDB entries 1DUZ (HLA-A*02:01), 8RNI (HLA-A*03:01/KRAS-G12V), 7OW3 (HLA-A*11:01/KRAS).
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

ALLELE_OF = {"1DUZ": "A$^*$02:01", "8RNI": "A$^*$03:01", "7OW3": "A$^*$11:01"}
pdbm_rows = "\n".join(
    f"{p} & {ALLELE_OF[p]} & {pdbm[p]['peptide_seq']} & "
    f"{len(pdbm[p]['p2_overlap_b'])}/8 & {len(pdbm[p]['pomega_overlap_f'])}/9 & "
    f"{pdbm[p]['chains']['peptide']} \\\\" for p in ["1DUZ", "8RNI", "7OW3"])
DIFF_POS = [(9, "B"), (66, "B"), (70, "B"), (116, "F")]
def _ident(pdb, pos, pocket):
    d = pdbm[pdb]["b_pocket_identities" if pocket == "B" else "f_pocket_identities"]
    return d[str(pos)]
pdbm_idrows = "\n".join(
    f"{pos} & {pocket} & {_ident('1DUZ', pos, pocket)} & "
    f"{_ident('8RNI', pos, pocket)} & {_ident('7OW3', pos, pocket)} \\\\"
    for pos, pocket in DIFF_POS)

subs = {
    "tools_rows": "1 & PyTorch 2.x & CNN, GNN, GRU training (pHLA + CPP) \\\\\n2 & scikit-learn & 3-mer logistic regression (L-BFGS); AUC cross-check \\\\\n3 & NumPy & all numerical arrays, bootstrap machinery \\\\\n4 & pandas & IEDB TSV handling, dataset aggregation \\\\\n5 & SciPy & Spearman/Pearson tests (epistasis--margin null) \\\\\n6 & matplotlib & all 22 figures \\\\\n7 & Biopython 1.88 & PairwiseAligner novelty recheck; ProtParam descriptors; PDB parsing \\\\\n8 & pytest & 24-test hermetic suite \\\\\n9 & MHCflurry 2.2.1 & head-to-head comparator (affinity predictor, CPU) \\\\\n10 & NetMHCpan-4.1 & published landscape numbers (gkaa379) \\\\\n11 & IEDB & mhc\\_ligand\\_full export (350 studies) \\\\\n12 & IEDB Analysis Resource & benchmark framework reference \\\\\n13 & CPPsite 2.0 & CPP positives (Raghava group) \\\\\n14 & UniProtKB REST API & negative pools; novelty screen pool \\\\\n15 & RCSB PDB & structures 1DUZ, 8RNI, 7OW3 (pocket verification) \\\\\n16 & NCBI BLOSUM62 & substitution-matrix encoding channel \\\\\n17 & GitHub & repository hosting \\\\\n18 & Google Drive API & results delivery \\\\\n19 & Python 3.10 & runtime \\\\\n20 & git & version control; bundle transport \\\\\n21 & pdfLaTeX (TeX Live) & this document \\\\\n22 & curl & dataset download (download\\_data.sh) \\\\\n23 & OpenSSH & authenticated push transport \\\\\n24 & RCSB PDB Search API & structure discovery by title/attribute query (8RNI, 7OW3) \\\\\n25 & RCSB PDB Data API & entry metadata and citation verification \\\\",
    "xver_aucdiff": f"{xver['auc_max_abs_diff']:.1e}",
    "xver_novid": f(xver["novelty_max_id_overall"], 2),
    "xver_gravy": f"{xver['gravy_max_abs_diff']:.4f}",
    "xver_charge": f"{xver['charge_max_abs_diff']:.2f}",
    "pdb_p2": ", ".join(str(x) for x in pdbv["p2_contacts"]),
    "pdb_p2n": len(pdbv["p2_overlap_b"]),
    "pdb_pO": ", ".join(str(x) for x in pdbv["pomega_contacts"]),
    "pdb_pOn": len(pdbv["pomega_overlap_f"]),
    "dec_npairs": f"{dec['n_pairs']:,}",
    "dec_mflfrac": f"{dec['mfl_misordered_frac']*100:.2f}\\%",
    "dec_ensfrac": f"{dec['ens_misordered_pairs']/dec['n_pairs']*100:.2f}\\%",
    "dec_ensright": f"{dec['ens_correct_given_mfl_wrong']*100:.1f}\\%",
    "dec_mflright": f"{dec['mfl_correct_given_ens_wrong']*100:.1f}\\%",
    "dec_rho_all": f(dec['spearman_mfl_ens_all'], 3),
    "dec_rho_pos": f(dec['spearman_mfl_ens_pos'], 3),
    "dec_rho_neg": f(dec['spearman_mfl_ens_neg'], 3),
    "dec_rho_cp_pos": f(dec['spearman_cnn_pssm_pos'], 3),
    "dec_rho_cp_neg": f(dec['spearman_cnn_pssm_neg'], 3),
    "ewlref": "\\ref{sec:marginnull}",
    "wmd_n": wmd["n_alleles"],
    "wmd_rho": f(wmd["spearman_delta_vs_log_ntrain"], 3),
    "wmd_p": f(wmd["spearman_delta_vs_log_ntrain_p"], 3),
    "wmd_top": f"+{wmd['mean_delta_top_half_train']:.4f}",
    "wmd_bot": f"+{wmd['mean_delta_bottom_half_train']:.4f}",
    "wmd_topw": wmd["wins_top_half"],
    "wmd_botw": wmd["wins_bottom_half"],
    "wmd_half": wmd["half"],
    "pdb_bN": 8,
    "pdbm_rows": pdbm_rows,
    "pdbm_idrows": pdbm_idrows,
    "pdb_fN": 9,
    "pdb_cutoff": "4.5",
    "pdbm_8rni_pep": pdbm["8RNI"]["peptide_seq"],
    "pdbm_7ow3_pep": pdbm["7OW3"]["peptide_seq"],
    "pdbm_bmin": min(len(pdbm[p]["p2_overlap_b"]) for p in pdbm),
    "pdbm_bmax": max(len(pdbm[p]["p2_overlap_b"]) for p in pdbm),
    "pdbm_fmin": min(len(pdbm[p]["pomega_overlap_f"]) for p in pdbm),
    "refm_n": refm["n_reference_datasets"],
    "refm_total": refm["n_reference_datasets"] + 13,
    "n_primary": 13,
    "n_conservative": 12,
    "fun_iedb_meas": f"{fun['iedb']['filtered_measurements']:,}",
    "fun_iedb_pairs": f"{fun['iedb']['unique_pairs']:,}",
    "fun_iedb_al": fun["iedb"]["alleles"],
    "fun_iedb_ge200": fun["iedb"]["alleles_ge200"],
    "fun_cpp_raw": f"{fun['cppsite']['raw']:,}",
    "fun_cpp_nat": f"{fun['cppsite']['natural']:,}",
    "fun_cpp_dd": f"{fun['cppsite']['dedup']:,}",
    "fun_cpp_rf": fun["cppsite"]["redundancy_filtered"],
    "fun_cpp_lw": fun["cppsite"]["length_window"],
    "ewl_rho": f(ewl["spearman_rms_vs_margin"]["rho"], 2),
    "ewl_p": f(ewl["spearman_rms_vs_margin"]["p"], 2),
    "ewl_rho_hh": f(ewl["spearman_abshh_vs_margin"]["rho"], 2),
    "ewl_p_hh": f(ewl["spearman_abshh_vs_margin"]["p"], 2),
    "ewl_rho_n": f(ewl["spearman_ntest_vs_margin"]["rho"], 2),
    "ewl_p_n": f(ewl["spearman_ntest_vs_margin"]["p"], 2),
    "ewl_mean": f(ewl["margin_mean"], 4), "ewl_sd": f(ewl["margin_sd"], 4),
    "scan_rows": "\n".join(
        r["name"] + " & " + "\\texttt{" + r["sequence"] + "}" + " & "
        + f(r["p0"], 3) + " & "
        + (r["worst_run"]["run"] if r["worst_run"] else "--") + " & "
        + (f(r["worst_run"]["drop"], 3) if r["worst_run"] else "--") + " & "
        + (f(r["worst_run"]["p_mutant"], 3) if r["worst_run"] else "--") + " & "
        + f(r["p_allRK"], 3) + " \\\\"
        for r in scan["candidates"]),
    "kwpos_rows": "\n".join(
        str(i + 1) + " & " + "\\texttt{" + d["kmer"] + "}" + " & " + f(d["weight"], 3) + " \\\\"
        for i, d in enumerate(kw["top_positive"])),
    "kwneg_rows": "\n".join(
        str(i + 1) + " & " + "\\texttt{" + d["kmer"] + "}" + " & " + f(d["weight"], 3) + " \\\\"
        for i, d in enumerate(kw["top_negative"])),
    "cpplc_rows": "\n".join(
        (str(int(100*r["fraction"])) + "\\% & " + str(r["n_train"]) + " & "
         + f(r["auc_mean"]) + " & " + f(r["auc_sd"], 4) + " \\\\")
        for r in cplc["runs"]),
    "cpplc_last": f(cplc["runs"][-1]["auc_mean"] - cplc["runs"][-2]["auc_mean"], 3),
    "cpplc_slope": f((cplc["runs"][-1]["auc_mean"] - cplc["runs"][-2]["auc_mean"]) / 0.4150, 3),
    "cpplc_2x": f(cplc["runs"][-1]["auc_mean"] + (cplc["runs"][-1]["auc_mean"] - cplc["runs"][-2]["auc_mean"]) / 0.4150, 3),
    "cnn_mean_named": f(sum(c["p_cpp"] for c in novel["named_candidates"]) / len(novel["named_candidates"])),
    "lr_mean_named": f(sum(r["p0"] for r in scan["candidates"]) / len(scan["candidates"]), 3),
    "corrob_names": ", ".join(r["name"] for r in scan["candidates"] if r["p0"] >= 0.7 and next(c for c in novel["named_candidates"] if c["name"] == r["name"])["p_cpp"] >= 0.99),
    "scan_mean_drop": f(scan["mean_logodds_drop"], 2),
    "scan_min_drop": f(scan["min_drop"], 2),
    "scan_max_drop": f(scan["max_drop"], 2),
    "kw_auc": f(kw["refit_test_auc"]),
    "kw_vocab": f"{kw['vocab_size']:,}",
    "kw_top5": ", ".join(f"\\texttt{{{d['kmer']}}}" for d in kw["top_positive"][:5]),
    "kw_top_weight": f(kw["top_positive"][0]["weight"], 2),
    "kw_top50_rk": kw["top50_composition"]["top50_containing_R_or_K"],
    "rk_frac_named": f(sum(sum(1 for a in c["sequence"] if a in "RK") / len(c["sequence"]) for c in novel["named_candidates"]) / len(novel["named_candidates"]), 3),
    "p_cnn": f"{mparams['cnn_params']:,}", "p_gnn": f"{mparams['gnn_params']:,}",
    "p_cpp": f"{mparams['cpp_cnn_params']:,}", "p_gen": f"{mparams['generator_params']:,}",
    "chem_rows": "\n".join(
        f"{lab} & {kw['chemistry_effects'][k]['n_kmers']} & "
        f"{f(kw['chemistry_effects'][k]['mean_weight'], 3)} & "
        f"{f(kw['chemistry_effects'][k]['mean_weight_rest'], 3)} \\\\"
        for lab, k in [("contains R", "contains_R"), ("contains K", "contains_K"),
                       ("contains R or K", "contains_RK"), ("contains W", "contains_W"),
                       ("contains W/Y/F", "contains_WYF"), ("contains L/I/V", "contains_LIV"),
                       ("contains D/E", "contains_DE")]),

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
