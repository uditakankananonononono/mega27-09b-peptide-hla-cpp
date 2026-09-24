# Appendix A: Identifiability of additive binding models

**Setup.** Peptides of length L over alphabet A (|A|=20). The additive
(position-specific scoring) model assigns free energy
  G_add(x; w) = b + Σ_{i=1..L} w(i, x_i),  w ∈ R^{L×20}.
Data: pairs (x^(n), g^(n)) with g^(n) = log10 IC50 (linear in ΔG by §2.1).

**Proposition 1 (identifiability up to gauge).**
Let X ∈ {0,1}^{n×20L} be the one-hot design matrix of the data. Then:
(i) any two additive parameter vectors w, w' giving identical predictions on
    all peptides differ by a *gauge*: w'(i,a) = w(i,a) + c_i for constants
    c_1..c_L (absorbed into b);
(ii) if X has full column rank (guaranteed when n ≥ 20L and every amino acid
    occurs at every position in the training set), then modulo the gauge,
    the least-squares solution is unique and equals the projection of the
    true energy function onto the additive subspace.

*Proof.* (i) Suppose w - w' = Δ predicts zero on every sequence. Take any
position i and residues a, a'. Consider a background sequence x and its
single-mutant x' at position i: predictions differ by Δ(i,a) - Δ(i,a') = 0,
so Δ(i,·) is constant in the residue: Δ(i,a) = c_i. Hence w' = w + gauge.
(ii) X^T X is block-diagonal in positions after the gauge quotient; each
block is the covariance of a categorical variable taking all 20 values with
positive probability, hence positive definite. Uniqueness of the ridge/OLS
solution follows; the solution is the orthogonal projection onto
span(additive) by normal equations. ∎

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
the measured PSSM ceiling (test AUC 0.9156 vs CNN 0.9057/GNN (see §5)) is a
*prediction target* for non-additive models.

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
