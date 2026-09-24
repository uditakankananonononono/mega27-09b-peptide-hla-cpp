"""Additive position-specific scoring model (PSSM) - the linear baseline.

Implements the classical stabilized matrix method (SMM-style ridge fit):
    log10(IC50) ~ sum_i w[pos_i][aa_i]
solved in closed form per allele. Serves as the additive-model baseline and
as the empirical anchor for the identifiability analysis in the paper.
"""
from __future__ import annotations

import numpy as np

from ..features import AA_INDEX, MAX_LEN_PHLA


class AllelePSSM:
    def __init__(self, ridge: float = 1.0, max_len: int = MAX_LEN_PHLA):
        self.ridge = ridge
        self.max_len = max_len
        self.W: np.ndarray | None = None  # (max_len, 20)
        self.b: float = 0.0

    def _design(self, seqs: list[str]) -> np.ndarray:
        n = len(seqs)
        X = np.zeros((n, self.max_len * 20), dtype=np.float32)
        for r, s in enumerate(seqs):
            for i, a in enumerate(s[: self.max_len]):
                X[r, i * 20 + AA_INDEX[a]] = 1.0
        return X

    def fit(self, seqs: list[str], log_ic50: np.ndarray) -> "AllelePSSM":
        X = self._design(seqs)
        y = np.asarray(log_ic50, dtype=np.float64)
        Xb = np.hstack([X, np.ones((X.shape[0], 1))]).astype(np.float64)
        A = Xb.T @ Xb
        A[np.diag_indices_from(A)] += self.ridge
        beta = np.linalg.solve(A, Xb.T @ y)
        self.W = beta[:-1].reshape(self.max_len, 20)
        self.b = float(beta[-1])
        return self

    def predict(self, seqs: list[str]) -> np.ndarray:
        if self.W is None:
            raise RuntimeError("model not fit")
        X = self._design(seqs)
        return X.astype(np.float64) @ self.W.ravel() + self.b
