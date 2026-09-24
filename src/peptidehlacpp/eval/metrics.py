"""Evaluation metrics standard in the peptide-HLA binding literature."""
from __future__ import annotations

import numpy as np
from scipy import stats


def auc(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order].astype(bool)
    npos, nneg = y.sum(), (~y).sum()
    if npos == 0 or nneg == 0:
        return float("nan")
    tps = np.cumsum(y)
    fps = np.cumsum(~y)
    return float(np.sum((1 - fps / nneg) * np.diff(np.concatenate([[0], tps / npos]))))


def auc_top_frac(labels: np.ndarray, scores: np.ndarray, frac: float = 0.1) -> float:
    """AUC restricted to the top-scoring `frac` of peptides (AUC0.1)."""
    n = max(1, int(round(frac * len(labels))))
    order = np.argsort(-scores, kind="mergesort")[:n]
    return auc(labels[order], scores[order])


def ppv(labels: np.ndarray, scores: np.ndarray, top_n: int | None = None) -> float:
    """Positive predictive value among the top-N predictions (N = #true binders)."""
    y = labels.astype(bool)
    n = top_n or max(1, int(y.sum()))
    order = np.argsort(-scores, kind="mergesort")[:n]
    return float(y[order].mean())


def srcc(values_true: np.ndarray, values_pred: np.ndarray) -> float:
    return float(stats.spearmanr(values_true, values_pred).statistic)


def rmse(values_true: np.ndarray, values_pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((values_true - values_pred) ** 2)))
