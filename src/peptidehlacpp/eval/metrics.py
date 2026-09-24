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
    """Partial AUC over FPR in [0, frac] (the literature's AUC0.1), normalized
    to [0, 1]: TPR integrated against FPR up to `frac`, divided by `frac`."""
    y = labels.astype(bool)
    npos, nneg = y.sum(), (~y).sum()
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y = y[order]
    tpr = np.concatenate([[0.0], np.cumsum(y) / npos])
    fpr = np.concatenate([[0.0], np.cumsum(~y) / nneg])
    # integrate TPR dFPR from 0 to frac with linear interpolation at the cut
    mask = fpr <= frac
    area = np.trapezoid(tpr[mask], fpr[mask])
    if fpr[-1] > frac and mask.any():
        i = mask.sum()
        if i < len(fpr) and fpr[i] > fpr[i - 1]:
            t_at = tpr[i - 1] + (tpr[i] - tpr[i - 1]) * (frac - fpr[i - 1]) / (fpr[i] - fpr[i - 1])
            area += 0.5 * (tpr[i - 1] + t_at) * (frac - fpr[i - 1])
    return float(area / frac)


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
