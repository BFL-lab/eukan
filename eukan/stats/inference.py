"""Statistical primitives used by inter-prediction comparison.

Thin scipy adapters that normalize edge cases (empty inputs, zero marginals)
and return plain Python types, so callers don't deal with scipy result
classes or empty-input quirks. Cohen's kappa is the only piece scipy
doesn't ship; it's a few lines on a confusion matrix.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import nan

import numpy as np
from scipy.stats import chi2_contingency as _chi2_contingency
from scipy.stats import false_discovery_control as _false_discovery_control
from scipy.stats import ks_2samp as _ks_2samp


def ks_2samp(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    """Two-sample Kolmogorov-Smirnov test.

    Returns ``(D, p_value)``. If either side has length 0, returns
    ``(0.0, 1.0)`` — nothing to distinguish.
    """
    if len(x) == 0 or len(y) == 0:
        return 0.0, 1.0
    res = _ks_2samp(x, y, method="auto")
    return float(res.statistic), float(res.pvalue)


def chi2_contingency(table: Sequence[Sequence[float]]) -> tuple[float, float, int]:
    """Pearson chi-squared test on a 2-D contingency table.

    Returns ``(chi2, p_value, dof)``. If any row or column has a zero
    marginal (so a category was absent in both groups), returns
    ``(0.0, 1.0, 0)`` rather than letting scipy raise on the degenerate
    table.
    """
    arr = np.asarray(table, dtype=float)
    if arr.ndim != 2 or arr.size == 0:
        return 0.0, 1.0, 0
    if np.any(arr.sum(axis=1) == 0) or np.any(arr.sum(axis=0) == 0):
        return 0.0, 1.0, 0
    res = _chi2_contingency(arr)
    return float(res.statistic), float(res.pvalue), int(res.dof)


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values, as a 1-D NumPy array."""
    pv = np.asarray(p_values, dtype=float)
    if pv.size == 0:
        return pv
    return np.asarray(_false_discovery_control(pv, method="bh"), dtype=float)


def cohen_kappa(labels_a: Sequence[str], labels_b: Sequence[str]) -> float:
    """Cohen's kappa from two equal-length label sequences.

    Returns ``nan`` for empty input (kappa is undefined). Returns ``1.0``
    when both labelers used a single shared category for every observation
    — the formula's 0/0 case, conventionally treated as perfect agreement.
    """
    if len(labels_a) != len(labels_b):
        raise ValueError("labels_a and labels_b must have the same length")
    n = len(labels_a)
    if n == 0:
        return nan

    categories = sorted(set(labels_a) | set(labels_b))
    cat_to_idx = {c: i for i, c in enumerate(categories)}
    k = len(categories)
    cm = np.zeros((k, k), dtype=float)
    for a, b in zip(labels_a, labels_b, strict=True):
        cm[cat_to_idx[a], cat_to_idx[b]] += 1.0

    total = cm.sum()
    p_o = float(np.trace(cm) / total)
    row_marg = cm.sum(axis=1) / total
    col_marg = cm.sum(axis=0) / total
    p_e = float((row_marg * col_marg).sum())
    if p_e >= 1.0:
        return 1.0
    return (p_o - p_e) / (1.0 - p_e)
