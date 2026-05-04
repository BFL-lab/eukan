"""Tests for eukan.stats.inference — statistical primitives."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import chi2_contingency as scipy_chi2
from scipy.stats import false_discovery_control as scipy_bh

from eukan.stats.inference import bh_fdr, chi2_contingency, cohen_kappa, ks_2samp


class TestKs2samp:
    def test_identical_arrays(self):
        d, p = ks_2samp([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0])
        assert d == 0.0
        assert p == pytest.approx(1.0, rel=1e-6)

    def test_empty_either_side(self):
        assert ks_2samp([], [1.0, 2.0]) == (0.0, 1.0)
        assert ks_2samp([1.0, 2.0], []) == (0.0, 1.0)
        assert ks_2samp([], []) == (0.0, 1.0)

    def test_clearly_different(self):
        rng = np.random.default_rng(0)
        x = rng.uniform(0.0, 0.4, 100).tolist()
        y = rng.uniform(0.6, 1.0, 100).tolist()
        d, p = ks_2samp(x, y)
        assert d == pytest.approx(1.0)  # complete separation
        assert p < 1e-20

    def test_returns_python_floats(self):
        d, p = ks_2samp([1.0, 2.0, 3.0], [1.5, 2.5, 3.5])
        assert isinstance(d, float)
        assert isinstance(p, float)


class TestChi2Contingency:
    def test_known_2x2(self):
        chi2, p, dof = chi2_contingency([[10, 20], [30, 40]])
        e_chi2, e_p, e_dof, _ = scipy_chi2(np.array([[10, 20], [30, 40]]))
        assert chi2 == pytest.approx(float(e_chi2))
        assert p == pytest.approx(float(e_p))
        assert dof == int(e_dof)

    def test_2x6_classification_table(self):
        # Realistic gene-classification contingency: 2 preds by 6 categories
        # (exact, inexact, missing, merged, fragmented, novel).
        table = [[120, 30, 40, 5, 5, 100], [100, 50, 50, 10, 10, 80]]
        chi2, p, dof = chi2_contingency(table)
        e = scipy_chi2(np.array(table))
        assert chi2 == pytest.approx(float(e.statistic))
        assert p == pytest.approx(float(e.pvalue))
        assert dof == int(e.dof)

    def test_zero_row_marginal(self):
        chi2, p, dof = chi2_contingency([[0, 0, 0], [10, 20, 30]])
        assert (chi2, p, dof) == (0.0, 1.0, 0)

    def test_zero_col_marginal(self):
        chi2, p, dof = chi2_contingency([[10, 0, 30], [20, 0, 40]])
        assert (chi2, p, dof) == (0.0, 1.0, 0)

    def test_empty_table(self):
        chi2, p, dof = chi2_contingency([])
        assert (chi2, p, dof) == (0.0, 1.0, 0)


class TestBhFdr:
    def test_matches_scipy(self):
        p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205]
        adj = bh_fdr(p)
        expected = scipy_bh(np.asarray(p), method="bh")
        np.testing.assert_allclose(adj, expected)

    def test_empty(self):
        adj = bh_fdr([])
        assert adj.size == 0

    def test_monotone_after_adjustment(self):
        # Adjusted p-values must be monotone non-decreasing in input rank.
        p = [0.5, 0.01, 0.05, 0.005, 0.2]
        adj = bh_fdr(p)
        order = np.argsort(p)
        diffs = np.diff(adj[order])
        assert np.all(diffs >= -1e-12)

    def test_returns_ndarray_of_floats(self):
        adj = bh_fdr([0.01, 0.5])
        assert isinstance(adj, np.ndarray)
        assert adj.dtype == np.float64


class TestCohenKappa:
    def test_perfect_agreement(self):
        a = ["match", "missing", "exact", "match"]
        assert cohen_kappa(a, a) == 1.0

    def test_complete_disagreement_two_categories(self):
        # Uniform marginals, zero diagonal -> kappa = -1
        a = ["x"] * 5 + ["y"] * 5
        b = ["y"] * 5 + ["x"] * 5
        assert cohen_kappa(a, b) == pytest.approx(-1.0)

    def test_single_category_agreement_returns_one(self):
        a = ["match"] * 10
        b = ["match"] * 10
        assert cohen_kappa(a, b) == 1.0

    def test_empty_returns_nan(self):
        assert math.isnan(cohen_kappa([], []))

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            cohen_kappa(["a"], ["a", "b"])

    def test_partial_agreement_in_unit_interval(self):
        # 4/5 agree on x, 4/5 agree on y -> kappa positive but < 1
        a = ["x", "x", "x", "x", "x", "y", "y", "y", "y", "y"]
        b = ["x", "x", "x", "x", "y", "y", "y", "y", "y", "x"]
        k = cohen_kappa(a, b)
        assert 0.0 < k < 1.0
