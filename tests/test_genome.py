"""Regression tests for ContigIndex (eukan.infra.genome)."""

from __future__ import annotations

from pathlib import Path

import pytest

from eukan.infra.genome import ContigIndex


def _fa(tmp_path: Path) -> Path:
    p = tmp_path / "g.fa"
    p.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTTGGGGCC\n")
    return p


def test_hit_and_single_record_cache(tmp_path):
    with ContigIndex(_fa(tmp_path)) as c:
        assert str(c["chr1"].seq) == "ACGTACGTAC"
        assert str(c["chr2"].seq) == "TTTTGGGGCC"
        assert str(c["chr1"].seq) == "ACGTACGTAC"  # re-hit after cache evicted


def test_get_missing_returns_default(tmp_path):
    with ContigIndex(_fa(tmp_path)) as c:
        assert "chr1" in c
        assert "absent" not in c
        assert c.get("absent") is None


def test_repeated_missing_lookup_does_not_corrupt_cache(tmp_path):
    # Regression: a missing contig must raise KeyError without leaving
    # _cache_id set to the missing id while _cache_record stays None — that
    # turned a second same-id lookup into an AssertionError that get() does
    # not catch (the bug that crashed the segemehl-on-wrong-genome walk).
    with ContigIndex(_fa(tmp_path)) as c:
        assert c.get("absent") is None
        assert c.get("absent") is None              # second miss must not raise
        assert str(c["chr1"].seq) == "ACGTACGTAC"   # real lookups still work
        with pytest.raises(KeyError):
            _ = c["absent"]                          # KeyError, not AssertionError
