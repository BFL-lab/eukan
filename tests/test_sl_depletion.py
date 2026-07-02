"""Unit tests for eukan.assembly.sl_depletion SL-motif primitives."""

from __future__ import annotations

from eukan.assembly.sl_depletion import (
    _find_sites,
    _revcomp,
    _variants,
)

SL = "GTACTTTATT"  # 10 bp, meets _MIN_MOTIF_LEN


def _patterns(motif, k=0):
    return _variants(motif, k) | _variants(_revcomp(motif), k)


# --- low-level helpers -----------------------------------------------------


def test_revcomp():
    assert _revcomp("GTACTTTATT") == "AATAAAGTAC"


def test_variants_counts():
    assert _variants("ACGT", 0) == {"ACGT"}
    # 1 substitution: original + 3 alternatives per position
    assert len(_variants("ACGT", 1)) == 1 + 3 * 4


def test_find_sites_exact_and_merge():
    seq = "A" * 20 + SL + "C" * 20
    sites = _find_sites(seq, _patterns(SL, 0), len(SL))
    assert sites == [(20, 30)]


def test_find_sites_reverse_complement():
    rc = _revcomp(SL)
    seq = "G" * 15 + rc + "T" * 15
    sites = _find_sites(seq, _patterns(SL), len(SL))
    assert sites == [(15, 25)]


def test_find_sites_tolerates_one_mismatch():
    mutated = "G" + "A" + SL[2:]  # one substitution at position 1
    seq = "A" * 20 + mutated + "C" * 20
    assert _find_sites(seq, _patterns(SL, 1), len(SL)) == [(20, 30)]
    # exact-only search misses it
    assert _find_sites(seq, _patterns(SL, 0), len(SL)) == []
