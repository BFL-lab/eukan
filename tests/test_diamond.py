"""Unit tests for eukan.assembly.diamond."""

from __future__ import annotations

from eukan.assembly.diamond import coding_strand


def test_coding_strand():
    assert coding_strand("+", 2) == "+"
    assert coding_strand("+", -1) == "-"
    assert coding_strand("-", 1) == "-"   # reverse-labelled query, forward frame
    assert coding_strand("-", -2) == "+"  # antisense-assembled rescue
    assert coding_strand(".", 3) == "+"   # '.' stitched forward, like '+'
    assert coding_strand(".", -3) == "-"
