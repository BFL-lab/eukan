"""Tests for eukan.annotation.evidence — GFF3 source-token extraction.

The consensus engine matches a staged evidence file to its weights.txt entry by
the GFF3 source column (col 2) of its first data line, so the TRANSCRIPT weight
tracks whatever source token the file actually carries.
"""

from __future__ import annotations

from pathlib import Path

from eukan.annotation.evidence import _first_source_token


def _make_gff(path: Path, source: str) -> Path:
    path.write_text(
        "##gff-version 3\n"
        f"chr1\t{source}\tgene\t1\t300\t.\t+\t.\tID=g1\n"
        f"chr1\t{source}\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        f"chr1\t{source}\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=g1.t1\n"
    )
    return path


class TestFirstSourceToken:
    def test_returns_column_two(self, tmp_path):
        gff = _make_gff(tmp_path / "x.gff3", "combinr-assembly")
        assert _first_source_token(gff) == "combinr-assembly"

    def test_skips_comments_and_blank_lines(self, tmp_path):
        path = tmp_path / "x.gff3"
        path.write_text(
            "##gff-version 3\n"
            "# some comment\n"
            "\n"
            "chr1\tgenemark\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        )
        assert _first_source_token(path) == "genemark"

    def test_returns_none_for_header_only(self, tmp_path):
        path = tmp_path / "x.gff3"
        path.write_text("##gff-version 3\n# no data\n")
        assert _first_source_token(path) is None
