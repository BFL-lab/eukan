"""Unit tests for eukan.gff.rows (line-level GFF3/GTF reading)."""

from __future__ import annotations

from eukan.gff.rows import iter_gff3_rows, parse_attrs


def test_parse_attrs():
    assert parse_attrs("ID=m1;Parent=g1;contains=a,b") == {
        "ID": "m1", "Parent": "g1", "contains": "a,b"
    }


def test_parse_attrs_gtf():
    assert parse_attrs('gene_id "g1"; transcript_id "t1"; cov "2.5";') == {
        "gene_id": "g1", "transcript_id": "t1", "cov": "2.5"
    }


def test_parse_attrs_strips_and_skips_empty_parts():
    assert parse_attrs(" ID = a ;;Note=x=y; \n") == {"ID": "a", "Note": "x=y"}


def test_iter_gff3_rows_skips_comments_blanks_and_short_rows(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "# comment\n"
        "\n"
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g\n"
        "chr1\tsrc\tgene\t1\t9\n"
        "   \n"
        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tParent=g\textra"
    )
    assert list(iter_gff3_rows(gff)) == [
        ["chr1", "src", "gene", "1", "9", ".", "+", ".", "ID=g"],
        ["chr1", "src", "exon", "1", "9", ".", "+", ".", "Parent=g", "extra"],
    ]
