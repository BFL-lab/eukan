"""Tests for eukan.stats.format — terminal formatting and TSV export."""

from __future__ import annotations

from pathlib import Path

import pytest

from eukan.stats import (
    compare_annotations,
    compare_multiple,
    format_multi_results,
    format_results,
    write_details_tsv,
    write_stats_tsv,
)
from eukan.stats.models import PAIR_TEST_TSV_COLUMNS, TSV_COLUMNS

# Reused fixtures from test_compare; duplicated here to keep test files
# independent.
_REF_GFF = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
"""

_P_MISSING = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
"""


def _write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


class TestFormatResultsSinglePred:
    def test_contains_expected_sections(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        result = compare_annotations(ref, ref)
        out = format_results(result)
        assert "ANNOTATION QUALITY ASSESSMENT" in out
        assert "GENE Level" in out
        assert "mRNA Level" in out
        assert "CDS Level" in out
        assert "Intron Level" in out
        assert "SUMMARY" in out

    def test_paths_appear_in_output(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        pred = _write(tmp_path / "pred.gff3", _REF_GFF)
        result = compare_annotations(ref, pred)
        out = format_results(result)
        assert str(ref) in out
        assert str(pred) in out


class TestFormatMultiResults:
    @pytest.fixture
    def two_preds_result(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _REF_GFF)
        p2 = _write(tmp_path / "p2.gff3", _P_MISSING)
        return compare_multiple(ref, [("p1", p1), ("p2", p2)])

    def test_contains_top_header(self, two_preds_result):
        out = format_multi_results(two_preds_result)
        assert "ANNOTATION QUALITY ASSESSMENT (multi-prediction)" in out

    def test_lists_each_prediction_in_header(self, two_preds_result):
        out = format_multi_results(two_preds_result)
        # Predictions are listed at the top
        first_section = out.split("### PREDICTION:")[0]
        assert "p1" in first_section
        assert "p2" in first_section

    def test_per_prediction_banners(self, two_preds_result):
        out = format_multi_results(two_preds_result)
        assert "### PREDICTION: p1 ###" in out
        assert "### PREDICTION: p2 ###" in out

    def test_comparative_section_present(self, two_preds_result):
        out = format_multi_results(two_preds_result)
        assert "COMPARATIVE SUMMARY" in out
        assert "F1 by level" in out
        assert "Cohen's kappa" in out
        assert "Powerset of gene-level matches" in out
        assert "Significance tests" in out

    def test_f1_table_has_each_pred_label(self, two_preds_result):
        out = format_multi_results(two_preds_result)
        f1_section = out.split("F1 by level")[1].split("Cohen's kappa")[0]
        assert "p1" in f1_section
        assert "p2" in f1_section


class TestWriteDetailsTsvSingle:
    def test_columns_match_schema(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        result = compare_annotations(ref, ref)
        out_path = tmp_path / "details.tsv"
        write_details_tsv(result, out_path)
        first_line = out_path.read_text().splitlines()[0]
        assert first_line.split("\t") == list(TSV_COLUMNS)

    def test_one_row_per_record(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        result = compare_annotations(ref, ref)
        out_path = tmp_path / "details.tsv"
        write_details_tsv(result, out_path)
        n_lines = len(out_path.read_text().splitlines())
        assert n_lines == 1 + len(result.records)


class TestWriteDetailsTsvMulti:
    def test_prediction_column_prepended(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _REF_GFF)
        p2 = _write(tmp_path / "p2.gff3", _P_MISSING)
        result = compare_multiple(ref, [("p1", p1), ("p2", p2)])
        out_path = tmp_path / "details.tsv"
        write_details_tsv(result, out_path)
        first_line = out_path.read_text().splitlines()[0]
        assert first_line.split("\t") == ["prediction", *TSV_COLUMNS]

    def test_concatenates_records_per_pred(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _REF_GFF)
        p2 = _write(tmp_path / "p2.gff3", _P_MISSING)
        result = compare_multiple(ref, [("p1", p1), ("p2", p2)])
        out_path = tmp_path / "details.tsv"
        write_details_tsv(result, out_path)
        lines = out_path.read_text().splitlines()
        rows = [line.split("\t") for line in lines[1:]]
        # Each pred contributed records, all tagged with the right label
        labels = {row[0] for row in rows}
        assert labels == {"p1", "p2"}
        expected = sum(len(p.records) for p in result.per_prediction)
        assert len(rows) == expected


class TestWriteStatsTsv:
    def test_columns_match_schema(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _REF_GFF)
        p2 = _write(tmp_path / "p2.gff3", _P_MISSING)
        result = compare_multiple(ref, [("p1", p1), ("p2", p2)])
        out_path = tmp_path / "stats.tsv"
        write_stats_tsv(result, out_path)
        first_line = out_path.read_text().splitlines()[0]
        assert first_line.split("\t") == list(PAIR_TEST_TSV_COLUMNS)

    def test_one_row_per_pair_test(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _REF_GFF)
        p2 = _write(tmp_path / "p2.gff3", _P_MISSING)
        result = compare_multiple(ref, [("p1", p1), ("p2", p2)])
        out_path = tmp_path / "stats.tsv"
        write_stats_tsv(result, out_path)
        n_data_rows = len(out_path.read_text().splitlines()) - 1
        assert n_data_rows == len(result.pair_tests)
