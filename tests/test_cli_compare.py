"""End-to-end CLI tests for the ``eukan compare`` command."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from eukan.cli import cli

_REF_GFF = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
chr1\t.\tgene\t3000\t3500\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t3000\t3500\t.\t+\t.\tID=m3;Parent=g3
chr1\t.\tCDS\t3000\t3500\t.\t+\t0\tID=c3;Parent=m3
"""

_P_PARTIAL = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
"""


@pytest.fixture
def gffs(tmp_path) -> dict[str, Path]:
    paths = {
        "ref": tmp_path / "ref.gff3",
        "p1": tmp_path / "p1.gff3",
        "p2": tmp_path / "p2.gff3",
    }
    paths["ref"].write_text(_REF_GFF)
    paths["p1"].write_text(_REF_GFF)
    paths["p2"].write_text(_P_PARTIAL)
    return paths


class TestSinglePrediction:
    def test_single_pred_runs_and_uses_legacy_layout(self, gffs):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]), "-p", str(gffs["p1"]),
        ])
        assert result.exit_code == 0, result.output
        # Legacy single-pred header (no "(multi-prediction)" suffix).
        assert "ANNOTATION QUALITY ASSESSMENT" in result.output
        assert "(multi-prediction)" not in result.output
        # No comparative section in single-pred output.
        assert "COMPARATIVE SUMMARY" not in result.output

    def test_single_pred_with_output_file_writes_tsv(self, gffs, tmp_path):
        out_path = tmp_path / "details.tsv"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]), "-p", str(gffs["p1"]),
            "-o", str(out_path),
        ])
        assert result.exit_code == 0, result.output
        assert out_path.exists()
        # Single-pred TSV does NOT have a leading "prediction" column.
        first_line = out_path.read_text().splitlines()[0]
        assert not first_line.startswith("prediction\t")

    def test_label_with_n1_emits_warning(self, gffs):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]), "-p", str(gffs["p1"]),
            "-L", "ignored",
        ])
        assert result.exit_code == 0, result.output
        # Click 8.3 captures stderr separately on result.
        assert "ignored" in result.stderr.lower()

    def test_stats_file_with_n1_errors(self, gffs, tmp_path):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]), "-p", str(gffs["p1"]),
            "--stats-file", str(tmp_path / "stats.tsv"),
        ])
        assert result.exit_code != 0
        assert "two or more" in result.output


class TestMultiPrediction:
    def test_two_preds_produces_comparative_section(self, gffs):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]),
            "-p", str(gffs["p1"]), "-p", str(gffs["p2"]),
        ])
        assert result.exit_code == 0, result.output
        assert "(multi-prediction)" in result.output
        assert "COMPARATIVE SUMMARY" in result.output
        assert "Cohen's kappa" in result.output
        assert "Powerset of gene-level matches" in result.output

    def test_labels_appear_in_output(self, gffs):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]),
            "-p", str(gffs["p1"]), "-p", str(gffs["p2"]),
            "-L", "alpha", "-L", "beta",
        ])
        assert result.exit_code == 0, result.output
        assert "### PREDICTION: alpha ###" in result.output
        assert "### PREDICTION: beta ###" in result.output

    def test_label_count_mismatch_errors(self, gffs):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]),
            "-p", str(gffs["p1"]), "-p", str(gffs["p2"]),
            "-L", "only-one",
        ])
        assert result.exit_code != 0
        assert "counts must match" in result.output

    def test_writes_multi_details_and_stats_tsv(self, gffs, tmp_path):
        out_path = tmp_path / "details.tsv"
        stats_path = tmp_path / "stats.tsv"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]),
            "-p", str(gffs["p1"]), "-p", str(gffs["p2"]),
            "-o", str(out_path), "--stats-file", str(stats_path),
        ])
        assert result.exit_code == 0, result.output
        assert out_path.exists() and stats_path.exists()
        # Multi-pred TSV has the "prediction" column prepended.
        first_line = out_path.read_text().splitlines()[0]
        assert first_line.startswith("prediction\t")
        # Stats TSV has the long-form schema.
        stats_header = stats_path.read_text().splitlines()[0]
        for col in ("pred_a", "pred_b", "level", "test", "p_adj"):
            assert col in stats_header

    def test_ecdf_metrics_all_three(self, gffs, tmp_path):
        stats_path = tmp_path / "stats.tsv"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(gffs["ref"]),
            "-p", str(gffs["p1"]), "-p", str(gffs["p2"]),
            "--ecdf-metrics", "sn",
            "--ecdf-metrics", "sp",
            "--ecdf-metrics", "f1",
            "--stats-file", str(stats_path),
        ])
        assert result.exit_code == 0, result.output
        # 4 levels * 1 pair * (3 KS + 1 chi2) = 16 stats rows
        n_data = len(stats_path.read_text().splitlines()) - 1
        assert n_data == 16

    def test_label_collision_warning_on_same_stem(self, tmp_path):
        # Two pred files with the same stem in different directories.
        a = tmp_path / "a"
        b = tmp_path / "b"
        a.mkdir()
        b.mkdir()
        ref = tmp_path / "ref.gff3"
        ref.write_text(_REF_GFF)
        (a / "pred.gff3").write_text(_REF_GFF)
        (b / "pred.gff3").write_text(_P_PARTIAL)

        runner = CliRunner()
        result = runner.invoke(cli, [
            "compare", "-r", str(ref),
            "-p", str(a / "pred.gff3"),
            "-p", str(b / "pred.gff3"),
        ])
        assert result.exit_code == 0, result.output
        assert "auto-numbered" in result.stderr
        # Auto-numbered labels (pred_1, pred_2) appear in the output.
        assert "### PREDICTION: pred_1 ###" in result.output
        assert "### PREDICTION: pred_2 ###" in result.output
