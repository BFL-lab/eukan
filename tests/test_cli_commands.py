"""CLI command bodies: what each command checks before handing its driver a config."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

import eukan.assembly
from eukan.cli import cli


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    """An empty run dir with stub inputs, as the current directory."""
    for name in ("g.fa", "l.fq", "r.fq", "s.fq"):
        (tmp_path / name).write_text(">x\nACGT\n")
    monkeypatch.chdir(tmp_path)
    return tmp_path


class TestAssembleChecks:
    @pytest.mark.parametrize(("args", "message"), [
        (["-l", "l.fq"], "Paired-end mode requires both --left and --right."),
        (["-s", "s.fq", "-S", "RF"],
         "Paired-end strand types (RF/FR) cannot be used with single-end reads."),
        (["-l", "l.fq", "-r", "r.fq", "-S", "R"],
         "Single-end strand types (R/F) cannot be used with paired-end reads."),
    ])
    def test_inconsistent_read_pairing_stops_before_the_pipeline(
        self, run_dir, monkeypatch, args, message,
    ):
        monkeypatch.setattr(eukan.assembly, "run_assembly", lambda *a, **k: pytest.fail("ran"))
        result = CliRunner().invoke(cli, ["assemble", "-g", "g.fa", *args])
        assert result.exit_code == 1
        assert "Error: invalid configuration" in result.output
        assert message in result.output

    @pytest.mark.parametrize(("args", "message"), [
        ([], "Provide --left/--right (paired) or --single reads."),
        (["-s", "s.fq", "--memory-gb", "0"], "--memory-gb must be at least 1 GiB."),
        (["-s", "s.fq", "--defuse"], "--defuse needs a protein DB; supply --uniprot."),
    ])
    def test_usage_errors(self, run_dir, args, message):
        result = CliRunner().invoke(cli, ["assemble", "-g", "g.fa", *args])
        assert result.exit_code == 2
        assert message in result.output
