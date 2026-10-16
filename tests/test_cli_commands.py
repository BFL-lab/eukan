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


class TestPathsReachTheCommandAbsolute:
    """Relative and symlinked paths reach the configs and helpers resolved."""

    @pytest.fixture
    def tree(self, run_dir):
        real = run_dir / "real"
        real.mkdir()
        for name in ("g.fa", "p.faa", "x.sbt", "lib.fa", "ref.gff3"):
            (real / name).write_text(">x\nACGT\n")
            (run_dir / f"ln-{name}").symlink_to(real / name)
        return real.resolve()

    @pytest.mark.parametrize(("args", "target", "pick", "expected"), [
        (["annotate", "-g", "ln-g.fa", "-p", "real/../ln-p.faa"],
         ("eukan.annotation", "run_annotation_pipeline"),
         lambda a, k: (a[0].genome, a[0].proteins[0]), ("g.fa", "p.faa")),
        (["assemble", "-g", "ln-g.fa", "-s", "ln-p.faa", "--uniprot", "real/lib.fa"],
         ("eukan.assembly", "run_assembly"),
         lambda a, k: (a[0].genome, a[0].single_reads, a[0].uniprot_db),
         ("g.fa", "p.faa", "lib.fa")),
        (["func-annot", "-p", "ln-p.faa", "--pfam", "ln-lib.fa"],
         ("eukan.functional", "run_functional_annotation"),
         lambda a, k: (a[0].proteins, a[0].pfam_db), ("p.faa", "lib.fa")),
        (["mask-repeats", "-g", "ln-g.fa", "--lib", "ln-lib.fa"],
         ("eukan.repeats", "run_repeats"),
         lambda a, k: (a[0].genome, a[0].lib), ("g.fa", "lib.fa")),
        (["prep-submission", "-t", "ln-x.sbt", "--organism", "Foo bar", "-g", "ln-g.fa",
          "-i", "ln-ref.gff3", "-d", "real/out"],
         ("eukan.submission", "run_prep_submission"),
         lambda a, k: (a[0].template, a[0].genome, a[0].gff3, a[0].output_dir),
         ("x.sbt", "g.fa", "ref.gff3", "out")),
        (["db-fetch", "-o", "real/db"], ("eukan.functional.dbfetch", "fetch_databases"),
         lambda a, k: (a[0],), ("db",)),
        (["check", "--db-dir", "real/db"], ("eukan.infra.health", "run_checks"),
         lambda a, k: (k["db_dir"],), ("db",)),
        (["status", "-d", "real"], ("eukan.infra.manifest", "load_manifest"),
         lambda a, k: (a[0] / "x",), ("x",)),
    ])
    def test_resolved(self, tree, monkeypatch, args, target, pick, expected):
        import importlib

        seen = []

        def stub(*a, **k):
            seen.append(pick(a, k))
            return ([], [], [], []) if args[0] == "check" else None

        monkeypatch.setattr(importlib.import_module(target[0]), target[1], stub)
        CliRunner().invoke(cli, args)
        assert seen == [tuple(tree / name for name in expected)]
