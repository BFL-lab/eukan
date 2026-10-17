"""CLI command bodies: what each command checks before handing its driver a config."""

from __future__ import annotations

import os
import subprocess
import sys

import pytest
from click.testing import CliRunner

import eukan.assembly
from eukan.cli import cli
from eukan.infra.utils import available_cpus


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


class TestOptionsAreConfigFields:
    """Each command hands its options to its config by name, so every option
    except the re-run flags and output modifiers must be named after a field."""

    @pytest.mark.parametrize(("command", "config"), [
        ("annotate", "PipelineConfig"),
        ("assemble", "AssemblyConfig"),
        ("func-annot", "FunctionalConfig"),
        ("mask-repeats", "RepeatsConfig"),
        ("prep-submission", "SubmissionConfig"),
    ])
    def test_destinations_are_field_names(self, command, config):
        import eukan.settings

        fields = getattr(eukan.settings, config).model_fields
        # click-option-group adds a hidden parameter per group heading
        params = [p for p in cli.commands[command].params if not getattr(p, "hidden", False)]
        names = {p.name for p in params} - {"force", "print_command", "dry_run"}
        assert {n for n in names if not n.startswith("run_")} - set(fields) == set()

    def test_renamed_options_reach_the_config(self, run_dir, monkeypatch):
        seen = []
        monkeypatch.setattr(eukan.assembly, "run_assembly", lambda c, **k: seen.append(c))
        result = CliRunner().invoke(cli, [
            "assemble", "-g", "g.fa", "-l", "l.fq", "-r", "r.fq", "-n", "3", "-c", "6",
            "-m", "30", "-M", "900", "--phred", "64", "--adapter-sequence", "ACGT",
        ])
        assert result.exit_code == 0, result.output
        (config,) = seen
        assert (config.left_reads, config.right_reads) == (run_dir / "l.fq", run_dir / "r.fq")
        assert (config.num_cpu, config.genetic_code) == (3, "6")
        assert (config.min_intron_len, config.max_intron_len) == (30, 900)
        assert (config.phred_quality, config.adapter_sequences) == (64, ["ACGT"])


class TestCodeTableEpilog:
    @pytest.mark.parametrize(("flag", "shown"), [("--help", True), ("-h", False)])
    def test_only_long_help_lists_the_genetic_codes(self, flag, shown):
        result = CliRunner().invoke(cli, ["gff3toseq", flag])
        assert result.exit_code == 0
        assert ("Genetic codes (NCBI translation tables):" in result.output) is shown


def test_gff3toseq_rejects_an_unknown_genetic_code(run_dir):
    """-c 99 is a configuration error, not a KeyError out of the translation."""
    (run_dir / "chr1.fa").write_text(">chr1\nATGAAATAG\n")
    (run_dir / "genes.gff3").write_text(
        "chr1\tt\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\tt\tmRNA\t1\t9\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\tt\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=m1\n"
    )
    result = CliRunner().invoke(
        cli, ["gff3toseq", "-g", "chr1.fa", "-i", "genes.gff3", "-c", "99"],
    )
    assert result.exit_code == 1
    assert "Error: Unknown NCBI genetic code: 99" in result.output
    assert "Hint: Valid codes: 1, 2, 3," in result.output


def test_db_fetch_help_has_no_options_header():
    """db-fetch lays out its option groups like every other command."""
    result = CliRunner().invoke(cli, ["db-fetch", "-h"])
    assert "Options:" not in result.output
    assert "\nPipeline parameters: \n-o, --output-dir PATH" in result.output


@pytest.mark.parametrize(("args", "target"), [
    (["db-fetch", "--homology-db", "KOFAM", "-d", "PFam"],
     ("eukan.functional.dbfetch", "fetch_databases")),
    (["check", "--homology-db", "KOFAM", "--db-dir", "."], ("eukan.infra.health", "run_checks")),
])
def test_choices_arrive_as_declared(run_dir, monkeypatch, args, target):
    """A case-insensitive Choice hands over the choice as declared (lower case)."""
    import importlib

    seen = []

    def stub(*a, **k):
        seen.append((k["homology_db"], k.get("databases")))
        return ([], [], [], [])

    monkeypatch.setattr(importlib.import_module(target[0]), target[1], stub)
    CliRunner().invoke(cli, args)
    assert seen == [("kofam", ["pfam"] if args[0] == "db-fetch" else None)]


@pytest.mark.skipif(not hasattr(os, "sched_setaffinity"), reason="needs CPU affinity")
def test_numcpu_default_follows_cpu_affinity():
    """A job pinned to some CPUs (SLURM cpuset, docker --cpuset-cpus) defaults to
    that many threads, not to every CPU in the machine."""
    code = (
        "import os\n"
        "os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})\n"
        "from eukan.cli.annotate import annotate\n"
        "print(next(p.default for p in annotate.params if p.name == 'num_cpu'))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True,
    ).stdout
    assert out.strip() == "1"


def test_available_cpus_falls_back_to_cpu_count(monkeypatch):
    """Platforms without sched_getaffinity (macOS) use the machine's count."""
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    assert available_cpus() == (os.cpu_count() or 1)
