"""Integration tests — verify install paths, imports, CLI wiring, and pyhmmer.

These tests don't require external bioinformatics tools. They verify that:
- All modules import cleanly (no circular deps, no missing deps)
- The CLI entry point works and all subcommands are registered
- pyhmmer is importable and functional
- Settings load correctly from various sources
- The tool registry loads from tools.toml
"""

import importlib
import subprocess
import sys
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Module imports — catches circular deps and missing packages
# ---------------------------------------------------------------------------


MODULES = [
    "eukan",
    "eukan.cli",
    "eukan.settings",
    "eukan.check",
    "eukan.annotation.orf",
    "eukan.infra",
    "eukan.infra.runner",
    "eukan.infra.manifest",
    "eukan.infra.steps",
    "eukan.infra.logging",
    "eukan.annotation",
    "eukan.annotation.orchestrator",
    "eukan.assembly",
    "eukan.assembly.orchestrator",
    "eukan.functional",
    "eukan.functional.orchestrator",
    "eukan.functional.dbfetch",
    "eukan.gff",
    "eukan.gff.parser",
    "eukan.gff.intersecter",
    "eukan.gff.io",
]


@pytest.mark.parametrize("module_name", MODULES)
def test_module_imports(module_name):
    """Every module should import without error."""
    mod = importlib.import_module(module_name)
    assert mod is not None


# ---------------------------------------------------------------------------
# CLI wiring — verify all subcommands are registered
# ---------------------------------------------------------------------------


EXPECTED_COMMANDS = ["annotate", "assemble", "check", "db-fetch", "func-annot", "gff3toseq", "status"]


def test_cli_help():
    """eukan --help should exit 0 and list all commands."""
    result = subprocess.run(
        ["eukan", "--help"], capture_output=True, text=True,
    )
    assert result.returncode == 0, f"eukan --help failed: {result.stderr}"
    for cmd in EXPECTED_COMMANDS:
        assert cmd in result.stdout, f"Missing subcommand: {cmd}"


def test_cli_version():
    """eukan --version should print a version string."""
    result = subprocess.run(["eukan", "--version"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "0.2.0" in result.stdout


@pytest.mark.parametrize("subcmd", ["annotate", "assemble", "check", "db-fetch", "func-annot", "gff3toseq", "status"])
def test_subcommand_help(subcmd):
    """Every subcommand should respond to --help."""
    result = subprocess.run(["eukan", subcmd, "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "Usage:" in result.stdout


def test_cli_short_help():
    """eukan -h should work (not just --help)."""
    result = subprocess.run(["eukan", "-h"], capture_output=True, text=True)
    assert result.returncode == 0


# ---------------------------------------------------------------------------
# pyhmmer — verify the HMMER replacement works
# ---------------------------------------------------------------------------


def test_pyhmmer_imports():
    """pyhmmer and its key submodules should be importable."""
    import pyhmmer
    import pyhmmer.easel
    import pyhmmer.plan7
    import pyhmmer.hmmer
    assert hasattr(pyhmmer.hmmer, "phmmer")
    assert hasattr(pyhmmer.hmmer, "hmmscan")
    assert hasattr(pyhmmer.hmmer, "hmmpress")


def test_pyhmmer_alphabet():
    """Should be able to create an amino acid alphabet."""
    import pyhmmer.easel
    alpha = pyhmmer.easel.Alphabet.amino()
    assert alpha is not None


def test_pyhmmer_digital_sequence():
    """Should be able to create a digital sequence."""
    import pyhmmer.easel
    alpha = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name=b"test", sequence="MKFLILLFNILCLFPVLAADNH")
    dseq = seq.digitize(alpha)
    # pyhmmer >=0.12 returns str, older versions return bytes
    assert dseq.name in (b"test", "test")
    assert len(dseq) == 22


def test_pyhmmer_phmmer_minimal():
    """phmmer should run on a trivial input (one query vs one target)."""
    import pyhmmer.easel
    import pyhmmer.hmmer

    alpha = pyhmmer.easel.Alphabet.amino()
    query = pyhmmer.easel.TextSequence(name=b"q1", sequence="MKFLILLFNILCLFPVLAADNH").digitize(alpha)
    target = pyhmmer.easel.TextSequence(name=b"t1", sequence="MKFLILLFNILCLFPVLAADNH").digitize(alpha)

    results = list(pyhmmer.hmmer.phmmer([query], [target], cpus=1))
    assert len(results) == 1
    # Should find itself as a hit
    assert len(results[0]) >= 1


# ---------------------------------------------------------------------------
# Tool registry — verify tools.toml loads
# ---------------------------------------------------------------------------


def test_tool_registry_loads():
    """tools.toml should load and contain expected tools."""
    from eukan.check import get_tools
    tools = get_tools()
    assert len(tools) > 10
    names = {t.name for t in tools}
    assert "augustus" in names
    assert "samtools" in names
    assert "genemark" in names


def test_tool_registry_no_hmmer():
    """HMMER tools should NOT be in the registry (replaced by pyhmmer)."""
    from eukan.check import get_tools
    tools = get_tools()
    binaries = {t.binary for t in tools}
    assert "hmmscan" not in binaries
    assert "phmmer" not in binaries
    assert "hmmpress" not in binaries


# ---------------------------------------------------------------------------
# Settings — verify pydantic models work
# ---------------------------------------------------------------------------


def test_pipeline_config_from_env(tmp_path, monkeypatch):
    """PipelineConfig should read EUKAN_ env vars."""
    from eukan.settings import PipelineConfig

    genome = tmp_path / "g.fa"
    genome.touch()
    monkeypatch.setenv("EUKAN_NUM_CPU", "2")

    config = PipelineConfig(genome=genome, proteins=[genome])
    assert config.num_cpu == 2


def test_functional_config_defaults(tmp_path):
    """FunctionalConfig should have sensible database defaults."""
    from eukan.settings import FunctionalConfig

    proteins = tmp_path / "p.faa"
    proteins.touch()
    config = FunctionalConfig(proteins=proteins)
    assert "databases" in str(config.uniprot_db)
    assert "databases" in str(config.pfam_db)


# ---------------------------------------------------------------------------
# Environment file generation
# ---------------------------------------------------------------------------


def test_generate_env_no_hmmer():
    """Generated environment.yml should not include hmmer."""
    from eukan.check import generate_environment_yml
    content = generate_environment_yml()
    assert "hmmer" not in content
    assert "bioconda" in content
    assert "augustus" in content
