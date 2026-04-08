"""Pre-flight checks for external tool availability.

Reads the tool registry from tools.toml (repo root) and verifies that
each tool is installed, on PATH, and responds to a version/help probe.
Also checks database integrity via the manifest.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path


# ---------------------------------------------------------------------------
# Tool registry (loaded from tools.toml)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Tool:
    """An external tool dependency."""

    name: str
    binary: str
    version_cmd: list[str]
    required_by: list[str]
    env_var: str | None = None
    conda_package: str | None = None
    min_version: str | None = None
    install_hint: str | None = None


def _find_tools_toml() -> Path | None:
    """Locate tools.toml by checking multiple locations."""
    candidates = [
        Path(__file__).resolve().parent.parent / "tools.toml",  # relative to package
        Path.cwd() / "tools.toml",                               # current directory
    ]
    # Also check EUKAN_ROOT env var (set in Docker)
    eukan_root = os.environ.get("EUKAN_ROOT")
    if eukan_root:
        candidates.insert(0, Path(eukan_root) / "tools.toml")

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_tool_registry() -> list[Tool]:
    """Load the tool registry from tools.toml.

    Falls back to an empty list if the file is not found.
    """
    toml_path = _find_tools_toml()
    if not toml_path:
        return []

    if sys.version_info >= (3, 11):
        import tomllib
    else:
        try:
            import tomllib
        except ImportError:
            import tomli as tomllib  # type: ignore[no-redef]

    with open(toml_path, "rb") as f:
        data = tomllib.load(f)

    tools: list[Tool] = []
    for name, cfg in data.items():
        if not isinstance(cfg, dict) or "binary" not in cfg:
            continue
        tools.append(
            Tool(
                name=name,
                binary=cfg["binary"],
                version_cmd=cfg.get("version_cmd", [cfg["binary"]]),
                required_by=cfg.get("required_by", []),
                env_var=cfg.get("env_var"),
                conda_package=cfg.get("conda_package"),
                min_version=cfg.get("min_version"),
                install_hint=cfg.get("install_hint"),
            )
        )
    return tools


@lru_cache(maxsize=1)
def get_tools() -> tuple[Tool, ...]:
    """Load and cache the tool registry. Thread-safe via lru_cache."""
    return tuple(load_tool_registry())


# ---------------------------------------------------------------------------
# Check logic
# ---------------------------------------------------------------------------


@dataclass
class CheckResult:
    """Result of checking a single tool."""

    tool: Tool
    found: bool
    on_path: str | None
    version_ok: bool
    version_output: str
    env_ok: bool


def check_tool(tool: Tool) -> CheckResult:
    """Check if a tool is available and responds to its version command."""
    binary_path = shutil.which(tool.binary)

    if not binary_path:
        return CheckResult(
            tool=tool,
            found=False,
            on_path=None,
            version_ok=False,
            version_output="",
            env_ok=_check_env(tool),
        )

    try:
        result = subprocess.run(
            tool.version_cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        # Many bioinformatics tools exit non-zero on --help or with no args
        # but still produce useful output. Accept non-zero if output doesn't
        # indicate a broken install (missing shared libraries, etc.)
        output = (result.stdout + result.stderr).strip()
        broken_indicators = ["cannot open shared object", "error while loading shared libraries"]
        is_broken = any(ind in output for ind in broken_indicators)
        version_ok = len(output) > 0 and not is_broken
    except FileNotFoundError:
        return CheckResult(
            tool=tool, found=False, on_path=binary_path,
            version_ok=False, version_output="binary not executable",
            env_ok=_check_env(tool),
        )
    except subprocess.TimeoutExpired:
        output = "timed out"
        version_ok = False
    except Exception as e:
        output = str(e)
        version_ok = False

    first_line = output.split("\n")[0][:120] if output else ""

    return CheckResult(
        tool=tool,
        found=True,
        on_path=binary_path,
        version_ok=version_ok,
        version_output=first_line,
        env_ok=_check_env(tool),
    )


def _check_env(tool: Tool) -> bool:
    """Check if the tool's required environment variable is set."""
    if tool.env_var is None:
        return True
    return tool.env_var in os.environ


# ---------------------------------------------------------------------------
# Python dependency checks
# ---------------------------------------------------------------------------


@dataclass
class PythonCheckResult:
    """Result of checking a Python dependency."""

    name: str
    ok: bool
    detail: str


def check_python_deps() -> list[PythonCheckResult]:
    """Verify that Python dependencies are importable and functional."""
    results: list[PythonCheckResult] = []

    # Module imports
    modules = [
        "eukan", "eukan.cli", "eukan.settings", "eukan.check",
        "eukan.infra.runner", "eukan.infra.manifest", "eukan.infra.steps", "eukan.infra.logging", "eukan.infra.environ",
        "eukan.annotation", "eukan.annotation.orchestrator",
        "eukan.assembly", "eukan.assembly.orchestrator",
        "eukan.functional", "eukan.functional.orchestrator", "eukan.functional.dbfetch",
        "eukan.gff.parser", "eukan.gff.intersecter", "eukan.gff.io",
    ]
    import_failures = []
    for mod in modules:
        try:
            __import__(mod)
        except Exception as e:
            import_failures.append(f"{mod}: {e}")
    if import_failures:
        results.append(PythonCheckResult(
            "module imports", False, "; ".join(import_failures)
        ))
    else:
        results.append(PythonCheckResult(
            "module imports", True, f"all {len(modules)} modules OK"
        ))

    # pyhmmer functional check
    try:
        import pyhmmer.easel
        import pyhmmer.hmmer
        alpha = pyhmmer.easel.Alphabet.amino()
        seq = pyhmmer.easel.TextSequence(
            name=b"test", sequence="MKFLILLFNILCLFPVLAADNH"
        ).digitize(alpha)
        hits = list(pyhmmer.hmmer.phmmer([seq], [seq], cpus=1))
        assert len(hits) == 1 and len(hits[0]) >= 1
        results.append(PythonCheckResult("pyhmmer", True, "phmmer search works"))
    except Exception as e:
        results.append(PythonCheckResult("pyhmmer", False, str(e)))

    # gffutils functional check
    try:
        import gffutils
        gff = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1"
        db = gffutils.create_db(gff, ":memory:", from_string=True)
        genes = list(db.features_of_type("gene"))
        assert len(genes) == 1
        results.append(PythonCheckResult("gffutils", True, "in-memory DB works"))
    except Exception as e:
        results.append(PythonCheckResult("gffutils", False, str(e)))

    # BioPython functional check
    try:
        import io
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        buf = io.StringIO()
        SeqIO.write(record, buf, "fasta")
        buf.seek(0)
        parsed = list(SeqIO.parse(buf, "fasta"))
        assert len(parsed) == 1
        results.append(PythonCheckResult("biopython", True, "sequence I/O works"))
    except Exception as e:
        results.append(PythonCheckResult("biopython", False, str(e)))

    # pydantic-settings functional check
    try:
        import tempfile
        from eukan.settings import PipelineConfig
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            config = PipelineConfig(genome=Path(f.name), proteins=[Path(f.name)])
            assert config.num_cpu >= 1
        results.append(PythonCheckResult("pydantic-settings", True, "config loads OK"))
    except Exception as e:
        results.append(PythonCheckResult("pydantic-settings", False, str(e)))

    # Tool registry
    try:
        tools = get_tools()
        assert len(tools) > 10
        results.append(PythonCheckResult("tool registry", True, f"{len(tools)} tools loaded from tools.toml"))
    except Exception as e:
        results.append(PythonCheckResult("tool registry", False, str(e)))

    return results


# ---------------------------------------------------------------------------
# Main check orchestration
# ---------------------------------------------------------------------------


def run_checks(
    subcommands: list[str] | None = None,
    db_dir: Path | None = None,
) -> tuple[list[CheckResult], list[CheckResult], list[tuple[str, str, bool]], list[PythonCheckResult]]:
    """Run all checks: Python deps, external tools, and databases.

    Returns:
        Tuple of (passed_tools, failed_tools, db_results, python_results).
    """
    # Python dependency checks (always run)
    python_results = check_python_deps()

    # External tool checks
    tools = get_tools()
    if subcommands:
        tools = [
            t for t in tools
            if any(s in t.required_by for s in subcommands)
        ]

    passed: list[CheckResult] = []
    failed: list[CheckResult] = []

    for tool in tools:
        result = check_tool(tool)
        if result.found and result.version_ok and result.env_ok:
            passed.append(result)
        else:
            failed.append(result)

    # Database checks
    db_results: list[tuple[str, str, bool]] = []
    db_relevant = not subcommands or any(
        s in ("func-annot", "annotate", "db-fetch") for s in subcommands
    )
    if db_relevant:
        from eukan.functional.dbfetch import check_databases
        db_dir = db_dir or Path("databases")
        db_results = check_databases(db_dir)

    return passed, failed, db_results, python_results


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def format_results(
    passed: list[CheckResult],
    failed: list[CheckResult],
    db_results: list[tuple[str, str, bool]] | None = None,
    python_results: list[PythonCheckResult] | None = None,
) -> str:
    """Format check results as a human-readable report."""
    lines: list[str] = []

    # --- Python dependencies ---
    if python_results:
        py_ok = [r for r in python_results if r.ok]
        py_fail = [r for r in python_results if not r.ok]

        if py_ok:
            lines.append(f"  {len(py_ok)} Python checks OK:")
            for r in py_ok:
                lines.append(f"    \u2713 {r.name:<30s} {r.detail}")
        if py_fail:
            lines.append(f"\n  {len(py_fail)} Python checks FAILED:")
            for r in py_fail:
                lines.append(f"    \u2717 {r.name:<30s} {r.detail}")
        lines.append("")

    # --- External tools ---
    if passed:
        lines.append(f"  {len(passed)} tools OK:")
        for r in passed:
            lines.append(f"    \u2713 {r.tool.name:<30s} {r.version_output}")

    if failed:
        lines.append(f"\n  {len(failed)} tools MISSING or BROKEN:")
        for r in failed:
            issues = []
            if not r.found:
                issues.append(f"`{r.tool.binary}` not found on PATH")
            elif not r.version_ok:
                issues.append(f"version check failed: {r.version_output or 'no output'}")
            if not r.env_ok:
                issues.append(f"${r.tool.env_var} not set")
            lines.append(f"    \u2717 {r.tool.name:<30s} {'; '.join(issues)}")
            lines.append(f"      used by: {', '.join(r.tool.required_by)}")
            if r.tool.install_hint and not r.found:
                lines.append(f"      hint: {r.tool.install_hint}")

    total = len(passed) + len(failed)
    lines.append(f"\n  Checked {total} external tools total.")

    # --- Databases ---
    if db_results:
        lines.append("")
        db_ok = [r for r in db_results if r[2]]
        db_fail = [r for r in db_results if not r[2]]

        if db_ok:
            lines.append(f"  {len(db_ok)} databases OK:")
            for name, msg, _ in db_ok:
                lines.append(f"    \u2713 {name:<30s} {msg}")
        if db_fail:
            lines.append(f"\n  {len(db_fail)} databases MISSING or INVALID:")
            for name, msg, _ in db_fail:
                lines.append(f"    \u2717 {name:<30s} {msg}")
            lines.append("      run: eukan db-fetch")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Environment file generation
# ---------------------------------------------------------------------------


def generate_environment_yml() -> str:
    """Generate a conda environment.yml from the tool registry.

    Returns the YAML content as a string.
    """
    tools = get_tools()

    # Deduplicate conda packages (some tools share a package, e.g. hmmer)
    seen: set[str] = set()
    conda_deps: list[str] = []
    for tool in tools:
        if not tool.conda_package or tool.conda_package in seen:
            continue
        seen.add(tool.conda_package)
        pin = f">={tool.min_version}" if tool.min_version else ""
        conda_deps.append(f"  - {tool.conda_package}{pin}")

    # Tools requiring manual install
    manual: list[str] = []
    for tool in tools:
        if not tool.conda_package and tool.install_hint:
            manual.append(f"#   - {tool.name}: {tool.install_hint}")

    # Bundled tools (no conda package, no install hint)
    bundled: list[str] = []
    for tool in tools:
        if not tool.conda_package and not tool.install_hint:
            bundled.append(f"#   - {tool.binary}: bundled with parent package")

    lines = [
        "# Auto-generated from tools.toml — do not edit manually.",
        "# Regenerate with: eukan dev generate-env",
        "#",
        "# Tools requiring manual installation:",
        *manual,
        "#",
        "# Tools bundled with parent packages (no separate install):",
        *bundled,
        "",
        "name: eukan",
        "channels:",
        "  - bioconda",
        "  - conda-forge",
        "  - defaults",
        "",
        "dependencies:",
        "  # Python",
        "  - python>=3.10,<3.13",
        "",
        "  # External tools",
        *sorted(conda_deps),
        "",
        "  # Build dependencies (needed to compile fitild from source)",
        "  - gsl",
        "  - liblbfgs",
        "",
        "  # Perl (needed by AUGUSTUS, PASA, EVM, GeneMark helper scripts)",
        "  - perl",
        "  - perl-bioperl",
        "  - perl-dbd-sqlite",
        "  - perl-yaml",
        "  - perl-hash-merge",
        "  - perl-mce",
        "  - perl-parallel-forkmanager",
        "",
        "  # eukan + Python dependencies (installed via pip for version consistency)",
        "  - pip",
        "  - pip:",
        "    - .",
        "",
    ]

    return "\n".join(lines) + "\n"
