"""Pre-flight checks for external tool availability.

Reads the tool registry from tools.toml (via :mod:`tools_registry`) and
verifies that each tool is installed, on PATH, and responds to a
version/help probe. Also checks database integrity via the manifest.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from eukan.infra.tools_registry import Tool, load_tools

# ---------------------------------------------------------------------------
# Check logic
# ---------------------------------------------------------------------------

# Signals that indicate the binary was killed mid-run rather than exited
# normally — almost always a broken install (CPU baseline mismatch,
# corrupted binary, missing kernel feature). See _crash_signal().
_FATAL_SIGNALS = {"SIGILL", "SIGSEGV", "SIGBUS", "SIGABRT", "SIGFPE"}

# Substrings that bash/zsh emit when a child of a wrapper script is killed by
# a fatal signal, mapped to the signal name. The wrapper's own exit code is
# 128+N (signal-encoded, handled by _crash_signal), but its stderr also carries
# these human-readable strings — so wrappers like /bin/STAR (which pick a SIMD
# variant and exec it) get caught even when the wrapper exits with the encoded
# code.
_CRASH_STRING_TO_SIGNAL = {
    "Illegal instruction": "SIGILL",
    "Segmentation fault": "SIGSEGV",
    "Bus error": "SIGBUS",
}


def _crash_signal(returncode: int) -> str | None:
    """Return signal name (e.g. ``SIGILL``) if the process died from a signal.

    Two encodings to handle:
    - Direct child killed by signal N: Python's ``subprocess.run`` reports
      ``returncode == -N``.
    - Shell wrapper whose child was killed by signal N: the wrapper exits
      with ``128+N`` and Python sees that positive value.
    """
    if returncode < 0:
        sig_num = -returncode
    elif 128 < returncode < 160:  # 128 + 1..31
        sig_num = returncode - 128
    else:
        return None
    try:
        return signal.Signals(sig_num).name
    except ValueError:
        return f"signal {sig_num}"


@dataclass
class CheckResult:
    """Result of checking a single tool."""

    tool: Tool
    found: bool
    on_path: str | None
    version_ok: bool
    version_output: str
    env_ok: bool
    crash_signal: str | None = None


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

    crash_sig: str | None = None
    try:
        # Use subprocess_env() so tools that depend on LD_LIBRARY_PATH
        # (e.g. fitild needs liblbfgs from $CONDA_PREFIX/lib) load their
        # shared libs without requiring the user to export LD_LIBRARY_PATH.
        from eukan.infra.environ import subprocess_env
        result = subprocess.run(
            tool.version_cmd,
            capture_output=True,
            text=True,
            timeout=10,
            stdin=subprocess.DEVNULL,
            env=subprocess_env(),
        )
        output = (result.stdout + result.stderr).strip()

        # Detect fatal-signal exits both directly (returncode<0) and via
        # shell-wrapper encoding (returncode 128+N). Bash also prints the
        # signal name to stderr — catch that too in case the wrapper
        # masks the encoded exit code (e.g. STAR's `for SIMD; do …` loop
        # falls through and the parent shell ends up with $? from a
        # nested call).
        crash_sig = _crash_signal(result.returncode)
        if crash_sig is None:
            for crash_str, sig_name in _CRASH_STRING_TO_SIGNAL.items():
                if crash_str in output:
                    crash_sig = sig_name
                    break

        # Many bioinformatics tools exit non-zero on --help or with no args
        # but still produce useful output. Accept non-zero if output doesn't
        # indicate a broken install (missing shared libraries, fatal
        # signals, etc.).
        broken_indicators = (
            "cannot open shared object",
            "error while loading shared libraries",
        )
        is_broken = (
            crash_sig is not None
            or any(ind in output for ind in broken_indicators)
        )
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

    if crash_sig:
        version_output = f"crashed with {crash_sig}"
    else:
        version_output = output.split("\n")[0][:120] if output else ""

    return CheckResult(
        tool=tool,
        found=True,
        on_path=binary_path,
        version_ok=version_ok,
        version_output=version_output,
        env_ok=_check_env(tool),
        crash_signal=crash_sig,
    )


def _check_env(tool: Tool) -> bool:
    """Check that every env var declared by the tool is set."""
    return all(spec.var in os.environ for spec in tool.env_vars)


def _missing_env_vars(tool: Tool) -> list[str]:
    return [spec.var for spec in tool.env_vars if spec.var not in os.environ]


# ---------------------------------------------------------------------------
# CPU baseline detection — surfaced when binaries crash with SIGILL so the
# user can see at a glance whether their CPU lacks the instruction set the
# (bioconda) build expects.
# ---------------------------------------------------------------------------

# x86-64 microarchitecture levels per the System V psABI Supplement.
# Listed highest first so the first match wins.
_X86_64_LEVELS = (
    ("x86-64-v4", {"avx512f", "avx512bw", "avx512cd", "avx512dq", "avx512vl"}),
    ("x86-64-v3", {"avx2", "bmi1", "bmi2", "f16c", "fma", "abm", "movbe"}),
    ("x86-64-v2", {"sse4_2", "sse4_1", "ssse3", "popcnt", "cx16"}),
)


def cpu_baseline() -> tuple[str, set[str]] | None:
    """Detect the highest x86-64 microarch level supported by the running CPU.

    Returns (level_name, flags_set) on Linux x86_64, or None if the info
    isn't available (non-x86, missing /proc/cpuinfo, etc.). Used purely
    for surfacing diagnostics — the pipeline never gates on this.

    Note: Linux's /proc/cpuinfo names some flags differently from the
    psABI list (e.g. ``abm`` is implied by ``lzcnt`` + ``popcnt``;
    ``movbe`` is reported as-is). We accept reasonable substitutes.
    """
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("flags") and ":" in line:
                    flags = set(line.split(":", 1)[1].split())
                    break
            else:
                return None
    except OSError:
        return None

    # Substitutions: psABI flag → equivalent /proc/cpuinfo flag(s)
    if "lzcnt" in flags:
        flags.add("abm")  # abm = lzcnt + popcnt; popcnt checked separately

    for level, required in _X86_64_LEVELS:
        if required.issubset(flags):
            return level, flags
    return "x86-64", flags  # baseline always satisfied on amd64


# ---------------------------------------------------------------------------
# Python dependency checks
# ---------------------------------------------------------------------------


@dataclass
class PythonCheckResult:
    """Result of checking a Python dependency."""

    name: str
    ok: bool
    detail: str


def _check_python_dep(name: str, probe: Callable[[], str]) -> PythonCheckResult:
    """Run one dependency probe, turning any exception into a failed result.

    ``probe`` performs the functional check and returns the success detail
    string (or raises). This keeps each check to its happy-path body plus a
    returned detail, with the try/except handled once here.
    """
    try:
        return PythonCheckResult(name, True, probe())
    except Exception as e:
        return PythonCheckResult(name, False, str(e))


def _check_module_imports() -> PythonCheckResult:
    """Import every first-party module, collecting any failures."""
    modules = [
        "eukan", "eukan.cli", "eukan.settings", "eukan.infra.health",
        "eukan.infra.runner", "eukan.infra.manifest", "eukan.infra.steps", "eukan.infra.pipeline", "eukan.infra.logging", "eukan.infra.environ",
        "eukan.annotation", "eukan.annotation.pipeline",
        "eukan.assembly", "eukan.assembly.pipeline",
        "eukan.repeats", "eukan.repeats.pipeline", "eukan.repeats.modeler", "eukan.repeats.masker",
        "eukan.functional", "eukan.functional.pipeline", "eukan.functional.dbfetch",
        "eukan.submission",
        "eukan.gff.transforms", "eukan.gff.hierarchy", "eukan.gff.concordance", "eukan.gff.io",
        "eukan.compare", "eukan.compare.engine", "eukan.compare.format",
        "eukan.compare.models",
    ]
    failures = []
    for mod in modules:
        try:
            __import__(mod)
        except Exception as e:
            failures.append(f"{mod}: {e}")
    if failures:
        return PythonCheckResult("module imports", False, "; ".join(failures))
    return PythonCheckResult("module imports", True, f"all {len(modules)} modules OK")


def _probe_pyhmmer() -> str:
    import pyhmmer.easel
    import pyhmmer.hmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(
        name=b"test", sequence="MKFLILLFNILCLFPVLAADNH"
    ).digitize(alpha)
    hits = list(pyhmmer.hmmer.phmmer([seq], [seq], cpus=1))  # type: ignore[list-item]
    assert len(hits) == 1 and len(hits[0]) >= 1
    return "phmmer search works"


def _probe_gffutils() -> str:
    import gffutils
    gff = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1"
    db = gffutils.create_db(gff, ":memory:", from_string=True)
    assert len(list(db.features_of_type("gene"))) == 1
    return "in-memory DB works"


def _probe_biopython() -> str:
    import io

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    record = SeqRecord(Seq("ATGAAATAA"), id="test")
    buf = io.StringIO()
    SeqIO.write(record, buf, "fasta")
    buf.seek(0)
    assert len(list(SeqIO.parse(buf, "fasta"))) == 1
    return "sequence I/O works"


def _probe_pydantic_settings() -> str:
    import tempfile

    from eukan.settings import PipelineConfig
    with tempfile.NamedTemporaryFile(suffix=".fa") as f:
        config = PipelineConfig(genome=Path(f.name), proteins=[Path(f.name)])
        assert config.num_cpu >= 1
    return "config loads OK"


def _probe_scipy() -> str:
    from scipy.stats import (
        chi2_contingency,
        false_discovery_control,
        ks_2samp,
    )
    d, _ = ks_2samp([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    assert 0.0 <= d <= 1.0
    chi2_res = chi2_contingency([[10, 20], [30, 40]])
    assert int(chi2_res.dof) == 1
    adj = false_discovery_control([0.01, 0.5], method="bh")
    assert len(adj) == 2
    return "ks_2samp/chi2/BH work"


def _probe_tool_registry() -> str:
    tools = load_tools()
    assert len(tools) > 10
    return f"{len(tools)} tools loaded from tools.toml"


def check_python_deps() -> list[PythonCheckResult]:
    """Verify that Python dependencies are importable and functional."""
    return [
        _check_module_imports(),
        _check_python_dep("pyhmmer", _probe_pyhmmer),
        _check_python_dep("gffutils", _probe_gffutils),
        _check_python_dep("biopython", _probe_biopython),
        _check_python_dep("pydantic-settings", _probe_pydantic_settings),
        _check_python_dep("scipy", _probe_scipy),
        _check_python_dep("tool registry", _probe_tool_registry),
    ]


# ---------------------------------------------------------------------------
# Main check orchestration
# ---------------------------------------------------------------------------


def run_checks(
    subcommands: list[str] | None = None,
    db_dir: Path | None = None,
    homology_db: str | None = None,
) -> tuple[list[CheckResult], list[CheckResult], list[tuple[str, str, bool]], list[PythonCheckResult]]:
    """Run all checks: Python deps, external tools, and databases.

    Returns:
        Tuple of (passed_tools, failed_tools, db_results, python_results).
    """
    # Python dependency checks (always run)
    python_results = check_python_deps()

    # External tool checks
    tools: tuple[Tool, ...] | list[Tool] = load_tools()
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
        db_results = check_databases(db_dir, homology_db=homology_db)

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
        py_ok = [pr for pr in python_results if pr.ok]
        py_fail = [pr for pr in python_results if not pr.ok]

        if py_ok:
            lines.append(f"  {len(py_ok)} Python checks OK:")
            for pr in py_ok:
                lines.append(f"    \u2713 {pr.name:<30s} {pr.detail}")
        if py_fail:
            lines.append(f"\n  {len(py_fail)} Python checks FAILED:")
            for pr in py_fail:
                lines.append(f"    \u2717 {pr.name:<30s} {pr.detail}")
        lines.append("")

    # --- External tools ---
    if passed:
        lines.append(f"  {len(passed)} tools OK:")
        for tr in passed:
            lines.append(f"    \u2713 {tr.tool.name:<30s} {tr.version_output}")

    if failed:
        lines.append(f"\n  {len(failed)} tools MISSING or BROKEN:")
        for tr in failed:
            issues = []
            if not tr.found:
                issues.append(f"`{tr.tool.binary}` not found on PATH")
            elif tr.crash_signal in _FATAL_SIGNALS:
                issues.append(f"crashed with {tr.crash_signal}")
            elif not tr.version_ok:
                issues.append(f"version check failed: {tr.version_output or 'no output'}")
            if not tr.env_ok:
                missing = ", ".join(f"${v}" for v in _missing_env_vars(tr.tool))
                issues.append(f"env not set: {missing}")
            lines.append(f"    \u2717 {tr.tool.name:<30s} {'; '.join(issues)}")
            lines.append(f"      used by: {', '.join(tr.tool.required_by)}")
            if tr.tool.install_hint and not tr.found:
                lines.append(f"      hint: {tr.tool.install_hint}")

    total = len(passed) + len(failed)
    lines.append(f"\n  Checked {total} external tools total.")

    # When any tool died from a fatal signal, surface the CPU baseline so
    # the user can see whether the binary is built for a newer CPU than
    # they have. SIGILL on the same set of bioconda packages (augustus,
    # spaln, STAR) is almost always a CPU baseline mismatch \u2014 bioconda
    # ships x86-64-v3 (Haswell+) builds.
    crashed = [tr for tr in failed if tr.crash_signal in _FATAL_SIGNALS]
    if crashed:
        baseline = cpu_baseline()
        if baseline is not None:
            level, _ = baseline
            lines.append("")
            lines.append(f"  CPU baseline: {level}")
            sigill_tools = [tr.tool.name for tr in crashed if tr.crash_signal == "SIGILL"]
            if sigill_tools and level in ("x86-64", "x86-64-v2"):
                lines.append(
                    "    SIGILL crashes typically mean the binary requires a newer"
                )
                lines.append(
                    "    CPU baseline than this machine provides. Bioconda ships"
                )
                lines.append(
                    "    x86-64-v3 (AVX2/FMA/BMI2) builds for many tools. To fix:"
                )
                lines.append(
                    "    rebuild from source on this host, run inside Docker on a"
                )
                lines.append(
                    "    newer CPU, or pin to an older bioconda build."
                )

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
