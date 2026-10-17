"""eukan check — verify Python deps, external tools, and databases."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import click
from click_option_group import optgroup

from eukan.cli._framework import PreformattedEpilogCommand
from eukan.infra.databases import DEFAULT_DB_DIR, HOMOLOGY_DBS

if TYPE_CHECKING:
    from eukan.infra.health import CheckResult, PythonCheckResult

# Shown under the CPU baseline when a tool died of SIGILL on a pre-v3 CPU.
_SIGILL_ADVICE = (
    "    SIGILL crashes typically mean the binary requires a newer\n"
    "    CPU baseline than this machine provides. Bioconda ships\n"
    "    x86-64-v3 (AVX2/FMA/BMI2) builds for many tools. To fix:\n"
    "    rebuild from source on this host, run inside Docker on a\n"
    "    newer CPU, or pin to an older bioconda build."
)


def _section(what: str, broken: str, rows: list[tuple[str, str, bool, list[str]]]) -> list[str]:
    """Render ``(name, detail, ok, notes)`` rows as two lists.

    The passes go under "<n> <what> OK:", the failures, each followed by its
    note lines, under "<n> <what> <broken>:".
    """
    passes = [row for row in rows if row[2]]
    failures = [row for row in rows if not row[2]]
    lines: list[str] = []
    if passes:
        lines.append(f"  {len(passes)} {what} OK:")
        lines += [f"    \u2713 {name:<30s} {detail}" for name, detail, _, _ in passes]
    if failures:
        lines.append(f"\n  {len(failures)} {what} {broken}:")
        for name, detail, _, notes in failures:
            lines.append(f"    \u2717 {name:<30s} {detail}")
            lines += notes
    return lines


def _tool_notes(result: CheckResult) -> list[str]:
    notes = [f"      used by: {', '.join(result.tool.required_by)}"]
    if result.tool.install_hint and not result.found:
        notes.append(f"      hint: {result.tool.install_hint}")
    return notes


def format_results(
    passed: list[CheckResult],
    failed: list[CheckResult],
    db_results: list[tuple[str, str, bool]] | None = None,
    python_results: list[PythonCheckResult] | None = None,
) -> str:
    """Format check results as a human-readable report."""
    from eukan.infra.health import _FATAL_SIGNALS, cpu_baseline  # lazy: startup cost

    lines: list[str] = []
    if python_results:
        lines += _section(
            "Python checks", "FAILED", [(r.name, r.detail, r.ok, []) for r in python_results],
        )
        lines.append("")

    lines += _section("tools", "MISSING or BROKEN", [
        *((r.tool.name, r.version_output, True, []) for r in passed),
        *((r.tool.name, "; ".join(r.problems), False, _tool_notes(r)) for r in failed),
    ])
    lines.append(f"\n  Checked {len(passed) + len(failed)} external tools total.")

    # A tool killed by a fatal signal is usually a binary built for a newer
    # CPU than this one (bioconda ships x86-64-v3, Haswell+, builds), so show
    # the CPU baseline next to it.
    crashed = [r for r in failed if r.crash_signal in _FATAL_SIGNALS]
    baseline = cpu_baseline() if crashed else None
    if baseline is not None:
        level, _ = baseline
        lines += ["", f"  CPU baseline: {level}"]
        if level in ("x86-64", "x86-64-v2") and any(r.crash_signal == "SIGILL" for r in crashed):
            lines.append(_SIGILL_ADVICE)

    if db_results:
        lines.append("")
        lines += _section(
            "databases", "MISSING or INVALID", [(name, msg, ok, []) for name, msg, ok in db_results],
        )
        if not all(ok for _, _, ok in db_results):
            lines.append("      run: eukan db-fetch")

    return "\n".join(lines)


@click.command(cls=PreformattedEpilogCommand)
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--for", "subcommands", multiple=True,
    type=click.Choice(
        ["annotate", "assemble", "func-annot", "db-fetch", "mask-repeats", "prep-submission"],
        case_sensitive=False,
    ),
    help="Only check tools needed by these subcommands. If omitted, check all.",
)
@optgroup.option(
    "--db-dir", type=click.Path(resolve_path=True, path_type=Path), default=str(DEFAULT_DB_DIR),
    show_default=True, help="Database directory to check.",
)
@optgroup.option(
    "--homology-db", type=click.Choice(HOMOLOGY_DBS, case_sensitive=False),
    default=None,
    help="When set, only check the chosen homology DB (plus Pfam). Without "
         "it, every registered database is checked.",
)
def check(
    subcommands: tuple[str, ...], db_dir: Path, homology_db: str | None,
) -> None:
    """Verify Python deps, external tools, and databases."""
    from eukan.infra.health import run_checks  # lazy: startup cost

    passed, failed, db_results, python_results = run_checks(
        list(subcommands) if subcommands else None,
        db_dir=db_dir,
        homology_db=homology_db,
    )
    click.echo(format_results(passed, failed, db_results, python_results))

    py_failures = any(not r.ok for r in python_results)
    db_failures = any(not ok for _, _, ok in db_results)
    if failed or db_failures or py_failures:
        raise SystemExit(1)
