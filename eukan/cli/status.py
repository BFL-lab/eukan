"""eukan status — show pipeline run progress from eukan-run.json."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import click

if TYPE_CHECKING:
    from eukan.infra.manifest import RunManifest

# Keyed by StepStatus value.
_STATUS_ICONS = {
    "completed": "✓",
    "running": "▶",
    "failed": "✗",
    "pending": "·",
    "skipped": "-",
}


def pipeline_status(manifest: RunManifest, prefix: str) -> str:
    """Aggregate status of the steps keyed ``<prefix>/...``.

    'pending' when there are none; else 'failed' if any failed, 'running' if
    any is running, 'completed' if all completed, and 'running' otherwise.
    """
    statuses = {r.status.value for k, r in manifest.steps.items() if k.startswith(prefix + "/")}
    if not statuses:
        return "pending"
    if "failed" in statuses:
        return "failed"
    if "running" in statuses:
        return "running"
    if statuses == {"completed"}:
        return "completed"
    return "running"


def _tool_version_lines(versions: dict[str, str], indent: str) -> list[str]:
    if not versions:
        return []
    return [f"{indent}Tool versions:"] + [
        f"{indent}  {tool:<20s} {ver}" for tool, ver in versions.items()
    ]


def format_status(manifest: RunManifest) -> str:
    """Format the manifest as a human-readable status report: a block per
    pipeline that has run, then the steps."""
    from eukan.infra.manifest import PipelineName  # lazy: startup cost

    lines: list[str] = []
    if manifest.started_at is not None:
        # A version 1 manifest's run-level fields: those of whichever pipeline
        # created it.
        lines += [f"Run: {manifest.status}", f"Started: {manifest.started_at}"]
        if manifest.finished_at:
            lines.append(f"Finished: {manifest.finished_at}")
        lines += [
            f"Genome: {manifest.genome}",
            f"Kingdom: {manifest.kingdom or 'not set'}",
            f"CPUs: {manifest.num_cpu}",
            "",
        ]
        lines += _tool_version_lines(manifest.tool_versions or {}, "")
        lines += [""] if manifest.tool_versions else []

    for pipeline in PipelineName:
        run = manifest.pipelines.get(pipeline)
        if run is None:
            # No record (it ran before version 2): status from its steps alone.
            status = pipeline_status(manifest, pipeline.value)
            if status != "pending":
                lines += [f"{pipeline.value}: {status}", ""]
            continue
        lines.append(f"{pipeline.value}: {run.status}")
        transcripts = run.has_transcripts
        for label, value in (
            ("Started", run.started_at),
            ("Finished", run.finished_at),
            ("Genome", run.genome),
            ("Proteins", ", ".join(map(str, run.proteins)) if run.proteins else None),
            ("Kingdom", run.kingdom),
            ("Genetic code", run.genetic_code),
            ("Transcript evidence", None if transcripts is None else "yes" if transcripts else "no"),
            ("CPUs", run.num_cpu),
        ):
            if value is not None:
                lines.append(f"  {label}: {value}")
        lines += _tool_version_lines(run.tool_versions, "  ")
        lines.append("")

    lines.append("Steps:")
    for step_name, record in manifest.steps.items():
        icon = _STATUS_ICONS.get(record.status.value, "?")
        duration = f" ({record.duration_seconds}s)" if record.duration_seconds else ""
        lines.append(f"  {icon} {step_name:<25s} {record.status.value}{duration}")
        if record.error:
            lines.append(f"    error: {record.error[:200]}")
        if record.output_md5:
            lines.append(f"    output: {record.output_file} (md5:{record.output_md5[:12]}...)")

    return "\n".join(lines)


@click.command()
@click.option(
    "--work-dir", "-d", type=click.Path(exists=True, resolve_path=True, path_type=Path), default=".",
    show_default=True, help="Working directory containing eukan-run.json.",
)
def status(work_dir: Path) -> None:
    """Show the status of a pipeline run."""
    from eukan.infra.manifest import load_manifest  # lazy: startup cost

    manifest = load_manifest(work_dir)
    if not manifest:
        click.echo("No pipeline run found in this directory (no eukan-run.json).")
        raise SystemExit(1)

    click.echo(format_status(manifest))
