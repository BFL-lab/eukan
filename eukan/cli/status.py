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


def format_status(manifest: RunManifest) -> str:
    """Format the manifest as a human-readable status report."""
    from eukan.infra.manifest import PipelineName  # lazy: startup cost

    lines = [
        f"Run: {manifest.status}",
        f"Started: {manifest.started_at}",
    ]
    if manifest.finished_at:
        lines.append(f"Finished: {manifest.finished_at}")
    lines += [
        f"Genome: {manifest.genome}",
        f"Kingdom: {manifest.kingdom or 'not set'}",
        f"CPUs: {manifest.num_cpu}",
        "",
    ]

    if manifest.tool_versions:
        lines.append("Tool versions:")
        for tool, ver in manifest.tool_versions.items():
            lines.append(f"  {tool:<20s} {ver}")
        lines.append("")

    # Per-pipeline summary above the step detail.
    pipeline_summaries = [
        (p, pipeline_status(manifest, p.value)) for p in PipelineName
    ]
    if any(status != "pending" for _, status in pipeline_summaries):
        lines.append("Pipelines:")
        for pipeline, status in pipeline_summaries:
            if status == "pending":
                continue
            lines.append(f"  {pipeline.value:<12s} {status}")
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
