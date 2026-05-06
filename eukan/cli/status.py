"""eukan status — show pipeline run progress from eukan-run.json."""

from __future__ import annotations

from pathlib import Path

import click


@click.command()
@click.option(
    "--work-dir", "-d", type=click.Path(exists=True, path_type=Path), default=".",
    help="Working directory containing eukan-run.json.",
)
def status(work_dir: Path) -> None:
    """Show the status of a pipeline run."""
    from eukan.infra.manifest import format_status, load_manifest

    manifest = load_manifest(work_dir.resolve())
    if not manifest:
        click.echo("No pipeline run found in this directory (no eukan-run.json).")
        raise SystemExit(1)

    click.echo(format_status(manifest))
