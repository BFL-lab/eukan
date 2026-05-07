"""eukan check — verify Python deps, external tools, and databases."""

from __future__ import annotations

from pathlib import Path

import click


@click.command()
@click.option(
    "--for", "subcommands", multiple=True,
    type=click.Choice(
        ["annotate", "assemble", "func-annot", "db-fetch", "mask-repeats", "prep-submission"],
        case_sensitive=False,
    ),
    help="Only check tools needed by these subcommands. If omitted, check all.",
)
@click.option(
    "--db-dir", type=click.Path(path_type=Path), default="databases",
    show_default=True, help="Database directory to check.",
)
def check(subcommands: tuple[str, ...], db_dir: Path) -> None:
    """Verify Python deps, external tools, and databases."""
    from eukan.infra.health import format_results, run_checks

    passed, failed, db_results, python_results = run_checks(
        list(subcommands) if subcommands else None,
        db_dir=db_dir.resolve(),
    )
    click.echo(format_results(passed, failed, db_results, python_results))

    py_failures = any(not r.ok for r in python_results)
    db_failures = any(not ok for _, _, ok in db_results)
    if failed or db_failures or py_failures:
        raise SystemExit(1)
