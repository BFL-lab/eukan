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
@click.option(
    "--homology-db", type=click.Choice(["uniprot", "kofam"], case_sensitive=False),
    default=None,
    help="When set, only check the chosen homology DB (plus Pfam). Without "
         "it, every registered database is checked.",
)
def check(
    subcommands: tuple[str, ...], db_dir: Path, homology_db: str | None,
) -> None:
    """Verify Python deps, external tools, and databases."""
    from eukan.infra.health import format_results, run_checks

    passed, failed, db_results, python_results = run_checks(
        list(subcommands) if subcommands else None,
        db_dir=db_dir.resolve(),
        homology_db=homology_db.lower() if homology_db else None,
    )
    click.echo(format_results(passed, failed, db_results, python_results))

    py_failures = any(not r.ok for r in python_results)
    db_failures = any(not ok for _, _, ok in db_results)
    if failed or db_failures or py_failures:
        raise SystemExit(1)
