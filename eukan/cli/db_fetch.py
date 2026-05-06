"""eukan db-fetch — download UniProt + Pfam reference databases."""

from __future__ import annotations

from pathlib import Path

import click


@click.command("db-fetch")
@click.option(
    "--output-dir", "-o", type=click.Path(path_type=Path), default="databases",
    show_default=True, help="Directory to download databases into.",
)
@click.option(
    "--force", "-f", is_flag=True, help="Re-download even if databases are up to date.",
)
@click.option(
    "--database", "-d", multiple=True,
    type=click.Choice(["uniprot", "pfam"], case_sensitive=False),
    help="Specific database(s) to fetch. If omitted, fetch all.",
)
def db_fetch(output_dir: Path, force: bool, database: tuple[str, ...]) -> None:
    """Download reference databases (UniProt, Pfam)."""
    from eukan.functional.dbfetch import fetch_databases

    output_dir = output_dir.resolve()
    fetch_databases(
        output_dir,
        force=force,
        databases=list(database) if database else None,
    )
    click.echo("Done.")
