"""eukan db-fetch — download reference databases (UniProt or KOfam, plus Pfam)."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import force_option


@click.command("db-fetch")
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--output-dir", "-o", type=click.Path(resolve_path=True, path_type=Path), default="databases",
    show_default=True, help="Directory to download databases into.",
)
@optgroup.option(
    "--homology-db", type=click.Choice(["uniprot", "kofam"], case_sensitive=False),
    default="uniprot", show_default=True,
    help="Which homology DB to fetch alongside Pfam. 'uniprot' downloads "
         "SwissProt; 'kofam' downloads the KOfam HMM profiles + ko_list "
         "and presses an eukaryote-only HMM database.",
)
@optgroup.option(
    "--database", "-d", multiple=True,
    type=click.Choice(["uniprot", "pfam", "kofam", "ko_list"], case_sensitive=False),
    help="Specific database(s) to fetch. Overrides --homology-db when given.",
)
@force_option(help_text="Re-download even if databases are up to date.")
def db_fetch(
    output_dir: Path, homology_db: str, force: bool, database: tuple[str, ...],
) -> None:
    """Download reference databases (UniProt or KOfam, plus Pfam).

    \b
    Without --database, fetches Pfam plus the homology DB selected by
    --homology-db:
      uniprot  → uniprot_sprot.faa  (default; SwissProt phmmer search)
      kofam    → kofam_eukaryote.hmm + ko_list.tsv  (KEGG Orthology HMMs)

    \b
    Use --database to fetch an explicit subset (e.g. -d pfam to refresh
    only Pfam, or -d kofam -d ko_list to fetch KOfam without touching Pfam).
    """
    from eukan.functional.dbfetch import fetch_databases

    fetch_databases(
        output_dir,
        force=force,
        databases=list(database) if database else None,
        homology_db=homology_db.lower(),
    )
    click.echo("Done.")
