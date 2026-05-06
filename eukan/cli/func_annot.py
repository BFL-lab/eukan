"""eukan func-annot — UniProt + Pfam functional annotation."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PreformattedEpilogCommand,
    drop_none,
    force_option,
    numcpu_option,
)


@click.command("func-annot", cls=PreformattedEpilogCommand)
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option("--evalue", "-e", type=str, default="1e-1", show_default=True, help="E-value cutoff.")
@optgroup.group("Override options")
@optgroup.option(
    "--proteins", "-p", type=click.Path(exists=True, path_type=Path),
    help="Amino acid sequences in FASTA format.",
)
@optgroup.option(
    "--uniprot", type=click.Path(exists=True, path_type=Path),
    default=None, help="UniProt-SwissProt database FASTA.",
)
@optgroup.option(
    "--pfam", type=click.Path(exists=True, path_type=Path),
    default=None, help="Pfam HMM database.",
)
@optgroup.option(
    "--gff3", type=click.Path(exists=True, path_type=Path),
    default=None, help="GFF3 file to annotate with functional info.",
)
@force_option
def func_annot(
    proteins: Path,
    uniprot: Path | None,
    pfam: Path | None,
    gff3: Path | None,
    numcpu: int,
    evalue: str,
    force: bool,
) -> None:
    """Add functional annotations (UniProt + Pfam) to proteins.

    \b
    When run after `eukan annotate` and `eukan db-fetch`, the predicted
    protein sequences, UniProt, and Pfam databases are discovered
    automatically. Use the override options to point to different files
    or to run functional annotation independently of the main pipeline.
    """
    from eukan.functional import run_functional_annotation
    from eukan.settings import FunctionalConfig

    if proteins is None:
        raise click.UsageError(
            "No protein file found. Provide --proteins or run `eukan annotate` first."
        )

    config = FunctionalConfig(**drop_none(
        proteins=proteins.resolve(),
        work_dir=Path.cwd(),
        num_cpu=numcpu,
        evalue=evalue,
        uniprot_db=uniprot.resolve() if uniprot else None,
        pfam_db=pfam.resolve() if pfam else None,
        gff3_path=gff3.resolve() if gff3 else None,
    ))
    run_functional_annotation(config, force=force)
    click.echo("Done.")
