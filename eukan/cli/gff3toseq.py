"""eukan gff3toseq — extract protein/cDNA sequences from GFF3."""

from __future__ import annotations

from pathlib import Path
from typing import IO

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    FULL_CODE_TABLE,
    PreformattedEpilogCommand,
    code_option,
    genome_option,
)


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@optgroup.group("Required input")
@genome_option("Genome assembly in FASTA format.")
@optgroup.option(
    "--gff3", "-i", required=True, type=click.Path(exists=True, path_type=Path),
    help="GFF3 file with gene models.",
)
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--output-format", type=click.Choice(["protein", "cdna"], case_sensitive=False),
    default="protein", show_default=True, help="Output sequence type.",
)
@code_option(default=1)
@optgroup.group("Output options")
@optgroup.option(
    "--output-file", "-o", type=click.File("w"), default="-",
    show_default="stdout", help="Write FASTA to this file ('-' for stdout).",
)
def gff3toseq(
    genome: Path, gff3: Path, output_format: str, output_file: IO[str], genetic_code: int,
) -> None:
    """Extract protein or cDNA sequences from GFF3 + genome."""
    from eukan.gff.io import extract_sequences  # lazy: startup cost

    for record in extract_sequences(gff3, genome, extract_to=output_format, genetic_code=genetic_code):
        output_file.write(record.format("fasta"))
