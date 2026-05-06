"""eukan gff3toseq — extract protein/cDNA sequences from GFF3."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from eukan.cli._framework import FULL_CODE_TABLE, PreformattedEpilogCommand


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@click.option(
    "--genome", "-g", required=True, type=click.Path(exists=True, path_type=Path),
    help="Genome assembly in FASTA format.",
)
@click.option(
    "--gff3", "-i", required=True, type=click.Path(exists=True, path_type=Path),
    help="GFF3 file with gene models.",
)
@click.option(
    "--output", "-o", type=click.Choice(["protein", "cdna"]), default="protein",
    show_default=True, help="Output sequence type.",
)
@click.option(
    "--code", "-c", type=int, default=1, show_default=True,
    help="NCBI genetic code table number.",
)
def gff3toseq(genome: Path, gff3: Path, output: str, code: int) -> None:
    """Extract protein or cDNA sequences from GFF3 + genome."""
    from eukan.gff.io import extract_sequences

    for record in extract_sequences(gff3, genome, extract_to=output, genetic_code=code):
        sys.stdout.write(record.format("fasta"))
