"""Input validation and sanitization for the annotation pipeline."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from eukan.exceptions import FastaValidationError
from eukan.infra.logging import get_logger

log = get_logger(__name__)


def validate_fasta(path: Path) -> None:
    """Verify that a file is a parseable FASTA."""
    with open(path) as handle:
        records = SeqIO.parse(handle, "fasta")
        if not any(records):
            raise FastaValidationError(path, "file is empty or unparseable")


def sanitize_genome_fasta(genome: Path, work_dir: Path) -> Path:
    """Create a copy of the genome with clean FASTA headers.

    Strips descriptions from headers (everything after the first space),
    since tools like GeneMark use the full header as the sequence ID and
    spaces in GFF seqids break downstream parsing.

    Returns the path to the sanitized genome (or the original if already clean).
    """
    needs_cleaning = False
    with open(genome) as f:
        for line in f:
            if line.startswith(">") and " " in line.strip():
                needs_cleaning = True
                break

    if not needs_cleaning:
        return genome

    sanitized = work_dir / "genome.sanitized.fasta"
    if sanitized.exists():
        return sanitized

    log.info("Sanitizing genome FASTA headers (stripping descriptions)...")
    with open(genome) as fin, open(sanitized, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line.split()[0] + "\n")
            else:
                fout.write(line)

    return sanitized


def validate_gff(path: Path) -> None:
    """Verify that a file is valid GFF3.

    Uses streaming validation (no full DB load) for efficiency on large files.
    Raises GFFValidationError on invalid input.
    """
    from eukan.infra.logging import validate_gff as _validate_gff_streaming
    _validate_gff_streaming(path)
