"""Input validation and sanitization for FASTA and GFF3 files.

Cross-cutting boundary checks used by every pipeline at the input edge.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from eukan.exceptions import FastaValidationError, GFFValidationError
from eukan.infra.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# GFF3
# ---------------------------------------------------------------------------


def validate_gff(path: Path) -> bool:
    """Check if a file is valid GFF3 by streaming a few features.

    Uses gffutils.DataIterator to parse without building a full SQLite
    database, making this fast even for large files.

    Raises:
        GFFValidationError: If the file cannot be parsed or has invalid features.
    """
    import gffutils

    try:
        count = 0
        for f in gffutils.DataIterator(str(path)):
            if f.start is None or f.end is None or len(f.attributes) == 0:
                raise GFFValidationError(path, "feature has missing attributes or coordinates")
            count += 1
            if count >= 10:
                break
        if count == 0:
            raise GFFValidationError(path, "file contains no features")
        return True
    except GFFValidationError:
        raise
    except Exception as exc:
        raise GFFValidationError(path, str(exc)) from exc
