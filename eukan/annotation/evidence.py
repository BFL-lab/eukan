"""Shared evidence-weighting helpers for the consensus engine.

Maps each staged evidence file to its weights.txt class + GFF source token —
the single source of truth used by the combinr consensus engine
(:mod:`eukan.annotation.combinr_consensus`). Kept in its own module so combinr
can import it without pulling in any engine driver.
"""

from __future__ import annotations

from pathlib import Path

from eukan.gff.rows import iter_gff3_rows


def _first_source_token(gff3: Path) -> str | None:
    """Return the GFF3 source column (col 2) of the first data line, or None."""
    for cols in iter_gff3_rows(gff3):
        if cols[1]:
            return cols[1]
    return None


# Evidence basename -> (weights class, GFF source token). Single source of
# truth for how a staged evidence file maps to its weights.txt entry. PROTEIN
# uses weights[0], ABINITIO_PREDICTION uses weights[1] (see callers).
EVIDENCE_ROLES: dict[str, tuple[str, str]] = {
    "prot.gff3":         ("PROTEIN",             "prot_align"),
    "augustus.gff3":     ("ABINITIO_PREDICTION", "augustus"),
    "snap.gff3":         ("ABINITIO_PREDICTION", "snap"),
    "genemark.gff3":     ("ABINITIO_PREDICTION", "genemark"),
    "codingquarry.gff3": ("ABINITIO_PREDICTION", "codingquarry"),
}
