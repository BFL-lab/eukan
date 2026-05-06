"""GFF3 structural validation."""

from __future__ import annotations

from pathlib import Path

from eukan.exceptions import GFFValidationError


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
