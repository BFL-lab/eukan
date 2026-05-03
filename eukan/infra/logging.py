"""Logging configuration and shared utilities."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

CHUNK_SIZE = 8 * 1024 * 1024  # 8MB


def md5_file(path: Path) -> str:
    """Compute MD5 of a file in chunks (never loads entire file into memory)."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def get_logger(name: str) -> logging.Logger:
    """Get a logger for the given module name.

    All eukan loggers are children of the 'eukan' root logger,
    which is configured once by setup_logging().
    """
    return logging.getLogger(name)


def setup_logging(verbosity: int = 0) -> None:
    """Configure the eukan logging hierarchy.

    Args:
        verbosity: 0=INFO, 1=DEBUG, -1=WARNING (quiet).
    """
    level = {-1: logging.WARNING, 0: logging.INFO, 1: logging.DEBUG}.get(
        verbosity, logging.INFO
    )

    root = logging.getLogger("eukan")
    if root.handlers:
        return  # already configured

    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%H:%M:%S")
    )
    root.addHandler(handler)
    root.setLevel(level)


def count_gff3_features(gff3_path: Path, feature_type: str = "gene") -> int:
    """Count features of a given type in a GFF3 file by scanning column 3.

    Fast line-based parsing — does not load the file into a database.
    """
    count = 0
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 3 and cols[2] == feature_type:
                count += 1
    return count


def validate_gff(path: Path) -> bool:
    """Check if a file is valid GFF3 by streaming a few features.

    Uses gffutils.DataIterator to parse without building a full SQLite
    database, making this fast even for large files.

    Raises:
        GFFValidationError: If the file cannot be parsed or has invalid features.
    """
    import gffutils

    from eukan.exceptions import GFFValidationError

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
