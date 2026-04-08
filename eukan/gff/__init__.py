"""GFF3 format operations."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import gffutils

# Standard GFF3 dialect for gffutils — single source of truth.
# Import as: from eukan.gff import GFF3_DIALECT
GFF3_DIALECT: dict[str, object] = {
    "leading semicolon": False,
    "trailing semicolon": True,
    "quoted GFF2 values": False,
    "field separator": ";",
    "keyval separator": "=",
    "multival separator": ",",
    "fmt": "gff3",
    "repeated keys": False,
    "order": ["ID", "Name", "Parent"],
}


def create_gff_db(
    source: str | Path | gffutils.FeatureDB,
    *,
    transform: Callable | None = None,
    dialect: dict | None = None,
    merge_strategy: str = "create_unique",
    verbose: bool = False,
    **kwargs,
) -> gffutils.FeatureDB:
    """Create a gffutils in-memory database with standardised defaults.

    Wraps ``gffutils.create_db()`` with consistent defaults for dialect,
    verbose, and merge_strategy so callers don't need to repeat them.
    """
    data = str(source) if isinstance(source, Path) else source
    return gffutils.create_db(
        data,
        ":memory:",
        dialect=dialect or GFF3_DIALECT,
        transform=transform,
        merge_strategy=merge_strategy,
        verbose=verbose,
        **kwargs,
    )


def transform_db(
    db: gffutils.FeatureDB,
    transform: Callable,
    *,
    dialect: dict | None = None,
    merge_strategy: str = "create_unique",
    verbose: bool = False,
    **kwargs,
) -> gffutils.FeatureDB:
    """Re-create a database by applying a transform to an existing one.

    This is the common pattern of ``gffutils.create_db(existing_db, ...)``
    with a new transform function.
    """
    return gffutils.create_db(
        db,
        ":memory:",
        dialect=dialect or GFF3_DIALECT,
        transform=transform,
        merge_strategy=merge_strategy,
        verbose=verbose,
        **kwargs,
    )
