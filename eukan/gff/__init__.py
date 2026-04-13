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

    Wraps ``gffutils.create_db()`` with consistent defaults for verbose
    and merge_strategy so callers don't need to repeat them.  The dialect
    is auto-detected by default; pass ``dialect=GFF3_DIALECT`` explicitly
    when a specific dialect is needed.
    """
    data = str(source) if isinstance(source, Path) else source
    kw: dict = dict(
        merge_strategy=merge_strategy,
        verbose=verbose,
        **kwargs,
    )
    if transform is not None:
        kw["transform"] = transform
    if dialect is not None:
        kw["dialect"] = dialect
    return gffutils.create_db(data, ":memory:", **kw)


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

    Uses GFF3_DIALECT by default since the input is typically a
    previously parsed FeatureDB that should produce GFF3 output.
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
