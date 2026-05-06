"""Workarounds for gffutils API quirks.

A single home for the small "the library should support this directly,
but doesn't" routines so the reasoning lives in one place.
"""

from __future__ import annotations

import gffutils


def empty_db() -> gffutils.FeatureDB:
    """Create an empty in-memory FeatureDB.

    Workaround for gffutils raising ``EmptyInputError`` on empty input.
    Creates a DB with a placeholder feature, then deletes it. The
    resulting DB has the standard schema and supports every operation
    callers expect, just with no features.
    """
    db = gffutils.create_db(
        "chr1\t.\tgene\t1\t1\t.\t+\t.\tID=_empty_placeholder",
        ":memory:", from_string=True,
    )
    db.delete(db["_empty_placeholder"])
    return db
