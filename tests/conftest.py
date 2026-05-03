"""Shared pytest fixtures.

Anything used in more than one test module belongs here so the fixtures
stay in lockstep when their underlying APIs change.
"""

from __future__ import annotations

import gffutils
import pytest


@pytest.fixture
def db_from_string():
    """Build an in-memory gffutils.FeatureDB from a GFF3 string.

    Returns a callable so a single test can construct multiple dbs::

        def test_thing(db_from_string):
            db = db_from_string("chr1\\t.\\tgene\\t1\\t100\\t.\\t+\\t.\\tID=g1\\n")
    """
    def _build(gff: str) -> gffutils.FeatureDB:
        return gffutils.create_db(
            gff, ":memory:", from_string=True, merge_strategy="create_unique",
        )
    return _build
