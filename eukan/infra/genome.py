"""Shared utilities for accessing genome FASTA files.

Provides a small wrapper around :func:`Bio.SeqIO.index` that caches the
most recently accessed contig.  ``SeqIO.index`` uses lazy disk offsets
and re-reads each record on access, which is bad when callers walk many
features per contig (the common case here).  Caching one contig at a
time keeps RAM bounded but makes sequential same-contig access cheap.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class ContigIndex:
    """Lazy contig lookup with single-record caching.

    Behaves like ``dict[str, SeqRecord]`` but only one contig is held in
    memory at a time.  Designed for the common pattern of iterating
    features sorted by chromosome -- the cache hit rate is then ~100%
    and the working set is one chromosome instead of the whole genome.
    """

    __slots__ = ("_cache_id", "_cache_record", "_index")

    def __init__(self, fasta: str | Path) -> None:
        self._index = SeqIO.index(str(fasta), "fasta")
        self._cache_id: str | None = None
        self._cache_record: SeqRecord | None = None

    def __getitem__(self, contig_id: str) -> SeqRecord:
        if contig_id != self._cache_id:
            self._cache_id = contig_id
            self._cache_record = self._index[contig_id]
        # Either we just populated the cache, or contig_id matched
        assert self._cache_record is not None
        return self._cache_record

    def __contains__(self, contig_id: str) -> bool:
        return contig_id in self._index

    def get(self, contig_id: str, default: SeqRecord | None = None) -> SeqRecord | None:
        try:
            return self[contig_id]
        except KeyError:
            return default

    def __iter__(self):
        return iter(self._index)

    def close(self) -> None:
        self._index.close()
        self._cache_id = None
        self._cache_record = None

    def __enter__(self) -> ContigIndex:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()
