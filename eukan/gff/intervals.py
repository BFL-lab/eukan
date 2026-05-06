"""Strand-aware interval index with prefix-max-end overlap queries.

Used by both the consensus/concordance routines (``gff/intersecter.py``)
and the annotation comparison stats (``stats/compare.py``). Both modules
previously rolled their own copy of this data structure.

Scheme:

* Group intervals by ``(chrom, strand)``.
* Within each bucket, sort by ``start``.
* Build a parallel monotone running-max-end array
  (``max_end_to[i] = max(end[0..i])``).

Then bound the candidate set for an overlap query with two binary
searches:

* ``bisect_right(starts, query.end)`` → first interval whose start is
  past the query (right bound).
* ``bisect_left(max_end_to, query.start)`` → first interval that *might*
  reach the query (left bound). Earlier intervals strictly end before
  the query and can be skipped.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Iterable, Iterator
from typing import Generic, Protocol, TypeVar


class _HasInterval(Protocol):
    """Minimal interval-like protocol: chrom/strand/start/end attributes."""

    @property
    def chrom(self) -> str: ...

    @property
    def strand(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...


IT = TypeVar("IT", bound=_HasInterval)


def build_max_end_prefix(ends: Iterable[int]) -> list[int]:
    """Build a monotone running-max-end array.

    Given interval ends in start-sorted order, returns
    ``out[i] = max(ends[:i+1])``.  ``bisect_left(out, query.start)`` then
    gives the leftmost index whose end could reach the query, bounding
    the candidate set instead of scanning every interval.
    """
    out: list[int] = []
    running_max = 0
    for e in ends:
        if e > running_max:
            running_max = e
        out.append(running_max)
    return out


class IntervalIndex(Generic[IT]):
    """A spatial index of intervals keyed by ``(chrom, strand)``.

    Generic over any object exposing ``chrom``, ``strand``, ``start``,
    and ``end`` attributes — e.g. ``gffutils.Feature`` and
    ``stats.models.Interval`` both qualify.

    Use :meth:`overlapping` to enumerate hits along with their overlap
    bp counts.
    """

    __slots__ = ("_buckets",)

    def __init__(self, intervals: Iterable[IT]) -> None:
        groups: dict[tuple[str, str], list[IT]] = defaultdict(list)
        for iv in intervals:
            groups[(iv.chrom, iv.strand)].append(iv)

        self._buckets: dict[tuple[str, str], tuple[list[IT], list[int], list[int]]] = {}
        for key, ivs in groups.items():
            ivs.sort(key=lambda x: x.start)
            starts = [iv.start for iv in ivs]
            max_end_to = build_max_end_prefix(iv.end for iv in ivs)
            self._buckets[key] = (ivs, starts, max_end_to)

    def overlapping(self, query: _HasInterval) -> Iterator[tuple[IT, int]]:
        """Yield ``(target, overlap_bp)`` for each strand-aware overlap.

        Overlap bp is the inclusive intersection length, i.e.
        ``min(end_a, end_b) - max(start_a, start_b) + 1`` — at least 1
        for any returned target.
        """
        entry = self._buckets.get((query.chrom, query.strand))
        if entry is None:
            return
        ivs, starts, max_end_to = entry
        idx_right = bisect_right(starts, query.end)
        if idx_right == 0:
            return
        idx_left = bisect_left(max_end_to, query.start)
        for i in range(idx_left, idx_right):
            t = ivs[i]
            if t.end < query.start:
                continue
            ovl = min(query.end, t.end) - max(query.start, t.start) + 1
            yield t, ovl

    def has_overlap(self, query: _HasInterval) -> bool:
        """``True`` iff at least one indexed interval overlaps *query*."""
        for _ in self.overlapping(query):
            return True
        return False

    def __contains__(self, key: tuple[str, str]) -> bool:
        return key in self._buckets
