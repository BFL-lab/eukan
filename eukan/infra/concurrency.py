"""Shared concurrency helper: fan out, collect in input order, raise on the
first worker failure.

Used by the annotation pipeline to run independent steps side by side (e.g.
ORF finding with GeneMark) and by AUGUSTUS split prediction.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    fn: Callable[[T], R],
    items: Iterable[T],
    *,
    max_workers: int | None = None,
) -> list[R]:
    """Run *fn* over *items* concurrently and return results in input order.

    Equivalent to ``list(map(fn, items))`` but executed across a
    ``ThreadPoolExecutor``. Raises the first worker exception encountered.

    *max_workers* defaults to the number of items: one thread per item.
    """
    items_list = list(items)
    if not items_list:
        return []
    workers = max_workers if max_workers is not None else len(items_list)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fn, items_list))
