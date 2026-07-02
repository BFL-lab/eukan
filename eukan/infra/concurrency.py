"""Shared concurrency helper for the four ThreadPool sites in eukan.

The ``ThreadPoolExecutor`` pattern was reproduced ad-hoc in the
annotation pipeline, AUGUSTUS split-prediction, consensus partition
execution, and the concordance pass triad. Each site has a slightly
different shape but the same intent: fan out, collect, raise on first
worker failure.
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

    *max_workers* defaults to the number of items, capped by Python's
    default sensible limit when ``None`` is passed.
    """
    items_list = list(items)
    if not items_list:
        return []
    workers = max_workers if max_workers is not None else len(items_list)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fn, items_list))
