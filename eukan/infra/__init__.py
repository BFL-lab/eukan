"""Runtime infrastructure: subprocess execution, step lifecycle, manifest tracking, logging."""

from __future__ import annotations

import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib  # type: ignore[import-not-found,no-redef,unused-ignore]

__all__ = ["tomllib"]
