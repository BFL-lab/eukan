"""Runtime infrastructure: subprocess execution, step lifecycle, manifest tracking, logging."""

from __future__ import annotations

import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib  # type: ignore[no-redef]

__all__ = ["tomllib"]
