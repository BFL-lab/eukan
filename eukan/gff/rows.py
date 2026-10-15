"""Line-level GFF3/GTF reading: data rows and column-9 attributes.

The stdlib-only fast path for code that needs only a few columns, where a gffutils
database would be overkill.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path


def iter_gff3_rows(path: str | Path) -> Iterator[list[str]]:
    """Yield the tab-split columns of each data row in a GFF3/GTF file.

    Skips ``#`` lines (comments and directives), blank lines, and rows with
    fewer than 9 columns. Only the newline is stripped, so the last column keeps
    any trailing ``\\r`` or spaces.
    """
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 9:
                yield cols


def parse_attrs(col9: str) -> dict[str, str]:
    """Parse column 9 of either GFF3 (``key=value``) or GTF (``key "value"``)."""
    attrs: dict[str, str] = {}
    for part in col9.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:  # GFF3
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
        else:  # GTF: key "value"
            bits = part.split(None, 1)
            if len(bits) == 2:
                attrs[bits[0].strip()] = bits[1].strip().strip('"')
    return attrs
