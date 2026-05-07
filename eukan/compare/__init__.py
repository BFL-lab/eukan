"""Annotation quality assessment: compare predicted GFF3 against a reference."""

from eukan.compare.engine import compare_annotations, compare_multiple
from eukan.compare.format import (
    format_multi_results,
    format_results,
    write_details_tsv,
)

__all__ = [
    "compare_annotations",
    "compare_multiple",
    "format_multi_results",
    "format_results",
    "write_details_tsv",
]
