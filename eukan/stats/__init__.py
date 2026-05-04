"""Annotation quality assessment: compare predicted GFF3 against a reference."""

from eukan.stats.compare import compare_annotations, compare_multiple
from eukan.stats.format import (
    format_multi_results,
    format_results,
    write_details_tsv,
    write_stats_tsv,
)

__all__ = [
    "compare_annotations",
    "compare_multiple",
    "format_multi_results",
    "format_results",
    "write_details_tsv",
    "write_stats_tsv",
]
