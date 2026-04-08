"""Annotation quality assessment — re-exports from eukan.stats.

This module is kept for backwards compatibility with tests/run_pipeline.py.
The implementation now lives in eukan/stats/.
"""

from eukan.stats import (  # noqa: F401
    compare_annotations,
    format_results,
    write_details_tsv,
)
from eukan.stats.models import (  # noqa: F401
    ComparisonResult,
    FeatureRecord,
    GeneStats,
    Interval,
    SubfeatureStats,
)
