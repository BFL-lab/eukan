"""Annotation quality assessment — re-exports from eukan.compare.

This module is kept for backwards compatibility with tests/run_pipeline.py.
The implementation now lives in eukan/compare/.
"""

from eukan.compare import (  # noqa: F401
    compare_annotations,
    format_results,
    write_details_tsv,
)
from eukan.compare.models import (  # noqa: F401
    ComparisonResult,
    FeatureRecord,
    GeneStats,
    Interval,
    SubfeatureStats,
)
