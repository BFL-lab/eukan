"""Data structures for annotation quality assessment."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import NamedTuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MERGE_THRESHOLD = 0.50  # prediction must cover >50% of ref to count as merged
FRAG_THRESHOLD = 0.50   # prediction must overlap >=50% of its own length


# ---------------------------------------------------------------------------
# Interval
# ---------------------------------------------------------------------------


class Interval(NamedTuple):
    chrom: str
    strand: str
    start: int
    end: int
    feat_id: str
    parent_id: str | None


# ---------------------------------------------------------------------------
# Per-feature record
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class FeatureRecord:
    """One row in the per-feature detail output.

    Fields after ``classification`` default to ``None`` where they don't
    apply (e.g. boundary diffs only for gene-level inexact matches,
    pred fields empty for missing refs, ref fields empty for novel preds).
    """

    level: str              # "gene" | "mRNA" | "CDS" | "intron"
    classification: str     # gene: exact|inexact|missing|merged|fragmented|novel
                            # subfeature: match|missing|fp
    ref_id: str | None = None
    pred_id: str | None = None
    chrom: str | None = None
    strand: str | None = None
    ref_start: int | None = None
    ref_end: int | None = None
    pred_start: int | None = None
    pred_end: int | None = None
    overlap_bp: int | None = None
    sn: float | None = None
    sp: float | None = None
    f1: float | None = None
    boundary_5p: int | None = None
    boundary_3p: int | None = None

    # -- factory classmethods --------------------------------------------------

    @classmethod
    def from_ref(
        cls, level: str, classification: str, iv: Interval, **kwargs,
    ) -> FeatureRecord:
        """Create a record for a reference feature (missing/merged/fragmented)."""
        return cls(
            level=level, classification=classification,
            ref_id=iv.feat_id, chrom=iv.chrom, strand=iv.strand,
            ref_start=iv.start, ref_end=iv.end, **kwargs,
        )

    @classmethod
    def from_pred(
        cls, level: str, classification: str, iv: Interval, **kwargs,
    ) -> FeatureRecord:
        """Create a record for a prediction feature (novel/fp)."""
        return cls(
            level=level, classification=classification,
            pred_id=iv.feat_id, chrom=iv.chrom, strand=iv.strand,
            pred_start=iv.start, pred_end=iv.end, **kwargs,
        )

    @classmethod
    def from_match(
        cls,
        level: str,
        classification: str,
        ref: Interval,
        pred: Interval,
        overlap_bp: int,
        **kwargs,
    ) -> FeatureRecord:
        """Create a record for a matched ref/pred pair."""
        ref_len = ref.end - ref.start + 1
        pred_len = pred.end - pred.start + 1
        sn = overlap_bp / ref_len if ref_len else 0.0
        sp = overlap_bp / pred_len if pred_len else 0.0
        f1 = 2 * sn * sp / (sn + sp) if (sn + sp) else 0.0
        return cls(
            level=level, classification=classification,
            ref_id=ref.feat_id, pred_id=pred.feat_id,
            chrom=ref.chrom, strand=ref.strand,
            ref_start=ref.start, ref_end=ref.end,
            pred_start=pred.start, pred_end=pred.end,
            overlap_bp=overlap_bp, sn=sn, sp=sp, f1=f1,
            **kwargs,
        )


# Column names for TSV export, derived from the dataclass fields.
TSV_COLUMNS: tuple[str, ...] = tuple(f.name for f in fields(FeatureRecord))


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------


def _safe_mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


@dataclass
class GeneStats:
    ref_total: int = 0
    pred_total: int = 0
    exact: int = 0
    inexact: int = 0
    missing: int = 0
    merged: int = 0
    fragmented: int = 0
    novel: int = 0  # predictions with no reference overlap
    # Overlap metrics for matched (exact + inexact) genes
    sn_values: list[float] = field(default_factory=list)
    sp_values: list[float] = field(default_factory=list)
    f1_values: list[float] = field(default_factory=list)
    boundary_5p: list[int] = field(default_factory=list)  # 5' boundary diffs
    boundary_3p: list[int] = field(default_factory=list)  # 3' boundary diffs

    @property
    def tp(self) -> int:
        return self.exact + self.inexact

    @property
    def fn(self) -> int:
        return self.missing + self.fragmented + self.merged

    @property
    def fp(self) -> int:
        return self.novel

    @property
    def sensitivity(self) -> float:
        d = self.tp + self.fn
        return self.tp / d if d else 0.0

    @property
    def precision(self) -> float:
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def f1(self) -> float:
        s, p = self.sensitivity, self.precision
        return 2 * s * p / (s + p) if (s + p) else 0.0

    @property
    def mean_sn(self) -> float:
        return _safe_mean(self.sn_values)

    @property
    def mean_sp(self) -> float:
        return _safe_mean(self.sp_values)

    @property
    def mean_f1(self) -> float:
        return _safe_mean(self.f1_values)


@dataclass
class SubfeatureStats:
    """Statistics for mRNA, CDS, or intron level."""

    level_name: str
    ref_total: int = 0
    pred_total: int = 0
    match: int = 0
    missing: int = 0
    fp: int = 0
    sn_values: list[float] = field(default_factory=list)
    sp_values: list[float] = field(default_factory=list)
    f1_values: list[float] = field(default_factory=list)

    @property
    def sensitivity(self) -> float:
        d = self.match + self.missing
        return self.match / d if d else 0.0

    @property
    def precision(self) -> float:
        d = self.match + self.fp
        return self.match / d if d else 0.0

    @property
    def f1(self) -> float:
        s, p = self.sensitivity, self.precision
        return 2 * s * p / (s + p) if (s + p) else 0.0

    @property
    def mean_sn(self) -> float:
        return _safe_mean(self.sn_values)

    @property
    def mean_sp(self) -> float:
        return _safe_mean(self.sp_values)

    @property
    def mean_f1(self) -> float:
        return _safe_mean(self.f1_values)


# ---------------------------------------------------------------------------
# Comparison result
# ---------------------------------------------------------------------------


@dataclass
class ComparisonResult:
    gene_stats: GeneStats
    mrna_stats: SubfeatureStats
    cds_stats: SubfeatureStats
    intron_stats: SubfeatureStats
    ref_path: str
    pred_path: str
    records: list[FeatureRecord] = field(default_factory=list)
    # Short identifier used in multi-prediction output. Defaults to the
    # prediction file's stem in ``compare_annotations``; the multi-driver
    # may override with a user-supplied label.
    label: str = ""


# ---------------------------------------------------------------------------
# Multi-prediction result
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class PairTest:
    """One pairwise statistical test between two predictions at one feature level.

    For ``test="ks_*"``, ``n_a``/``n_b`` are the matched-feature sample sizes
    feeding the ECDFs. For ``test="chi2_classification"`` they are the
    per-prediction row totals of the contingency table.
    """

    pred_a: str
    pred_b: str
    level: str        # "gene" | "mRNA" | "CDS" | "intron"
    test: str         # "ks_sn" | "ks_sp" | "ks_f1" | "chi2_classification"
    statistic: float
    p_value: float
    p_adj: float
    n_a: int
    n_b: int


@dataclass
class MultiComparisonResult:
    """Per-prediction comparison results plus inter-prediction summaries."""

    ref_path: str
    per_prediction: list[ComparisonResult]
    # All pairwise tests, BH-adjusted within (level, test) families.
    pair_tests: list[PairTest] = field(default_factory=list)
    # Off-diagonal Cohen's kappa at gene level. Key is ``(a, b)`` with
    # ``a < b`` lexicographically.
    kappa_matrix: dict[tuple[str, str], float] = field(default_factory=dict)
    # Among ref genes (gene-level), the count matched (exact|inexact) by
    # exactly that subset of predictions. Subsets are sorted tuples of
    # labels; the empty tuple counts genes no prediction matched.
    powerset_matched: dict[tuple[str, ...], int] = field(default_factory=dict)


# Long-form TSV columns for the pairwise stats output.
PAIR_TEST_TSV_COLUMNS: tuple[str, ...] = tuple(f.name for f in fields(PairTest))


# ---------------------------------------------------------------------------
# Derive summary stats from per-feature records
# ---------------------------------------------------------------------------


def gene_stats_from_records(
    records: list[FeatureRecord],
    ref_total: int,
    pred_total: int,
) -> GeneStats:
    """Build a GeneStats from gene-level FeatureRecords."""
    stats = GeneStats(ref_total=ref_total, pred_total=pred_total)
    for r in records:
        cls = r.classification
        if cls == "exact":
            stats.exact += 1
        elif cls == "inexact":
            stats.inexact += 1
            if r.boundary_5p is not None:
                stats.boundary_5p.append(r.boundary_5p)
            if r.boundary_3p is not None:
                stats.boundary_3p.append(r.boundary_3p)
        elif cls == "missing":
            stats.missing += 1
        elif cls == "merged":
            stats.merged += 1
        elif cls == "fragmented":
            stats.fragmented += 1
        elif cls == "novel":
            stats.novel += 1

        # Overlap metrics for matched genes (exact + inexact).
        # sn/sp/f1 are set as a triple by FeatureRecord.from_match, so a
        # single None check on sn implies the others are also non-None.
        if cls in ("exact", "inexact") and r.sn is not None:
            stats.sn_values.append(r.sn)
            assert r.sp is not None and r.f1 is not None
            stats.sp_values.append(r.sp)
            stats.f1_values.append(r.f1)
    return stats


def subfeature_stats_from_records(
    records: list[FeatureRecord],
    level_name: str,
    ref_total: int,
    pred_total: int,
) -> SubfeatureStats:
    """Build a SubfeatureStats from level-specific FeatureRecords."""
    stats = SubfeatureStats(
        level_name=level_name, ref_total=ref_total, pred_total=pred_total,
    )
    for r in records:
        cls = r.classification
        if cls == "match":
            stats.match += 1
            if r.sn is not None:
                assert r.sp is not None and r.f1 is not None
                stats.sn_values.append(r.sn)
                stats.sp_values.append(r.sp)
                stats.f1_values.append(r.f1)
        elif cls == "missing":
            stats.missing += 1
        elif cls == "fp":
            stats.fp += 1
    return stats
