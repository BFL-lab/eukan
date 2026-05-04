"""Terminal formatting and TSV export for annotation comparison results."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import astuple
from pathlib import Path
from statistics import median

from eukan.stats.models import (
    PAIR_TEST_TSV_COLUMNS,
    TSV_COLUMNS,
    ComparisonResult,
    GeneStats,
    MultiComparisonResult,
    SubfeatureStats,
)

# Frame width matches the existing single-pred output.
_W = 72
# Column header prepended to the multi-pred details TSV.
_PREDICTION_COL = "prediction"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def _bar(value: float, width: int = 20) -> str:
    filled = int(value * width)
    return f"[{'#' * filled}{'.' * (width - filled)}]"


def _significance_marker(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ---------------------------------------------------------------------------
# Single-prediction terminal blocks (shared with multi-prediction output)
# ---------------------------------------------------------------------------


def _title_block(ref_path: str, pred_path: str, w: int) -> list[str]:
    return [
        "=" * w,
        "ANNOTATION QUALITY ASSESSMENT",
        "=" * w,
        f"  Reference:  {ref_path}",
        f"  Predicted:  {pred_path}",
    ]


def _format_gene_block(gs: GeneStats, w: int) -> list[str]:
    """Detailed gene-level breakdown. Begins with a blank separator line."""
    lines: list[str] = [
        "",
        "-" * w,
        "  GENE Level",
        "-" * w,
        f"  Reference genes: {gs.ref_total:,}    Predicted genes: {gs.pred_total:,}",
        "",
    ]

    lines.append("  Reference gene classification:")
    total = gs.ref_total or 1
    for label, count in [
        ("Exact", gs.exact), ("Inexact", gs.inexact),
        ("Missing (FN)", gs.missing), ("Merged", gs.merged),
        ("Fragmented", gs.fragmented),
    ]:
        lines.append(
            f"    {label:<16} {count:>5,}  "
            f"({count / total * 100:5.1f}%)  "
            f"{_bar(count / total)}"
        )

    lines.append("")
    lines.append(f"  Novel predictions (FP):  {gs.novel:,}")
    lines.append("")

    lines.append("  Count-based metrics:")
    lines.append(f"    Sensitivity:  {_pct(gs.sensitivity):>6}  (TP={gs.tp}, FN={gs.fn})")
    lines.append(f"    Precision:    {_pct(gs.precision):>6}  (TP={gs.tp}, FP={gs.fp})")
    lines.append(f"    F1:           {_pct(gs.f1):>6}")
    lines.append("")

    if gs.sn_values:
        n = len(gs.sn_values)
        lines.append(f"  Overlap-based metrics (n={n:,} matched genes):")
        lines.append(f"    Mean Sn (ovl/ref):   {_pct(gs.mean_sn):>6}  {_bar(gs.mean_sn)}")
        lines.append(f"    Mean Sp (ovl/pred):  {_pct(gs.mean_sp):>6}  {_bar(gs.mean_sp)}")
        lines.append(f"    Mean F1:             {_pct(gs.mean_f1):>6}  {_bar(gs.mean_f1)}")
        lines.append(
            f"    Median Sn: {_pct(median(gs.sn_values)):>6}  |  "
            f"Median Sp: {_pct(median(gs.sp_values)):>6}"
        )
        perfect_sn = sum(1 for v in gs.sn_values if v >= 0.99)
        perfect_sp = sum(1 for v in gs.sp_values if v >= 0.99)
        lines.append(
            f"    Perfect Sn (>=99%): {perfect_sn:,}/{n:,}  |  "
            f"Perfect Sp (>=99%): {perfect_sp:,}/{n:,}"
        )

    if gs.boundary_5p:
        abs_5p = [abs(v) for v in gs.boundary_5p]
        abs_3p = [abs(v) for v in gs.boundary_3p]
        lines.append("")
        lines.append(f"  Boundary differences (inexact matches, n={len(gs.boundary_5p):,}):")
        lines.append(
            f"    5' mean |diff|: {sum(abs_5p) / len(abs_5p):.0f} bp  "
            f"median: {median(abs_5p):.0f} bp"
        )
        lines.append(
            f"    3' mean |diff|: {sum(abs_3p) / len(abs_3p):.0f} bp  "
            f"median: {median(abs_3p):.0f} bp"
        )

    return lines


def _format_subfeat_block(ss: SubfeatureStats, w: int) -> list[str]:
    """Subfeature-level (mRNA/CDS/intron) breakdown. Begins with a blank line."""
    lines: list[str] = [
        "",
        "-" * w,
        f"  {ss.level_name} Level (within matched parent pairs)",
        "-" * w,
        f"  Reference: {ss.ref_total:,}    Predicted: {ss.pred_total:,}",
        "",
    ]

    total_ref = ss.ref_total or 1
    total_pred = ss.pred_total or 1
    lines.append("  Classification:")
    lines.append(
        f"    Match:        {ss.match:>5,}  "
        f"({ss.match / total_ref * 100:5.1f}% of ref)  "
        f"{_bar(ss.match / total_ref)}"
    )
    lines.append(
        f"    Missing (FN): {ss.missing:>5,}  "
        f"({ss.missing / total_ref * 100:5.1f}% of ref)"
    )
    lines.append(
        f"    FP:           {ss.fp:>5,}  "
        f"({ss.fp / total_pred * 100:5.1f}% of pred)"
    )
    lines.append("")

    lines.append("  Count-based metrics:")
    lines.append(f"    Sensitivity:  {_pct(ss.sensitivity):>6}")
    lines.append(f"    Precision:    {_pct(ss.precision):>6}")
    lines.append(f"    F1:           {_pct(ss.f1):>6}")

    if ss.sn_values:
        n = len(ss.sn_values)
        lines.append("")
        lines.append(f"  Overlap-based metrics (n={n:,} matched pairs):")
        lines.append(f"    Mean Sn (ovl/ref):   {_pct(ss.mean_sn):>6}  {_bar(ss.mean_sn)}")
        lines.append(f"    Mean Sp (ovl/pred):  {_pct(ss.mean_sp):>6}  {_bar(ss.mean_sp)}")
        lines.append(f"    Mean F1:             {_pct(ss.mean_f1):>6}  {_bar(ss.mean_f1)}")
        lines.append(
            f"    Median Sn: {_pct(median(ss.sn_values)):>6}  |  "
            f"Median Sp: {_pct(median(ss.sp_values)):>6}"
        )
        perfect_sn = sum(1 for v in ss.sn_values if v >= 0.99)
        perfect_sp = sum(1 for v in ss.sp_values if v >= 0.99)
        lines.append(
            f"    Perfect Sn (>=99%): {perfect_sn:,}/{n:,}  |  "
            f"Perfect Sp (>=99%): {perfect_sp:,}/{n:,}"
        )

    return lines


def _format_summary_table(result: ComparisonResult, w: int) -> list[str]:
    """Per-prediction summary table. Begins with a blank line."""
    gs = result.gene_stats
    lines: list[str] = ["", "=" * w, "  SUMMARY", "=" * w]
    header = (
        f"  {'Level':<8} {'Sn':>7} {'Prec':>7} {'F1':>7}"
        f"  |  {'OvlSn':>7} {'OvlSp':>7} {'OvlF1':>7}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))

    lines.append(
        f"  {'Gene':<8} "
        f"{_pct(gs.sensitivity):>7} "
        f"{_pct(gs.precision):>7} "
        f"{_pct(gs.f1):>7}"
        f"  |  "
        f"{_pct(gs.mean_sn):>7} "
        f"{_pct(gs.mean_sp):>7} "
        f"{_pct(gs.mean_f1):>7}"
    )
    for ss in [result.mrna_stats, result.cds_stats, result.intron_stats]:
        lines.append(
            f"  {ss.level_name:<8} "
            f"{_pct(ss.sensitivity):>7} "
            f"{_pct(ss.precision):>7} "
            f"{_pct(ss.f1):>7}"
            f"  |  "
            f"{_pct(ss.mean_sn):>7} "
            f"{_pct(ss.mean_sp):>7} "
            f"{_pct(ss.mean_f1):>7}"
        )
    lines.append("=" * w)
    return lines


def format_results(result: ComparisonResult) -> str:
    """Format a single-prediction comparison as a fixed-width terminal report."""
    lines: list[str] = []
    lines.extend(_title_block(result.ref_path, result.pred_path, _W))
    lines.extend(_format_gene_block(result.gene_stats, _W))
    for ss in [result.mrna_stats, result.cds_stats, result.intron_stats]:
        lines.extend(_format_subfeat_block(ss, _W))
    lines.extend(_format_summary_table(result, _W))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Multi-prediction terminal report
# ---------------------------------------------------------------------------


_LEVEL_DISPLAY: dict[str, str] = {
    "gene": "Gene", "mrna": "mRNA", "cds": "CDS", "intron": "Intron",
}


def _multi_title_block(result: MultiComparisonResult, w: int) -> list[str]:
    lines: list[str] = [
        "=" * w,
        "ANNOTATION QUALITY ASSESSMENT (multi-prediction)",
        "=" * w,
        f"  Reference:  {result.ref_path}",
        "  Predictions:",
    ]
    label_w = max((len(p.label) for p in result.per_prediction), default=1)
    for p in result.per_prediction:
        lines.append(f"    {p.label:<{label_w}}  {p.pred_path}")
    return lines


def _per_prediction_section(result: ComparisonResult, w: int) -> list[str]:
    """Per-prediction breakdown for one prediction inside a multi-pred report."""
    lines: list[str] = ["", "", f"### PREDICTION: {result.label} ###"]
    lines.extend(_format_gene_block(result.gene_stats, w))
    for ss in [result.mrna_stats, result.cds_stats, result.intron_stats]:
        lines.extend(_format_subfeat_block(ss, w))
    lines.extend(_format_summary_table(result, w))
    return lines


def _stats_for_level(result: ComparisonResult, level: str) -> GeneStats | SubfeatureStats:
    if level == "gene":
        return result.gene_stats
    if level == "mrna":
        return result.mrna_stats
    if level == "cds":
        return result.cds_stats
    if level == "intron":
        return result.intron_stats
    raise ValueError(f"unknown level: {level!r}")


def _comparative_f1_table(result: MultiComparisonResult) -> list[str]:
    """Headline F1 score by (level x prediction)."""
    if not result.per_prediction:
        return []
    labels = [p.label for p in result.per_prediction]
    cell_w = max(8, *(len(la) + 2 for la in labels))
    lines = ["", "  F1 by level (count-based) / prediction:"]
    header = "  " + " " * 8 + "".join(f"{la:>{cell_w}}" for la in labels)
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for level in ("gene", "mrna", "cds", "intron"):
        row = f"  {_LEVEL_DISPLAY[level]:<8}"
        for pred in result.per_prediction:
            stats = _stats_for_level(pred, level)
            row += f"{_pct(stats.f1):>{cell_w}}"
        lines.append(row)
    return lines


def _kappa_matrix_table(result: MultiComparisonResult) -> list[str]:
    """Cohen's kappa matrix at gene level."""
    if not result.kappa_matrix:
        return []
    labels = [p.label for p in result.per_prediction]
    cell_w = max(8, *(len(la) + 2 for la in labels))
    lines = ["", "  Cohen's kappa (gene-level classification, off-diagonal):"]
    header = "  " + " " * cell_w + "".join(f"{la:>{cell_w}}" for la in labels)
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for la in labels:
        row = f"  {la:<{cell_w}}"
        for lb in labels:
            if la == lb:
                row += f"{'—':>{cell_w}}"
            else:
                lo, hi = sorted([la, lb])
                kap = result.kappa_matrix.get((lo, hi))
                cell = "—" if kap is None or kap != kap else f"{kap: .3f}"
                row += f"{cell:>{cell_w}}"
        lines.append(row)
    return lines


def _powerset_table(result: MultiComparisonResult) -> list[str]:
    """Powerset of prediction subsets that matched each ref gene.

    Full enumeration up to N=6 predictions; condensed (all/none/uniques)
    above. Subsets are grouped and sorted by descending cardinality.
    """
    if not result.powerset_matched:
        return []
    n = len(result.per_prediction)
    total = sum(result.powerset_matched.values()) or 1
    lines = ["", "  Powerset of gene-level matches (exact|inexact):"]
    if n <= 6:
        by_size: dict[int, list[tuple[tuple[str, ...], int]]] = defaultdict(list)
        for subset, count in result.powerset_matched.items():
            by_size[len(subset)].append((subset, count))
        for size in sorted(by_size.keys(), reverse=True):
            for subset, count in sorted(by_size[size], key=lambda x: x[0]):
                label = "(none)" if not subset else "{" + ", ".join(subset) + "}"
                lines.append(
                    f"    {label:<40} {count:>6,}  ({count / total * 100:5.1f}%)"
                )
    else:
        labels_set = {p.label for p in result.per_prediction}
        all_set = tuple(sorted(labels_set))
        all_count = result.powerset_matched.get(all_set, 0)
        none_count = result.powerset_matched.get((), 0)
        lines.append(
            f"    Matched by all {n}:".ljust(44)
            + f" {all_count:>6,}  ({all_count / total * 100:5.1f}%)"
        )
        lines.append(
            "    Matched by none:".ljust(44)
            + f" {none_count:>6,}  ({none_count / total * 100:5.1f}%)"
        )
        lines.append("    Uniquely matched:")
        for label in sorted(labels_set):
            uniq = result.powerset_matched.get((label,), 0)
            lines.append(
                f"      Only by {label}".ljust(44)
                + f" {uniq:>6,}  ({uniq / total * 100:5.1f}%)"
            )
        accounted = all_count + none_count + sum(
            result.powerset_matched.get((la,), 0) for la in labels_set
        )
        other = sum(result.powerset_matched.values()) - accounted
        lines.append(
            "    Other combinations:".ljust(44)
            + f" {other:>6,}  ({other / total * 100:5.1f}%)"
        )
    lines.append(f"    Total ref genes: {sum(result.powerset_matched.values()):,}")
    return lines


def _significance_table(result: MultiComparisonResult) -> list[str]:
    """Pairwise significance tests; show only BH-adjusted p < 0.05."""
    if not result.pair_tests:
        return []
    sig = [t for t in result.pair_tests if t.p_adj < 0.05]
    lines = [
        "",
        "  Significance tests (BH-adjusted; * <0.05, ** <0.01, *** <0.001):",
    ]
    if not sig:
        lines.append("    No pairwise tests reached significance after BH adjustment.")
        return lines
    header = (
        f"    {'Level':<8} {'Test':<22} {'Pred A':<14} {'Pred B':<14} "
        f"{'Stat':>10} {'p_adj':>10}"
    )
    lines.append(header)
    lines.append("    " + "-" * (len(header) - 4))
    for t in sorted(sig, key=lambda x: (x.level, x.test, x.pred_a, x.pred_b)):
        marker = _significance_marker(t.p_adj)
        lines.append(
            f"    {t.level:<8} {t.test:<22} {t.pred_a:<14} {t.pred_b:<14} "
            f"{t.statistic:>10.4f} {t.p_adj:>10.4g} {marker}"
        )
    return lines


def _comparative_section(result: MultiComparisonResult, w: int) -> list[str]:
    lines: list[str] = ["", "=" * w, "  COMPARATIVE SUMMARY", "=" * w]
    lines.extend(_comparative_f1_table(result))
    lines.extend(_kappa_matrix_table(result))
    lines.extend(_powerset_table(result))
    lines.extend(_significance_table(result))
    lines.append("")
    lines.append("=" * w)
    return lines


def format_multi_results(result: MultiComparisonResult) -> str:
    """Format a multi-prediction comparison as a terminal report."""
    lines: list[str] = []
    lines.extend(_multi_title_block(result, _W))
    for pred in result.per_prediction:
        lines.extend(_per_prediction_section(pred, _W))
    lines.extend(_comparative_section(result, _W))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# TSV detail export
# ---------------------------------------------------------------------------


def _format_field(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def _write_single_details_tsv(result: ComparisonResult, output_file: Path) -> Path:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for record in result.records:
            fh.write("\t".join(_format_field(v) for v in astuple(record)) + "\n")
    return output_file


def _write_multi_details_tsv(result: MultiComparisonResult, output_file: Path) -> Path:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as fh:
        fh.write(_PREDICTION_COL + "\t" + "\t".join(TSV_COLUMNS) + "\n")
        for pred in result.per_prediction:
            for record in pred.records:
                fh.write(
                    pred.label
                    + "\t"
                    + "\t".join(_format_field(v) for v in astuple(record))
                    + "\n"
                )
    return output_file


def write_details_tsv(
    result: ComparisonResult | MultiComparisonResult,
    output_file: Path,
) -> Path:
    """Write per-feature detail records as a TSV file.

    For ``MultiComparisonResult``, prepends a ``prediction`` column carrying
    the per-row prediction label. Creates parent directories as needed.
    """
    if isinstance(result, MultiComparisonResult):
        return _write_multi_details_tsv(result, output_file)
    return _write_single_details_tsv(result, output_file)


def write_stats_tsv(result: MultiComparisonResult, output_file: Path) -> Path:
    """Write pairwise statistical-test results as a long-form TSV."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as fh:
        fh.write("\t".join(PAIR_TEST_TSV_COLUMNS) + "\n")
        for t in result.pair_tests:
            row = [
                t.pred_a,
                t.pred_b,
                t.level,
                t.test,
                f"{t.statistic:.6f}",
                f"{t.p_value:.6g}",
                f"{t.p_adj:.6g}",
                str(t.n_a),
                str(t.n_b),
            ]
            fh.write("\t".join(row) + "\n")
    return output_file
