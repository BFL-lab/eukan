"""Terminal formatting and TSV export for annotation comparison results."""

from __future__ import annotations

from dataclasses import astuple
from pathlib import Path
from statistics import median

from eukan.stats.models import ComparisonResult, FeatureRecord, TSV_COLUMNS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def _bar(value: float, width: int = 20) -> str:
    filled = int(value * width)
    return f"[{'#' * filled}{'.' * (width - filled)}]"


# ---------------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------------


def format_results(result: ComparisonResult) -> str:
    lines: list[str] = []
    w = 72

    lines.append("=" * w)
    lines.append("ANNOTATION QUALITY ASSESSMENT")
    lines.append("=" * w)
    lines.append(f"  Reference:  {result.ref_path}")
    lines.append(f"  Predicted:  {result.pred_path}")

    # ---- Gene level ----
    gs = result.gene_stats
    lines.append("")
    lines.append("-" * w)
    lines.append("  GENE Level")
    lines.append("-" * w)
    lines.append(f"  Reference genes: {gs.ref_total:,}    Predicted genes: {gs.pred_total:,}")
    lines.append("")

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

    # ---- mRNA / CDS / Intron levels ----
    for ss in [result.mrna_stats, result.cds_stats, result.intron_stats]:
        lines.append("")
        lines.append("-" * w)
        lines.append(f"  {ss.level_name} Level (within matched parent pairs)")
        lines.append("-" * w)
        lines.append(f"  Reference: {ss.ref_total:,}    Predicted: {ss.pred_total:,}")
        lines.append("")

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

    # ---- Summary table ----
    lines.append("")
    lines.append("=" * w)
    lines.append("  SUMMARY")
    lines.append("=" * w)
    header = (
        f"  {'Level':<8} {'Sn':>7} {'Prec':>7} {'F1':>7}"
        f"  |  {'OvlSn':>7} {'OvlSp':>7} {'OvlF1':>7}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))

    # Gene row
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


def write_details_tsv(result: ComparisonResult, output_file: Path) -> Path:
    """Write per-feature detail records as a TSV file.

    Creates parent directories if needed. Returns the path written.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for record in result.records:
            fh.write(
                "\t".join(_format_field(v) for v in astuple(record)) + "\n"
            )

    return output_file
