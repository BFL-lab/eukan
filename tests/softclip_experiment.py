#!/usr/bin/env python3
"""Run the soft-clip / intron diagnostic on the hemistasia and flectonema BAMs.

Three views the diagnostic surfaces:

- **Soft-clip clusters** — anchored-substring keys with the count of
  distinct loci that support each. SL motifs show up here as keys with
  thousands of supporting loci.
- **Intron dinucleotides** — per-N-cigar-op donor/acceptor pair counts,
  to distinguish canonical from non-canonical splice landscapes.
- **Per-locus motif consolidation** — for each genomic locus that
  attracted soft-clipped reads, do all the clips consolidate to one
  motif? For consistent loci, how widely is that motif shared (SL-like
  → many other loci share it; locus-specific → only this locus)?

Usage::

    conda run -n eukan python tests/softclip_experiment.py run-all
    conda run -n eukan python tests/softclip_experiment.py one \\
        --bam tests/data/flectonema/STAR_Aligned.sortedByCoord.out.bam \\
        --genome tests/data/flectonema/GCA_964019425.1_peFleSpea1.1_genomic.fna \\
        --label flectonema

Outputs persisted next to each BAM:

- ``softclip_diagnostic.json``           — slim summary (a few KB)
- ``softclip_diagnostic_clusters.tsv``   — every cluster, sorted by n_loci
- ``softclip_diagnostic_loci.tsv``       — every non-singleton locus,
                                            sorted by n_clips
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import click

from eukan.assembly.bam_diagnostic import (
    DiagnosticParams,
    DiagnosticReport,
    IntronStats,
    LocusConsistencyStats,
    SoftClipStats,
    Verdict,
    compute_verdict,
    diagnose_bam,
    to_summary_dict,
)

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "tests" / "data"

KNOWN = {
    "flectonema": (
        DATA / "flectonema" / "STAR_Aligned.sortedByCoord.out.bam",
        DATA / "flectonema" / "GCA_964019425.1_peFleSpea1.1_genomic.fna",
    ),
    "hemistasia": (
        DATA / "hemistasia" / "STAR_Aligned.sortedByCoord.out.bam",
        DATA / "hemistasia" / "Hemistasia_assembly.fasta",
    ),
}

# Order the motif-share buckets are printed in, regardless of insertion
# order in the histogram dict.
_BUCKET_ORDER = ["1", "2-10", "11-100", "101-1000", ">1000"]


def _pct(n: int, d: int) -> str:
    return f"{(100.0 * n / d) if d else 0.0:.2f}%"


def _print_summary(label: str, soft: SoftClipStats, intron: IntronStats) -> None:
    bar = "=" * 70
    click.echo(f"\n{bar}\n {label}\n{bar}")
    click.echo(f"  reads scanned:           {soft.n_reads_scanned:>12,}")
    click.echo(
        f"  reads with soft-clip:    {soft.n_reads_with_clip:>12,} "
        f"({_pct(soft.n_reads_with_clip, soft.n_reads_scanned)})"
    )
    click.echo(f"  soft-clips total:        {soft.n_clips_total:>12,}")
    for side in ("5p", "3p"):
        n = soft.n_clips_by_side.get(side, 0)
        click.echo(f"    {side} (mRNA orient):     {n:>12,}")
    click.echo(f"  distinct loci:           {soft.n_loci:>12,}")
    click.echo(f"  distinct cluster keys:   {soft.n_clusters:>12,}")
    if soft.n_clusters:
        ratio = soft.n_loci / soft.n_clusters
        click.echo(f"  loci / cluster ratio:    {ratio:>12.2f}")

    click.echo("")
    click.echo(f"  introns (N cigar ops):   {intron.n_introns_total:>12,}")
    click.echo(f"  canonical (GT-AG+CT-AC): {intron.canonical_pct:>12.2f}%")
    if intron.by_dinucleotide:
        click.echo("  top dinucleotides:")
        for dinuc, n in list(intron.by_dinucleotide.items())[:8]:
            click.echo(f"    {dinuc:<6}  {n:>12,}  ({_pct(n, intron.n_introns_total)})")


def _print_verdict(verdict: Verdict) -> None:
    ts = verdict.trans_splicing
    nc = verdict.non_canonical_splice
    click.echo("\n  Empirical verdict:")
    click.echo(f"    Trans-splicing: {ts.call}")
    if ts.top_non_trivial_cluster_key:
        click.echo(
            f"      Top non-trivial cluster: {ts.top_non_trivial_cluster_key} "
            f"({ts.top_non_trivial_cluster_n_loci:,} loci, "
            f"{ts.top_non_trivial_cluster_n_reads:,} reads)"
        )
        if ts.top_non_trivial_cluster_consensus:
            click.echo(
                f"      Consensus motif:         {ts.top_non_trivial_cluster_consensus}"
            )
    else:
        click.echo("      Top non-trivial cluster: (none — all top clusters are low-complexity)")
    click.echo(
        f"      Motif-share bucket >1000: "
        f"{ts.sl_bucket_pct_of_consistent:.2f}% of consistent loci"
    )
    click.echo(f"    Non-canonical splice: {nc.call}")
    click.echo(f"      Canonical intron fraction: {nc.canonical_pct:.2f}%")
    click.echo(f"      Top non-canonical dinuc: {nc.top_non_canonical_dinuc}")


def _print_top_clusters(soft: SoftClipStats) -> None:
    if not soft.top_clusters:
        return
    click.echo("\n  top-15 clusters (key | n_loci | n_reads | example):")
    for key, n_loci, n_reads in soft.top_clusters:
        example = soft.cluster_examples.get(key, "")
        ex = example if len(example) <= 40 else example[:37] + "..."
        click.echo(f"    {key:<14}  {n_loci:>7,}  {n_reads:>9,}  {ex}")


def _print_locus_consistency(lc: LocusConsistencyStats) -> None:
    total = lc.n_loci_total
    click.echo("\n  locus consolidation:")
    click.echo(f"    singletons (n_clips=1):    {lc.n_loci_singleton:>10,} ({_pct(lc.n_loci_singleton, total)})")
    click.echo(f"    consistent (one motif):    {lc.n_loci_consistent:>10,} ({_pct(lc.n_loci_consistent, total)})")
    click.echo(f"    inconsistent:              {lc.n_loci_inconsistent:>10,} ({_pct(lc.n_loci_inconsistent, total)})")
    click.echo(f"    short-only (all clips <K): {lc.n_loci_short_only:>10,} ({_pct(lc.n_loci_short_only, total)})")

    if lc.motif_share_histogram:
        click.echo("\n  motif-share histogram (consistent non-singleton loci):")
        denom = lc.n_loci_consistent
        for bucket in _BUCKET_ORDER:
            n = lc.motif_share_histogram.get(bucket, 0)
            if n == 0 and bucket not in lc.motif_share_histogram:
                continue
            label = {
                "1": "1     (locus-specific)",
                "2-10": "2-10",
                "11-100": "11-100",
                "101-1000": "101-1000",
                ">1000": ">1000 (SL-like)",
            }[bucket]
            click.echo(f"    {label:<22} {n:>10,} ({_pct(n, denom)})")

    if lc.deepest_loci:
        click.echo("\n  top-15 deepest loci (chrom:pos:side | n_clips | status | motif | share | longest_clip):")
        for row in lc.deepest_loci:
            motif = row.motif_key or "-"
            longest = row.longest_clip if len(row.longest_clip) <= 38 else row.longest_clip[:35] + "..."
            click.echo(
                f"    {row.chrom}:{row.pos}:{row.side}  "
                f"{row.n_clips:>6,}  {row.status:<12}  {motif:<14}  "
                f"{row.motif_share:>7,}  {longest}"
            )


def _write_outputs(
    out_dir: Path, report: DiagnosticReport, verdict: Verdict,
) -> None:
    json_path = out_dir / "softclip_diagnostic.json"
    with open(json_path, "w") as f:
        json.dump(to_summary_dict(report, verdict), f, indent=2)
    click.echo(f"\n  wrote {json_path}")

    clusters_path = out_dir / "softclip_diagnostic_clusters.tsv"
    soft = report.softclip
    sorted_clusters = sorted(
        soft.cluster_to_loci.items(),
        key=lambda kv: (-kv[1], -soft.cluster_to_reads.get(kv[0], 0), kv[0]),
    )
    with open(clusters_path, "w") as f:
        f.write("cluster_key\tn_loci\tn_reads\texample\n")
        for key, n_loci in sorted_clusters:
            f.write(
                f"{key}\t{n_loci}\t{soft.cluster_to_reads.get(key, 0)}\t"
                f"{soft.cluster_examples.get(key, '')}\n"
            )
    click.echo(f"  wrote {clusters_path} ({len(sorted_clusters):,} rows)")

    loci_path = out_dir / "softclip_diagnostic_loci.tsv"
    sorted_rows = sorted(report.locus_consistency.all_rows, key=lambda r: -r.n_clips)
    with open(loci_path, "w") as f:
        f.write("chrom\tpos\tside\tn_clips\tstatus\tmotif_key\tmotif_share\tlongest_clip\n")
        for r in sorted_rows:
            f.write(
                f"{r.chrom}\t{r.pos}\t{r.side}\t{r.n_clips}\t{r.status}\t"
                f"{r.motif_key}\t{r.motif_share}\t{r.longest_clip}\n"
            )
    click.echo(f"  wrote {loci_path} ({len(sorted_rows):,} non-singleton loci)")


def _run_one(label: str, bam: Path, genome: Path, *, min_clip_len: int,
             cluster_key_len: int, min_mapq: int,
             hamming_tolerance: int, cluster_hamming_tolerance: int,
             min_consistency_fraction: float,
             consensus_min_majority_fraction: float) -> None:
    click.echo(f"\nDiagnosing {label}…")
    click.echo(f"  bam:    {bam}")
    click.echo(f"  genome: {genome}")
    click.echo(
        f"  per-locus H={hamming_tolerance}, "
        f"cross-locus H={cluster_hamming_tolerance}, "
        f"min-consistency-fraction={min_consistency_fraction}, "
        f"consensus-majority={consensus_min_majority_fraction}"
    )
    if not bam.exists():
        raise click.ClickException(f"BAM not found: {bam}")
    if not genome.exists():
        raise click.ClickException(f"Genome not found: {genome}")

    t0 = time.time()
    report = diagnose_bam(
        bam, genome,
        DiagnosticParams(
            min_clip_len=min_clip_len,
            cluster_key_len=cluster_key_len,
            min_mapq=min_mapq,
            hamming_tolerance=hamming_tolerance,
            cluster_hamming_tolerance=cluster_hamming_tolerance,
            min_consistency_fraction=min_consistency_fraction,
            consensus_min_majority_fraction=consensus_min_majority_fraction,
        ),
    )
    click.echo(f"  scan complete in {time.time() - t0:.1f}s")

    _print_summary(label, report.softclip, report.intron)
    _print_top_clusters(report.softclip)
    _print_locus_consistency(report.locus_consistency)
    verdict = compute_verdict(report)
    _print_verdict(verdict)
    _write_outputs(bam.parent, report, verdict)


@click.group()
def cli() -> None:
    """Soft-clip + intron diagnostic experiments."""


@cli.command("one")
@click.option("--bam", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--genome", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--label", required=True, type=str, help="Dataset label for the report header.")
@click.option("--min-clip-len", default=8, show_default=True, type=int)
@click.option("--cluster-key-len", default=12, show_default=True, type=int)
@click.option("--min-mapq", default=20, show_default=True, type=int)
@click.option("--hamming-tolerance", default=2, show_default=True, type=int,
              help="Per-locus Hamming tolerance: clip keys within this distance of the locus's dominant key consolidate.")
@click.option("--cluster-hamming-tolerance", default=1, show_default=True, type=int,
              help="Cross-locus canonicalization: raw cluster keys within this distance of an existing seed are absorbed into it.")
@click.option("--min-consistency-fraction", default=0.95, show_default=True, type=float,
              help="Fraction of tracked long-clip reads at a locus that must be within --hamming-tolerance of the dominant key for the locus to be marked consistent.")
@click.option("--consensus-min-majority-fraction", default=0.6, show_default=True, type=float,
              help="Per-column majority-vote threshold for the cluster consensus. Lower extends the consensus into noisier columns; higher terminates sooner.")
def cmd_one(
    bam: Path, genome: Path, label: str,
    min_clip_len: int, cluster_key_len: int, min_mapq: int,
    hamming_tolerance: int, cluster_hamming_tolerance: int,
    min_consistency_fraction: float,
    consensus_min_majority_fraction: float,
) -> None:
    """Run the diagnostic against one BAM + genome pair."""
    _run_one(
        label, bam, genome,
        min_clip_len=min_clip_len,
        cluster_key_len=cluster_key_len,
        min_mapq=min_mapq,
        hamming_tolerance=hamming_tolerance,
        cluster_hamming_tolerance=cluster_hamming_tolerance,
        min_consistency_fraction=min_consistency_fraction,
        consensus_min_majority_fraction=consensus_min_majority_fraction,
    )


@cli.command("run-all")
@click.option("--min-clip-len", default=8, show_default=True, type=int)
@click.option("--cluster-key-len", default=12, show_default=True, type=int)
@click.option("--min-mapq", default=20, show_default=True, type=int)
@click.option("--hamming-tolerance", default=2, show_default=True, type=int)
@click.option("--cluster-hamming-tolerance", default=1, show_default=True, type=int)
@click.option("--min-consistency-fraction", default=0.95, show_default=True, type=float)
@click.option("--consensus-min-majority-fraction", default=0.6, show_default=True, type=float)
def cmd_run_all(
    min_clip_len: int, cluster_key_len: int, min_mapq: int,
    hamming_tolerance: int, cluster_hamming_tolerance: int,
    min_consistency_fraction: float,
    consensus_min_majority_fraction: float,
) -> None:
    """Diagnose both hemistasia and flectonema using the known data paths."""
    for label, (bam, genome) in KNOWN.items():
        _run_one(
            label, bam, genome,
            min_clip_len=min_clip_len,
            cluster_key_len=cluster_key_len,
            min_mapq=min_mapq,
            hamming_tolerance=hamming_tolerance,
            cluster_hamming_tolerance=cluster_hamming_tolerance,
            min_consistency_fraction=min_consistency_fraction,
            consensus_min_majority_fraction=consensus_min_majority_fraction,
        )


if __name__ == "__main__":
    cli()
