"""Post-alignment processing shared across aligners (STAR, segemehl).

Derives the splice-junction / hint / diagnostic artifacts the rest of the
pipeline consumes, from an aligner's coordinate-sorted BAM plus a STAR-format
``SJ.out.tab``:

- ``splice_site_summary.json`` — drives AUGUSTUS non-canonical-splice allowance
- ``hints_introns.gff`` + ``hints_coverage.gff`` — AUGUSTUS / GeneMark hints
- ``softclip_diagnostic_summary.json`` — trans-splicing / non-canonical verdict

STAR emits ``SJ.out.tab`` natively; aligners that don't (segemehl) build an
equivalent from the BAM's N-CIGAR junctions via :func:`sj_table_from_bam`, then
reuse the same downstream functions so the artifacts are byte-for-byte the same
shape STAR produced.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from eukan.assembly.sl_depletion import is_adapter
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.utils import concat_files

if TYPE_CHECKING:
    from eukan.assembly.bam_diagnostic import TransSplicingCall
    from eukan.assembly.polya import PolyAStats

log = get_logger(__name__)

# Read-BAM poly-A is characterized at the same MAPQ floor diagnose_bam uses, so a
# backfill (poly-A-only) pass over the read BAM matches the full-walk numbers.
_DIAGNOSTIC_MIN_MAPQ = 20


def _log_read_polya(pa: PolyAStats, bam_name: str) -> None:
    """INFO-log the read-BAM poly-A soft-clip tally (no-op when there are none)."""
    if pa.n_polya or pa.n_polyt:
        log.info(
            "Poly-A soft-clips in %s: %d poly-A (3') + %d poly-T (5') of %d clips "
            "(%.3f%% poly-A; mean %.1f bp, max %d bp).",
            bam_name, pa.n_polya, pa.n_polyt, pa.n_clips_examined,
            pa.polya_pct_of_clips, pa.polya_mean_len, pa.polya_len_max,
        )

# STAR motif codes → canonical/semi-canonical splice site dinucleotide pairs.
_MOTIF_NAMES: dict[int, str] = {
    1: "GT-AG",
    2: "CT-AC",
    3: "GC-AG",
    4: "CT-GC",
    5: "AT-AC",
    6: "GT-AT",
}
# Inverse: dinucleotide pair → (STAR motif code, strand). strand 1 = '+',
# 2 = '-', 0 = undefined (non-canonical).
_MOTIF_CODE_STRAND: dict[str, tuple[int, int]] = {
    "GT-AG": (1, 1), "CT-AC": (2, 2), "GC-AG": (3, 1),
    "CT-GC": (4, 2), "AT-AC": (5, 1), "GT-AT": (6, 2),
}


def _motif_and_strand(dinuc: str | None) -> tuple[int, int]:
    """Map a donor-acceptor pair to a STAR ``(motif_code, strand)`` tuple."""
    return _MOTIF_CODE_STRAND.get(dinuc or "", (0, 0))


def sj_table_from_bam(
    bam: Path,
    genome: Path,
    wd: Path,
    *,
    min_intron: int,
    max_intron: int,
    out_name: str,
) -> Path:
    """Derive a STAR-format ``SJ.out.tab`` from a BAM's N-CIGAR junctions.

    For aligners without a native junction table (segemehl): walk primary
    alignments, tally per-junction unique (``NH==1``) / multi (``NH>1``) read
    support, classify the donor/acceptor dinucleotide into a STAR motif code +
    strand, and emit the 9-column ``SJ.out.tab`` STAR would. Reuses the BAM-walk
    primitives validated by the splice-wobble analysis. ``max_intron <= 0``
    disables the upper length filter. Returns the written path.

    Output columns match STAR's OutSJ.cpp (1-based, inclusive intron bounds):
    ``chrom, intron_start, intron_end, strand, motif, annotated, n_unique,
    n_multi, max_overhang``.
    """
    import pysam

    from eukan.assembly.bam_diagnostic import _dinucleotide, _walk_introns
    from eukan.infra.genome import ContigIndex

    # (chrom, start, end) 0-based half-open intron -> [n_unique, n_multi]
    counts: dict[tuple[str, int, int], list[int]] = defaultdict(lambda: [0, 0])
    with pysam.AlignmentFile(str(bam), "rb") as af:
        for read in af:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.cigartuples is None:
                continue
            nh = read.get_tag("NH") if read.has_tag("NH") else 1
            slot = 0 if nh == 1 else 1
            for chrom, start, end in _walk_introns(read):
                ilen = end - start
                if ilen < min_intron or (0 < max_intron < ilen):
                    continue
                counts[(chrom, start, end)][slot] += 1

    out = wd / out_name
    with ContigIndex(genome) as contigs, open(out, "w") as f:
        for chrom, start, end in sorted(counts):
            n_unique, n_multi = counts[(chrom, start, end)]
            motif, strand = _motif_and_strand(
                _dinucleotide(contigs, chrom, start, end)
            )
            # STAR SJ.out.tab is 1-based inclusive: col2 = start+1, col3 = end.
            f.write(
                f"{chrom}\t{start + 1}\t{end}\t{strand}\t{motif}\t0\t"
                f"{n_unique}\t{n_multi}\t0\n"
            )
    return out


def analyze_splice_sites(sj_file: Path, genome: Path, wd: Path) -> None:
    """Extract splice site dinucleotides from a junction table and summarize.

    For each junction in an SJ.out.tab, extracts the donor and acceptor
    dinucleotides from the genome FASTA.  Writes ``splice_site_summary.json``
    with per-type counts and read support.

    SJ.out.tab columns (STAR OutSJ.cpp convention):
      col2 = first base of intron (1-based)
      col3 = last base of intron (1-based)
      col5 = motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, ...)
      col7 = unique reads, col8 = multi-mapping reads
    """
    from eukan.infra.genome import ContigIndex

    # Tally: splice_type → {"count": int, "unique_reads": int}
    tallies: dict[str, dict[str, int]] = defaultdict(
        lambda: {"count": 0, "unique_reads": 0}
    )

    with ContigIndex(genome) as contigs, open(sj_file) as fin:
        reader = csv.reader(fin, delimiter="\t")
        for row in reader:
            chrom = row[0]
            intron_start = int(row[1])  # 1-based, first base of intron
            intron_end = int(row[2])    # 1-based, last base of intron
            motif = int(row[4])
            unique = int(row[6])

            if motif != 0:
                # Use the motif classification for canonical/semi-canonical.
                splice_type = _MOTIF_NAMES[motif]
            else:
                # Extract actual dinucleotides from the genome.
                seq = contigs.get(chrom)
                if seq is None or seq.seq is None or intron_end > len(seq):
                    splice_type = "unknown"
                else:
                    genome_seq = seq.seq
                    donor = str(genome_seq[intron_start - 1 : intron_start + 1]).upper()
                    acceptor = str(genome_seq[intron_end - 2 : intron_end]).upper()
                    splice_type = f"{donor}-{acceptor}"

            tallies[splice_type]["count"] += 1
            tallies[splice_type]["unique_reads"] += unique

    summary = dict(sorted(tallies.items(), key=lambda kv: -kv[1]["count"]))
    with open(wd / Artifact.SPLICE_SUMMARY, "w") as f:
        json.dump(summary, f, indent=2)

    # Log summary (skip canonical — they dominate).
    for stype, counts in summary.items():
        if stype in ("GT-AG", "CT-AC"):
            continue
        log.info(
            "Splice sites (%s): %d junctions, %d unique reads",
            stype, counts["count"], counts["unique_reads"],
        )


def _log_trans_splicing_verdict(ts: TransSplicingCall) -> None:
    """Log the trans-splicing verdict, flagging when the dominant soft-clip motif
    is actually residual sequencing adapter rather than a genuine spliced leader.

    Reads aren't adapter-trimmed upstream, so Illumina/Nextera read-through can be
    the dominant soft-clip cluster and masquerade as trans-splicing. Surfacing it
    here (in the soft-clip analysis) tells the user to adapter-trim rather than
    trust a phantom spliced leader; SL detection separately excludes it.
    """
    label = ts.top_non_trivial_cluster_consensus or ts.top_non_trivial_cluster_key
    if label and (
        is_adapter(ts.top_non_trivial_cluster_consensus)
        or is_adapter(ts.top_non_trivial_cluster_key)
    ):
        note = (
            f"; the {ts.call} trans-splicing signal is attributable to it"
            if ts.call in ("STRONG", "MODERATE")
            else ""
        )
        log.warning(
            "Soft-clip analysis: the dominant soft-clip motif %s matches a known "
            "sequencing adapter%s. This is residual adapter read-through, not a "
            "spliced leader; adapter-trim the reads before assembly. SL detection "
            "excludes it by default (--no-sl-adapter-filter to override).",
            label, note,
        )
        return
    if ts.call in ("STRONG", "MODERATE"):
        log.warning(
            "Trans-splicing signal %s: top motif %s spans %d loci (%d reads). "
            "Reads may need splice-leader trimming before annotation.",
            ts.call,
            label,
            ts.top_non_trivial_cluster_n_loci,
            ts.top_non_trivial_cluster_n_reads,
        )
    else:
        log.info("Trans-splicing signal: ABSENT")


def run_softclip_diagnostic(bam: Path, genome: Path, wd: Path) -> None:
    """Walk the aligner BAM for soft-clip + intron motifs and log a verdict.

    Idempotent: if the summary JSON already exists, this is a no-op. The
    verdict surfaces trans-splicing and non-canonical splice prevalence so the
    user knows whether downstream gene prediction will need special handling
    (read pre-processing for trans-splicing; ``--splice-permissive`` for
    non-canonical splice landscapes).
    """
    from eukan.assembly.bam_diagnostic import (
        compute_verdict,
        diagnose_bam,
        to_summary_dict,
    )
    from eukan.assembly.polya import (
        characterize_polya_bam,
        has_section,
        stats_to_dict,
        write_polya_section,
    )

    summary_path = wd / Artifact.SOFTCLIP_DIAGNOSTIC.value
    if not bam.exists():
        return
    if summary_path.exists():
        log.info("Soft-clip diagnostic already produced %s, skipping", summary_path.name)
        # The poly-A "reads" section is a later addition and is written separately
        # from the SL summary, so backfill it (cheap poly-A-only pass over the read
        # BAM) when an existing or pre-feature summary would otherwise skip it — else
        # the section is silently never produced on resume / force / upgrade-in-place.
        if not has_section(wd, "reads"):
            pa = characterize_polya_bam(bam, "reads", min_mapq=_DIAGNOSTIC_MIN_MAPQ)
            write_polya_section(wd, "reads", stats_to_dict(pa))
            _log_read_polya(pa, bam.name)
        return

    log.info("Running soft-clip / intron diagnostic over %s...", bam.name)
    report = diagnose_bam(bam, genome)
    verdict = compute_verdict(report)

    with open(summary_path, "w") as f:
        json.dump(to_summary_dict(report, verdict), f, indent=2)

    # Poly-A characterization of the read soft-clips (separate from the SL verdict),
    # piggy-backing on the single diagnose_bam walk just performed.
    write_polya_section(wd, "reads", stats_to_dict(report.polya))
    _log_read_polya(report.polya, bam.name)

    _log_trans_splicing_verdict(verdict.trans_splicing)

    nc = verdict.non_canonical_splice
    if nc.call in ("EXTENSIVE", "MODERATE"):
        log.warning(
            "Non-canonical splice signal %s: canonical fraction %.2f%% "
            "(top non-canonical %s). Consider --splice-permissive on the assemble step.",
            nc.call, nc.canonical_pct, nc.top_non_canonical_dinuc,
        )
    else:
        log.info("Canonical splice site usage typical: %.2f%%", nc.canonical_pct)


def generate_rnaseq_hints(
    sj_file: Path,
    bam: Path,
    genome: Path,
    wd: Path,
    *,
    diagnose: bool = True,
    source_label: str = "RNASEQ",
) -> None:
    """Generate AUGUSTUS hints + splice summary from an aligner's SJ + BAM.

    ``source_label`` is the GFF column-2 source written on intron hints
    (cosmetic; AUGUSTUS keys on the ``src=`` attribute). Coverage hints are
    derived from the BAM and are aligner-agnostic.
    """
    if sj_file.exists():
        analyze_splice_sites(sj_file, genome, wd)
    if diagnose:
        run_softclip_diagnostic(bam, genome, wd)
    if sj_file.exists():
        write_intron_hints(sj_file, wd, source_label)

    # Generate coverage hints from the BAM.
    if bam.exists():
        run_cmd(
            ["samtools", "view", "-b", "-f", "0x10", str(bam)],
            cwd=wd, out_file="cov_reverse.bam", binary=True,
        )
        run_cmd(
            ["samtools", "view", "-b", "-F", "0x10", str(bam)],
            cwd=wd, out_file="cov_forward.bam", binary=True,
        )

        for direction, _strand, wig in [
            ("cov_reverse.bam", "-", "minus.wig"),
            ("cov_forward.bam", "+", "plus.wig"),
        ]:
            run_cmd(["bam2wig", direction], cwd=wd, out_file=wig)

        # wig2hints.pl reads from stdin, so we pipe the wig file in.
        for wig_file, strand, hints_file in [
            ("minus.wig", "-", "hints.ep.minus.gff"),
            ("plus.wig", "+", "hints.ep.plus.gff"),
        ]:
            _run_wig2hints(wd, wig_file, strand, hints_file)

        concat_files(
            [wd / hf for hf in ["hints.ep.minus.gff", "hints.ep.plus.gff"] if (wd / hf).exists()],
            wd / "hints_coverage.gff",
        )

        for f in ["cov_reverse.bam", "cov_forward.bam", "minus.wig", "plus.wig",
                  "hints.ep.minus.gff", "hints.ep.plus.gff"]:
            (wd / f).unlink(missing_ok=True)


def write_intron_hints(sj_file: Path, wd: Path, source_label: str) -> None:
    """Write ``hints_introns.gff`` from a STAR-format SJ.out.tab.

    GFF column 2 is ``source_label`` (cosmetic; AUGUSTUS keys on the ``src=``
    attribute). Score and ``mult=`` are unique + multi reads.
    """
    strand_map = {"0": ".", "1": "+", "2": "-"}
    with open(sj_file) as fin, open(wd / "hints_introns.gff", "w") as fout:
        reader = csv.reader(fin, delimiter="\t")
        for row in reader:
            chrom, start, end = row[0], row[1], row[2]
            strand = strand_map.get(row[3], ".")
            unique = int(row[6]) + int(row[7])
            fout.write(
                f"{chrom}\t{source_label}\tintron\t{start}\t{end}\t{unique}\t"
                f"{strand}\t.\tmult={unique};pri=4;src=E\n"
            )


def _run_wig2hints(wd: Path, wig_file: str, strand: str, out_file: str) -> None:
    """Run wig2hints.pl, reading the wig from stdin and writing GFF to stdout."""
    run_cmd(
        [
            "wig2hints.pl",
            "--width=10", "--margin=10", "--minthresh=2",
            "--minscore=4", "--prune=0.1", "--src=W",
            "--type=exonpart", "--radius=4.5", "--pri=4",
            f"--strand={strand}",
        ],
        cwd=wd,
        in_file=wig_file,
        out_file=out_file,
    )
