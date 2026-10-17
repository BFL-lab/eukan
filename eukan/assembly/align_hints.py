"""Post-alignment processing of the RNA-seq read alignments (minimap2).

Derives the splice-junction and hint artifacts the rest of the pipeline
consumes, from an aligner's coordinate-sorted BAM, one job per function:

- :func:`sj_table_from_bam` — the STAR-format ``SJ.out.tab``
- :func:`analyze_splice_sites` — ``splice_site_summary.json``, which drives the
  AUGUSTUS non-canonical-splice allowance
- :func:`write_intron_hints` + :func:`write_coverage_hints` — ``hints_introns.gff``
  and ``hints_coverage.gff``, the AUGUSTUS / GeneMark hints

The soft-clip diagnostic over the same BAM (``softclip_diagnostic_summary.json``)
is :func:`eukan.assembly.bam_diagnostic.run_softclip_diagnostic`.

minimap2 emits no junction table, so :func:`sj_table_from_bam` builds the
``SJ.out.tab`` equivalent from the BAM's N-CIGAR junctions and the downstream
functions consume it unchanged; the artifacts keep the shape they had when STAR
wrote the table natively.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import pysam

from eukan.assembly.bam_reads import dinucleotide, walk_introns
from eukan.infra.artifacts import Artifact
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.utils import concat_files

log = get_logger(__name__)


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

    For aligners without a native junction table (minimap2): walk primary
    alignments, tally per-junction unique (``NH==1``) / multi (``NH>1``) read
    support, classify the donor/acceptor dinucleotide into a STAR motif code +
    strand, and emit the 9-column ``SJ.out.tab`` STAR would. Reuses the BAM-walk
    primitives validated by the splice-wobble analysis. ``max_intron <= 0``
    disables the upper length filter. Returns the written path.

    Output columns match STAR's OutSJ.cpp (1-based, inclusive intron bounds):
    ``chrom, intron_start, intron_end, strand, motif, annotated, n_unique,
    n_multi, max_overhang``.
    """
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
            for chrom, start, end in walk_introns(read):
                ilen = end - start
                if ilen < min_intron or (0 < max_intron < ilen):
                    continue
                counts[(chrom, start, end)][slot] += 1

    out = wd / out_name
    with ContigIndex(genome) as contigs, open(out, "w") as f:
        for chrom, start, end in sorted(counts):
            n_unique, n_multi = counts[(chrom, start, end)]
            motif, strand = _motif_and_strand(
                dinucleotide(contigs, chrom, start, end)
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


def write_coverage_hints(bam: Path, wd: Path) -> None:
    """Write ``hints_coverage.gff``: AUGUSTUS exonpart hints from the BAM's per-strand
    read coverage (samtools view -> bam2wig -> wig2hints.pl)."""
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
