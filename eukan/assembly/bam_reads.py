"""Per-read BAM primitives shared by the assembly steps' BAM walks.

Primary-alignment filtering, mRNA-oriented soft-clip extraction, intron
walking and splice-site dinucleotide lookup, used by the soft-clip diagnostic
(:mod:`eukan.assembly.bam_diagnostic`), the junction table
(:mod:`eukan.assembly.align_hints`), poly-A characterization
(:mod:`eukan.assembly.polya`), SL acceptor detection, the SL cut and strand
correction. CIGAR op codes are pysam's (``pysam.CMATCH``, ``pysam.CREF_SKIP``,
...).
"""

from __future__ import annotations

from collections.abc import Iterable

import pysam

from eukan.infra.genome import reverse_complement

# CIGAR ops that advance the reference position (N included) and the read position.
REF_CONSUMING = frozenset(
    [pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF]
)
QUERY_CONSUMING = frozenset(
    [pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF]
)


def iter_primary_alignments(
    bam: pysam.AlignmentFile, *, min_mapq: int,
) -> Iterable[pysam.AlignedSegment]:
    """Yield primary, mapped, MAPQ-passing reads with a usable cigar+seq."""
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.cigartuples is None or read.query_sequence is None:
            continue
        yield read


def extract_clips(
    read: pysam.AlignedSegment, min_clip_len: int,
) -> Iterable[tuple[str, str, int]]:
    """Yield ``(side, mRNA_orient_seq, anchor_pos)`` for each soft-clip end.

    BAM stores sequences in reference orientation, so a leading clip on a
    reverse-mapped read is the 3' end of the mRNA (and its sequence is
    the reverse complement of the original sequencer read). We
    normalize: the side label is in mRNA orientation, and the sequence is
    flipped to match.
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    ref_end = read.reference_end
    if cigar is None or seq is None or ref_end is None:
        return

    op0, len0 = cigar[0]
    if op0 == pysam.CSOFT_CLIP and len0 >= min_clip_len:
        clip_seq = seq[:len0]
        anchor_pos = read.reference_start
        if read.is_reverse:
            yield "3p", reverse_complement(clip_seq), anchor_pos
        else:
            yield "5p", clip_seq, anchor_pos

    op_last, len_last = cigar[-1]
    if op_last == pysam.CSOFT_CLIP and len_last >= min_clip_len:
        clip_seq = seq[-len_last:]
        anchor_pos = ref_end - 1
        if read.is_reverse:
            yield "5p", reverse_complement(clip_seq), anchor_pos
        else:
            yield "3p", clip_seq, anchor_pos


def walk_introns(
    read: pysam.AlignedSegment,
) -> Iterable[tuple[str, int, int]]:
    """Yield ``(chrom, intron_start, intron_end_exclusive)`` for each N op."""
    chrom = read.reference_name
    cigar = read.cigartuples
    if chrom is None or cigar is None:
        return
    pos = read.reference_start
    for op, length in cigar:
        if op == pysam.CREF_SKIP:
            yield chrom, pos, pos + length
        if op in REF_CONSUMING:
            pos += length


def dinucleotide(
    contigs, chrom: str, intron_start: int, intron_end: int,
) -> str | None:
    """Look up the donor/acceptor dinucleotide pair on the genome + strand.

    Returns ``"XX-YY"`` (e.g. ``"GT-AG"`` for + strand canonical,
    ``"CT-AC"`` for - strand canonical, matching STAR's motif
    convention), or ``None`` when the contig or coordinates fall outside
    the genome index. The reverse-strand canonical pair is reported
    distinctly so a downstream check can detect canonical sites without
    needing the read's transcribed strand.
    """
    rec = contigs.get(chrom)
    if rec is None or rec.seq is None:
        return None
    genome_seq = rec.seq
    if intron_end > len(genome_seq) or intron_start < 0:
        return None
    donor = str(genome_seq[intron_start : intron_start + 2]).upper()
    acceptor = str(genome_seq[intron_end - 2 : intron_end]).upper()
    return f"{donor}-{acceptor}"
