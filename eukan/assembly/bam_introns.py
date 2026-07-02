"""Split spliced alignments at over-long introns (a max-intron BAM filter).

segemehl has no native maximum-intron parameter, so its read→genome BAM can carry
arbitrarily long N-CIGAR introns. Fed to StringTie unchanged those bridge distant
loci into one fused transcript. This module rewrites such a BAM so every alignment
spanning an intron longer than ``max_intron_len`` is *split* into independent
single-end pieces at that intron — the long junction never reaches StringTie, while
short introns (real splice junctions) and per-piece coverage are preserved.

It is a deliberately narrow tool: it does not touch the shared read BAM that SL
acceptor detection and the non-canonical diagnostic read (they must see the true
alignments). The caller (:func:`eukan.assembly.stringtie.run_stringtie`) bounds a
StringTie-only copy. STAR's BAM is already bounded by ``--alignIntronMax``, so only
the segemehl path needs this.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from eukan.assembly.bam_diagnostic import _CIGAR_N, _CIGAR_REF_CONSUMING
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd

log = get_logger(__name__)

# CIGAR ops, mirroring pysam/SAM op codes (see bam_diagnostic for the shared set).
_CIGAR_INSERTION = 1
_CIGAR_SOFT_CLIP = 4
_CIGAR_HARD_CLIP = 5
_CIGAR_M = 0
_CIGAR_EQ = 7
_CIGAR_X = 8
# Query (read) consuming ops.
_CIGAR_QUERY_CONSUMING = frozenset(
    [_CIGAR_M, _CIGAR_INSERTION, _CIGAR_SOFT_CLIP, _CIGAR_EQ, _CIGAR_X]
)
# Ops that align a query base to the reference (a piece must contain ≥1 to be real).
_CIGAR_ALIGNED = frozenset([_CIGAR_M, _CIGAR_EQ, _CIGAR_X])

# Pairing/multiplicity flag bits cleared when emitting a split piece as an
# independent single-end primary read (the reverse-strand bit 0x10 is kept).
_PAIRING_FLAGS = 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80 | 0x100 | 0x800

# Tags invalidated by re-coordinating/clipping a piece; XS (strand) and RG are kept.
_STALE_TAGS = frozenset(["NM", "MD", "SA", "NH", "HI", "nM", "jM", "jI"])


def _split_read(
    read: pysam.AlignedSegment, header: pysam.AlignmentHeader, max_intron_len: int
) -> list[pysam.AlignedSegment] | None:
    """Split *read* into pieces at every N op > *max_intron_len*; ``None`` if none.

    Each piece keeps the full SEQ/QUAL with the out-of-piece bases soft-clipped, so
    the query-consuming CIGAR lengths still sum to ``len(SEQ)``. Pieces are emitted
    as independent single-end primary reads (distinct ``<name>.<i>`` query names) so
    StringTie counts each piece's local coverage and short-intron junctions without
    bridging the removed long intron.
    """
    cigar = read.cigartuples
    if cigar is None or max_intron_len <= 0:
        return None
    # Hard-clipped / secondary / supplementary reads are passed through untouched:
    # splitting them risks corrupting the SEQ accounting or duplicating extra loci.
    if read.is_secondary or read.is_supplementary:
        return None
    if any(op == _CIGAR_HARD_CLIP for op, _ in cigar):
        return None
    if not any(op == _CIGAR_N and length > max_intron_len for op, length in cigar):
        return None

    total_query = sum(length for op, length in cigar if op in _CIGAR_QUERY_CONSUMING)

    segments: list[dict] = []
    cur: dict | None = None
    ref_pos = read.reference_start
    q_pos = 0
    for op, length in cigar:
        if op == _CIGAR_N and length > max_intron_len:
            if cur is not None:
                segments.append(cur)
                cur = None
            ref_pos += length  # the long intron consumes reference but starts no piece
            continue
        if cur is None:
            cur = {"ref_start": ref_pos, "ops": [], "q_before": q_pos, "q_in": 0}
        cur["ops"].append((op, length))
        if op in _CIGAR_QUERY_CONSUMING:
            q_pos += length
            cur["q_in"] += length
        if op in _CIGAR_REF_CONSUMING:
            ref_pos += length
    if cur is not None:
        segments.append(cur)

    pieces: list[pysam.AlignedSegment] = []
    new_flag = read.flag & ~_PAIRING_FLAGS
    for seg in segments:
        if not any(op in _CIGAR_ALIGNED for op, _ in seg["ops"]):
            continue  # a soft-clip-only fragment beside the intron: nothing aligned
        q_after = total_query - seg["q_before"] - seg["q_in"]
        ops = [*seg["ops"]]
        if seg["q_before"] > 0:
            ops.insert(0, (_CIGAR_SOFT_CLIP, seg["q_before"]))
        if q_after > 0:
            ops.append((_CIGAR_SOFT_CLIP, q_after))
        pieces.append(_make_piece(read, header, new_flag, seg["ref_start"], ops, len(pieces)))
    return pieces or None


def _make_piece(
    read: pysam.AlignedSegment,
    header: pysam.AlignmentHeader,
    flag: int,
    ref_start: int,
    cigartuples: list[tuple[int, int]],
    index: int,
) -> pysam.AlignedSegment:
    """Build one split piece from *read* with a fresh name, flag, start, and CIGAR."""
    a = pysam.AlignedSegment(header)
    a.query_name = f"{read.query_name}.{index}"
    a.flag = flag
    a.reference_id = read.reference_id
    a.reference_start = ref_start
    a.mapping_quality = read.mapping_quality
    a.query_sequence = read.query_sequence
    a.query_qualities = read.query_qualities  # must follow query_sequence
    a.cigartuples = cigartuples
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0
    a.set_tags([t for t in read.get_tags() if t[0] not in _STALE_TAGS])
    return a


def split_long_introns(
    in_bam: Path, out_bam: Path, *, max_intron_len: int, num_cpu: int = 1
) -> int:
    """Rewrite *in_bam* into coordinate-sorted *out_bam*, splitting over-long introns.

    Returns the number of reads that were split. Reads without an over-long intron
    pass through unchanged. Pieces are written to a temp BAM then coordinate-sorted
    with ``samtools sort`` into *out_bam* (splitting can move a piece past later
    reads, so a re-sort is required). ``max_intron_len <= 0`` disables splitting
    (a plain re-sort copy).
    """
    tmp = out_bam.parent / f"{out_bam.name}.unsorted.bam"
    n_split = 0
    with pysam.AlignmentFile(str(in_bam), "rb") as inp:
        header = inp.header
        with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
            for read in inp:
                pieces = (
                    _split_read(read, header, max_intron_len)
                    if max_intron_len > 0
                    else None
                )
                if pieces is None:
                    out.write(read)
                else:
                    n_split += 1
                    for piece in pieces:
                        out.write(piece)
    run_cmd(
        ["samtools", "sort", "-@", str(num_cpu), "-o", out_bam.name, tmp.name],
        cwd=out_bam.parent,
    )
    tmp.unlink(missing_ok=True)
    return n_split
