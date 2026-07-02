"""Genomic spliced-leader cut: split transcript models at trans-splice acceptors.

Reuses the jaccard GFF3/GTF split logic (:mod:`eukan.assembly.jaccard`) — the cut
is the same exon-segregating split, fed SL acceptor *genomic* coordinates
(:mod:`eukan.assembly.sl_acceptors`) instead of read-coverage troughs. For each
transcript whose exons contain a same-strand acceptor, the model is cut so the
mature mRNA begins at the acceptor. An SL site imposes its strand on an
otherwise-unstranded (``.``) transcript.

Inputs cut: the two Trinity transcript→genome model sets (de novo + genome-guided),
already max-intron-split upstream by :mod:`eukan.assembly.max_intron` into
``{stem}.maxintron.gff3``. (The over-long-intron split is a separate, composable
step — this module now does the SL cut only.) All outputs are genome-coordinate
GFF3 that ``combinr assemble`` ingests directly; the cut **streams** one transcript
at a time — never materialising the full set — so a genome-wide model set stays
bounded.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pysam

from eukan.assembly.bam_diagnostic import (
    _CIGAR_D,
    _CIGAR_EQ,
    _CIGAR_M,
    _CIGAR_N,
    _CIGAR_X,
)
from eukan.assembly.jaccard import (
    _iter_transcript_models,
    _split_transcript,
    _Tx,
    _write_tx,
)
from eukan.assembly.sl_acceptors import AcceptorSite, load_sl_acceptors
from eukan.assembly.tracks import mapped_transcript_stems
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# Ref-consuming CIGAR ops that stay within one exon block (N splits exons).
_EXON_REF = frozenset([_CIGAR_M, _CIGAR_D, _CIGAR_EQ, _CIGAR_X])

_DENOVO_BAMS = ("trinity-denovo.genome.bam", "trinity-gg.genome.bam")
_GENOME_BAM_SUFFIX = ".genome.bam"


def _alignment_exons(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    """Genomic 1-based exon blocks for one alignment (split at CIGAR N)."""
    cigar = read.cigartuples
    if cigar is None:
        return []
    blocks: list[tuple[int, int]] = []
    cur_start = read.reference_start  # 0-based
    pos = cur_start
    for op, length in cigar:
        if op == _CIGAR_N:
            blocks.append((cur_start + 1, pos))
            pos += length
            cur_start = pos
        elif op in _EXON_REF:
            pos += length
    blocks.append((cur_start + 1, pos))
    return [b for b in blocks if b[1] >= b[0]]


def bam_to_transcript_gff3(bam_path: Path, out_gff: Path, source: str) -> int:
    """Stream one mRNA per mapped alignment (exons from CIGAR N) to *out_gff*.

    Each alignment gets a unique mRNA id (``<qname>.m<n>``) — an aligner reporting
    a transcript at every locus it maps to would otherwise collide ids, and combinr
    groups exons by ``Parent``, silently fusing distinct loci. Secondary alignments
    are kept (extra loci of multi-copy genes); supplementary/unmapped are skipped.
    Writes incrementally so a genome-wide BAM never builds a full in-memory list.
    """
    seen: Counter[str] = Counter()
    n = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for read in bam:
            if read.is_unmapped or read.is_supplementary or read.cigartuples is None:
                continue
            chrom = read.reference_name
            if chrom is None:
                continue
            exons = _alignment_exons(read)
            if not exons:
                continue
            qname = read.query_name or "tx"
            seen[qname] += 1
            strand = "-" if read.is_reverse else "+"
            _write_tx(fh, _Tx(f"{qname}.m{seen[qname]}", chrom, strand, source, exons))
            n += 1
    return n


def _project_genomic_to_spliced(
    exons_5to3: list[tuple[int, int]], strand: str, pos: int
) -> int | None:
    """1-based spliced offset of genomic *pos*, or ``None`` if intronic/outside.

    Inverse of :func:`jaccard._partition_exons`'s span accumulation: walk exons in
    5'→3' order (ascending genomic for ``+``, descending for ``-``) summing
    spliced length until *pos* lands inside an exon.
    """
    consumed = 0
    for gstart, gend in exons_5to3:
        if gstart <= pos <= gend:
            if strand == "-":
                return consumed + (gend - pos) + 1
            return consumed + (pos - gstart) + 1
        consumed += gend - gstart + 1
    return None


def _cut_one(
    tx: _Tx, sites: list[AcceptorSite], min_segment: int
) -> list[_Tx] | None:
    """Split *tx* at same-strand SL acceptors; ``None`` if nothing to cut.

    A ``.``-strand transcript is oriented by the acceptors when they agree on a
    single strand (the SL imposes strand); conflicting strands skip the cut. The
    over-long-intron split is a separate, upstream step
    (:mod:`eukan.assembly.max_intron`), so models reaching here are already
    intron-bounded.
    """
    usable: list[AcceptorSite] = []
    strand = tx.strand
    if sites:
        if tx.strand in ("+", "-"):
            usable = [s for s in sites if s.strand == tx.strand]
        else:
            strands = {s.strand for s in sites}
            if len(strands) == 1:
                strand = next(iter(strands))
                usable = sites
            # conflicting SL strands: no SL-imposed strand, so no cut

    exons_5to3 = tx.exons if strand != "-" else list(reversed(tx.exons))
    total_len = sum(e - s + 1 for s, e in tx.exons)

    clips: set[int] = set()
    for site in usable:
        off = _project_genomic_to_spliced(exons_5to3, strand, site.pos)
        # Cut *before* the acceptor base so the downstream piece starts at it.
        if off is not None and 1 < off <= total_len:
            clips.add(off - 1)

    valid = sorted(c for c in clips if 0 < c < total_len)
    if not valid:
        return None

    oriented = (
        tx if tx.strand == strand else _Tx(tx.tid, tx.chrom, strand, tx.source, tx.exons)
    )
    return _split_transcript(oriented, valid, min_segment)


def cut_models_at_sl(
    gff_or_gtf: str | Path,
    sites: list[AcceptorSite],
    out_gff: str | Path,
    *,
    min_segment: int,
) -> int:
    """Cut every transcript in *gff_or_gtf* at same-strand SL acceptors.

    Streams: reads one transcript model, cuts it, writes the result, and moves on —
    no full model list is held. Transcripts with no applicable acceptor pass through
    unchanged. (Over-long introns are split upstream by
    :mod:`eukan.assembly.max_intron`.) Returns the number of transcripts cut.
    """
    by_chrom: dict[str, list[AcceptorSite]] = {}
    for site in sites:
        by_chrom.setdefault(site.chrom, []).append(site)

    n_cut = 0
    n_conflict = 0
    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for tx in _iter_transcript_models(gff_or_gtf):
            span_lo, span_hi = tx.exons[0][0], tx.exons[-1][1]
            relevant = [s for s in by_chrom.get(tx.chrom, []) if span_lo <= s.pos <= span_hi]
            result = _cut_one(tx, relevant, min_segment)
            if result is None:
                # A multi-exon, stranded transcript whose only in-span SL acceptors
                # sit on the opposite strand: the splice-derived strand and the SL
                # disagree. Trust the introns — keep the strand, skip the cut.
                if (
                    relevant and len(tx.exons) > 1
                    and tx.strand in ("+", "-")
                    and all(s.strand != tx.strand for s in relevant)
                ):
                    n_conflict += 1
                _write_tx(fh, tx)
            else:
                for piece in result:
                    _write_tx(fh, piece)
                n_cut += 1
    if n_conflict:
        log.warning(
            "%d multi-exon transcript(s) had SL acceptors only on the opposite "
            "strand; kept the splice-derived strand and skipped those cuts.",
            n_conflict,
        )
    return n_cut


def run_sl_cut(config: AssemblyConfig) -> None:
    """Cut transcript models at SL trans-splice acceptors → ``{stem}.cut.gff3``.

    Reads each track's max-intron-split models (``{stem}.maxintron.gff3``, always
    written by :mod:`eukan.assembly.max_intron`) and cuts every transcript whose
    exons contain a same-strand SL acceptor so the mature mRNA begins at the
    acceptor. With no SL signal ``sl_acceptors.gff3`` is header-only, so this is a
    pass-through copy (and stays silent — sl_detect already reported the no-op).
    Outputs are the genome-coordinate GFF3 ``combinr assemble`` ingests directly;
    the cut streams one transcript at a time so a genome-wide model set stays bounded.
    """
    wd = config.work_dir
    acc_path = wd / Artifact.SL_ACCEPTORS.value
    sites = load_sl_acceptors(acc_path) if acc_path.exists() else []
    min_segment = config.min_sl_fragment

    for stem in mapped_transcript_stems():
        src = wd / f"{stem}.maxintron.gff3"
        if not src.exists():
            continue
        out = wd / f"{stem}.cut.gff3"
        n_cut = cut_models_at_sl(src, sites, out, min_segment=min_segment)
        # Only narrate when an SL signal actually drove a cut; a pass-through
        # (no acceptors) would just repeat sl_detect's no-op line per track.
        if sites:
            log.info("SL cut %s -> %s (%d transcripts cut).", src.name, out.name, n_cut)
