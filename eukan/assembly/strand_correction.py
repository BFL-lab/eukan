"""Homology-calibrated splice-strand correction for unstranded assemblies.

On a library mapped without ``-S/--strand-specific`` the strand is unknown.
The Trinity transcript->genome models then label each transcript
from a canonical ``GT-AG`` guess, or leave it ``.`` — wrong for organisms whose
true splice consensus is non-canonical, and unrecoverable from the BAM alone (the
unstranded junction set is strand-symmetric: every motif appears with its
reverse-complement twin).

This module breaks that symmetry with protein homology. It runs ``diamond blastx``
(``--strand both``) of the assembled transcripts against SwissProt; the query
*frame* of each hit gives the coding strand outright — a positive frame means the
labelled orientation is coding, a negative frame means its reverse complement is.
A hit therefore both confirms *or flips* the label (so an antisense-assembled de
novo contig is rescued, not just a forward-labelled one) and, read on that coding
strand, calibrates the organism's dominant splice consensus ``D``. Every other
multi-exon transcript is re-stranded by majority vote over its introns: a
``+``-genome donor-acceptor reading equal to ``D`` *or* to canonical ``GT-AG``
votes ``+``; equal to the RC-twin of ``D`` *or* to ``CT-AC`` votes ``-`` (the
canonical pair is always strand-informative, so it is counted even when ``D`` is
non-canonical). Mono-exonic and tied transcripts are left for the SL cut to orient.
Correction only ever rewrites the strand field
(coordinates are untouched), and combinr's FASTA emitter reverse-complements on
``-``, so a flipped strand yields correctly oriented evidence automatically.

The step is a no-op unless ``--uniprot`` is supplied *and* the library is
unstranded. It always emits ``<track>.genome.gff3`` (each Trinity transcript->genome
BAM converted to models) for the SL cut, and — when active — a
``<track>.genome.stranded.gff3`` per track plus a ``strand_correction.tsv`` audit
table.

Ported from ``strand_disambiguation.py`` (commit 52a8e63): the
homology-tool-agnostic hit parsing and the ``consensus_on_strand`` / ``introns_of``
splice-motif helpers. The per-locus *drop* is replaced by a per-transcript *flip*,
which fits the one-strand-per-transcript model output.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pysam

from eukan.assembly.bam_reads import REF_CONSUMING, dinucleotide
from eukan.assembly.diamond import coding_strand, run_blastx
from eukan.assembly.tracks import (
    GENOME_BAM_SUFFIX,
    RAW_GFF3,
    STRANDED_GFF3,
    TRANSCRIPT_SETS,
    mapped_transcript_stems,
)
from eukan.gff.transcripts import Transcript, read_transcripts, write_transcript, write_transcripts
from eukan.infra.genome import ContigIndex, reverse_complement
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_QUERY_FASTA = "strand_query.fasta"
_HITS_TSV = "strand_blastx.tsv"
_AUDIT_TSV = "strand_correction.tsv"

# Fallback consensus when too few confirmed introns calibrate one. The canonical
# pair is always strand-informative (GT-AG on +genome => + coding; its +genome
# RC-twin CT-AC => - coding), so the vote counts it alongside the learned dominant.
_CANONICAL = "GT-AG"
_CANONICAL_RC = "CT-AC"  # == _rc_swap(_CANONICAL)


# ---------------------------------------------------------------------------
# Pure helpers (ported / adapted from strand_disambiguation.py @ 52a8e63)
# ---------------------------------------------------------------------------


def parse_hits(path: Path) -> dict[str, int]:
    """Map each query id to the query frame of its best (highest-bitscore) hit.

    Tabular columns are ``qseqid sseqid bitscore qframe`` (diamond ``--outfmt 6``).
    The ``qframe`` sign encodes which query strand coded the hit: positive ⟹ the
    query as supplied (its labelled orientation) is coding; negative ⟹ its reverse
    complement is. diamond's ``--evalue`` has already filtered to significant hits,
    so any reported query is "confirmed". Comment/short rows skip.
    """
    best: dict[str, tuple[float, int]] = {}
    with path.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 4 or row[0].startswith("#"):
                continue
            try:
                bitscore, qframe = float(row[2]), int(row[3])
            except ValueError:
                continue
            qid = row[0]
            if qid not in best or bitscore > best[qid][0]:
                best[qid] = (bitscore, qframe)
    return {qid: qframe for qid, (_bs, qframe) in best.items()}


def introns_of(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """0-based ``[start, end)`` introns from 1-based inclusive exon blocks.

    Matches :func:`bam_reads.dinucleotide`'s convention: the first intron
    base is 0-based ``exon_end``; the exclusive end is 0-based ``next_start - 1``.
    """
    ex = sorted(exons)
    return [(ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)]


def consensus_on_strand(raw: str | None, strand: str) -> str | None:
    """Express a plus-genome donor-acceptor pair on the transcript's coding strand.

    On ``-`` the gene's real donor/acceptor are the reverse complement of the plus
    acceptor/donor, so ``"CT-AC"`` becomes ``"GT-AG"``. ``+`` and ``.`` (extracted
    forward) read as-is.
    """
    if raw is None:
        return None
    if strand == "-":
        return _rc_swap(raw)
    return raw


def _rc_swap(pair: str) -> str:
    """The reverse-complement twin of a donor-acceptor pair (``GT-AG`` <-> ``CT-AC``)."""
    donor, acceptor = pair.split("-")
    return f"{reverse_complement(acceptor)}-{reverse_complement(donor)}"


def _pick_consensus(tally: Counter[str], min_consensus: int) -> tuple[str, str]:
    """Return ``(dominant, rc_twin)`` coding-strand pair, falling back to GT-AG.

    Below *min_consensus* confirmed introns (or none at all) the calibration is too
    thin to trust, so the canonical ``GT-AG`` is used instead.
    """
    total = sum(tally.values())
    dominant = tally.most_common(1)[0][0] if total and total >= min_consensus else _CANONICAL
    return dominant, _rc_swap(dominant)


def _decide(
    tx: Transcript, hit_coding: str | None, dominant: str, rc: str, contigs: ContigIndex
) -> tuple[str, str]:
    """Return ``(new_strand, decision)`` for one transcript.

    A blastx hit gives the coding strand outright (*hit_coding*, derived from the
    query frame), so the label is set to it — homology can now *flip* an
    antisense-assembled contig, not just confirm a forward label (a ``.`` resolves
    to whichever strand coded the hit). Otherwise a multi-exon transcript is voted
    over its introns: a ``+``-genome donor-acceptor equal to *dominant* or to
    canonical ``GT-AG`` votes ``+``; equal to its RC-twin *rc* or to ``CT-AC`` votes
    ``-``; the majority wins. Mono-exonic / tied transcripts are left for the SL cut.
    """
    if hit_coding is not None:
        new = hit_coding
        if new == tx.strand:
            return new, "keep"
        return new, ("assign" if tx.strand not in ("+", "-") else "flip")
    introns = introns_of(tx.exons)
    if not introns:
        return tx.strand, "mono-exon"
    plus = minus = 0
    for istart, iend in introns:
        raw = dinucleotide(contigs, tx.chrom, istart, iend)
        if raw in (dominant, _CANONICAL):
            plus += 1
        elif raw in (rc, _CANONICAL_RC):
            minus += 1
    if plus == minus:
        return tx.strand, "ambiguous"
    new = "+" if plus > minus else "-"
    if new == tx.strand:
        return new, "keep"
    return new, ("assign" if tx.strand not in ("+", "-") else "flip")


# ---------------------------------------------------------------------------
# Trinity transcript->genome BAM -> models
# ---------------------------------------------------------------------------

# Ref-consuming CIGAR ops that stay within one exon block (N splits exons).
_EXON_REF = REF_CONSUMING - {pysam.CREF_SKIP}


def _alignment_exons(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    """Genomic 1-based exon blocks for one alignment (split at CIGAR N)."""
    cigar = read.cigartuples
    if cigar is None:
        return []
    blocks: list[tuple[int, int]] = []
    cur_start = read.reference_start  # 0-based
    pos = cur_start
    for op, length in cigar:
        if op == pysam.CREF_SKIP:
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
            write_transcript(fh, Transcript(f"{qname}.m{seen[qname]}", chrom, strand, source,
                                            exons))
            n += 1
    return n


# ---------------------------------------------------------------------------
# Step entry point
# ---------------------------------------------------------------------------


def _convert_track_bams(wd: Path) -> None:
    """Convert each Trinity transcript->genome BAM to ``<track>.genome.gff3`` models.

    Runs whether or not correction does: the SL cut consumes these (the conversion
    used to live in run_sl_cut).
    """
    for _, bam_name in TRANSCRIPT_SETS:
        bam = wd / bam_name
        if not bam.exists():
            continue
        out = wd / f"{bam.stem}{RAW_GFF3}"  # <track>.genome.gff3
        n = bam_to_transcript_gff3(bam, out, source=bam_name[: -len(GENOME_BAM_SUFFIX)])
        log.info("Converted %s -> %s (%d models).", bam.name, out.name, n)


def _blastx_coding_strands(
    config: AssemblyConfig, models_by_tag: dict[str, list[Transcript]], contigs: ContigIndex
) -> dict[str, str]:
    """Both-strand blastx -> per-hit query frame -> the coding strand of each
    confirmed transcript, keyed ``{tag}:{tid}``.

    A positive frame means the label is coding, a negative one that its reverse
    complement is, so a hit both confirms and flips labels.
    """
    wd = config.work_dir
    run_blastx(
        config, models_by_tag, contigs, wd / _QUERY_FASTA, wd / _HITS_TSV,
        evalue=config.strand_blastx_evalue, max_target_seqs=1,
        fields=("qseqid", "sseqid", "bitscore", "qframe"),
    )
    frames = parse_hits(wd / _HITS_TSV)
    coding: dict[str, str] = {}
    for tag, models in models_by_tag.items():
        for tx in models:
            qid = f"{tag}:{tx.tid}"
            if qid in frames:
                coding[qid] = coding_strand(tx.strand, frames[qid])
    return coding


def _learn_consensus(
    models_by_tag: dict[str, list[Transcript]],
    coding: dict[str, str],
    contigs: ContigIndex,
    min_consensus: int,
) -> tuple[str, str]:
    """The dominant splice consensus and its RC-twin, learned from the confirmed
    introns read on the homology-determined coding strand (not the original label)."""
    tally: Counter[str] = Counter()
    for tag, models in models_by_tag.items():
        for tx in models:
            cstrand = coding.get(f"{tag}:{tx.tid}")
            if cstrand is None:
                continue
            for istart, iend in introns_of(tx.exons):
                motif = consensus_on_strand(
                    dinucleotide(contigs, tx.chrom, istart, iend), cstrand
                )
                if motif:
                    tally[motif] += 1
    total = sum(tally.values())
    dominant, rc = _pick_consensus(tally, min_consensus)
    if total < min_consensus:
        log.warning(
            "Only %d confirmed introns (< %d); falling back to %s consensus.",
            total, min_consensus, dominant,
        )
    log.info(
        "Strand consensus %s (rc-twin %s); %d confirmed transcripts, %d introns.",
        dominant, rc, len(coding), total,
    )
    return dominant, rc


def _correct_and_write(
    wd: Path,
    models_by_tag: dict[str, list[Transcript]],
    coding: dict[str, str],
    dominant: str,
    rc: str,
    contigs: ContigIndex,
) -> Counter[str]:
    """Re-strand every transcript; write each track's ``<stem>.stranded.gff3`` and the
    audit TSV. Returns how often each decision was made."""
    counts: Counter[str] = Counter()
    with open(wd / _AUDIT_TSV, "w") as tsv:
        tsv.write("set\ttid\tn_introns\told_strand\tnew_strand\tdecision\thit\n")
        for tag, models in models_by_tag.items():
            for tx in models:
                hit_coding = coding.get(f"{tag}:{tx.tid}")
                n_introns = len(introns_of(tx.exons))
                new, decision = _decide(tx, hit_coding, dominant, rc, contigs)
                counts[decision] += 1
                old = tx.strand
                tx.strand = new
                tsv.write(
                    f"{tag}\t{tx.tid}\t{n_introns}\t{old}\t{new}\t{decision}"
                    f"\t{int(hit_coding is not None)}\n"
                )
            write_transcripts(models, wd / f"{tag}{STRANDED_GFF3}")
    return counts


def run_strand_correction(config: AssemblyConfig) -> None:
    """Convert the Trinity transcript->genome BAMs to models and (when enabled)
    homology-correct strand."""
    wd = config.work_dir
    _convert_track_bams(wd)

    # Clear any stranded models from a prior run before deciding whether to rewrite
    # them. When correction is now a no-op (stranded library, no --uniprot, or a
    # re-run that re-clipped a track), a stale *.stranded.gff3 would otherwise
    # shadow the fresh resolve_model_source() fallback in sl_cut and keep the
    # de-fused models out of combinr. The active path rewrites them below.
    for stem in mapped_transcript_stems():
        (wd / f"{stem}{STRANDED_GFF3}").unlink(missing_ok=True)

    # Gate: only for unstranded libraries with a protein DB supplied.
    if config.strand_specific is not None:
        log.info(
            "Strand-specific library (%s); skipping strand correction.",
            config.strand_specific,
        )
        return
    if config.uniprot_db is None:
        log.info("No --uniprot DB; skipping homology-based strand correction.")
        return

    # The raw genome GFF3 of each Trinity track present (from _convert_track_bams),
    # corrected into a per-track ``<stem>.stranded.gff3``.
    models_by_tag = {
        stem: read_transcripts(wd / f"{stem}{RAW_GFF3}")
        for stem in mapped_transcript_stems()
        if (wd / f"{stem}{RAW_GFF3}").exists()
    }
    if not models_by_tag:
        log.warning("No transcript models to strand-correct.")
        return

    with ContigIndex(config.genome) as contigs:
        coding = _blastx_coding_strands(config, models_by_tag, contigs)
        dominant, rc = _learn_consensus(
            models_by_tag, coding, contigs, config.min_strand_consensus
        )
        counts = _correct_and_write(wd, models_by_tag, coding, dominant, rc, contigs)

    log.info(
        "Strand correction: %d flipped, %d assigned, %d kept, %d mono-exon, "
        "%d ambiguous.",
        counts["flip"], counts["assign"], counts["keep"],
        counts["mono-exon"], counts["ambiguous"],
    )
