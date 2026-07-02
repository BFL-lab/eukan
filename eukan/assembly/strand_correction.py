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
import os
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

from eukan.assembly.bam_diagnostic import _dinucleotide, _reverse_complement
from eukan.assembly.jaccard import (
    _parse_transcript_models,
    _Tx,
    _write_transcript_models_gff3,
)
from eukan.assembly.sl_cut import _DENOVO_BAMS, _GENOME_BAM_SUFFIX, bam_to_transcript_gff3
from eukan.assembly.tracks import mapped_transcript_stems
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_QUERY_FASTA = "strand_query.fasta"
_HITS_TSV = "strand_blastx.tsv"
_AUDIT_TSV = "strand_correction.tsv"

# diamond memory caps for the 15 GB / 0-swap box (block-size in billions of query
# letters; index-chunks splits the reference load).
_DIAMOND_BLOCK_SIZE = "1.0"
_DIAMOND_INDEX_CHUNKS = "4"

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


def _coding_strand(label: str, qframe: int) -> str:
    """The coding strand implied by a blastx hit's query frame on a labelled query.

    The query was stitched in its labelled orientation (``.`` forward, like ``+``),
    so a positive ``qframe`` means that orientation is coding and a negative one
    means its reverse complement is. Returns ``"+"`` or ``"-"``.
    """
    base = label if label in ("+", "-") else "+"
    if qframe > 0:
        return base
    return "-" if base == "+" else "+"


def introns_of(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """0-based ``[start, end)`` introns from 1-based inclusive exon blocks.

    Matches :func:`bam_diagnostic._dinucleotide`'s convention: the first intron
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
    return f"{_reverse_complement(acceptor)}-{_reverse_complement(donor)}"


def _pick_consensus(tally: Counter[str], min_consensus: int) -> tuple[str, str]:
    """Return ``(dominant, rc_twin)`` coding-strand pair, falling back to GT-AG.

    Below *min_consensus* confirmed introns (or none at all) the calibration is too
    thin to trust, so the canonical ``GT-AG`` is used instead.
    """
    total = sum(tally.values())
    dominant = tally.most_common(1)[0][0] if total and total >= min_consensus else _CANONICAL
    return dominant, _rc_swap(dominant)


def _stitch(tx: _Tx, contigs: ContigIndex) -> str:
    """Spliced transcript sequence in the transcript's labelled orientation."""
    seq = "".join(str(contigs[tx.chrom][s - 1 : e].seq) for s, e in tx.exons)
    return str(Seq(seq).reverse_complement()) if tx.strand == "-" else seq


def _decide(
    tx: _Tx, hit_coding: str | None, dominant: str, rc: str, contigs: ContigIndex
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
        raw = _dinucleotide(contigs, tx.chrom, istart, iend)
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
# diamond DB resolution
# ---------------------------------------------------------------------------


def _resolve_diamond_db(config: AssemblyConfig) -> str:
    """Return the diamond ``--db`` base (no ``.dmnd``); build it from a FASTA if needed.

    A ``.dmnd`` path is used directly. A FASTA is compiled with ``diamond makedb``,
    cached as ``<stem>.dmnd`` beside the FASTA when that dir is writable, else in the
    work dir, and reused while it is at least as new as the source.
    """
    src = config.uniprot_db
    assert src is not None  # gated by the caller
    if src.suffix == ".dmnd":
        return str(src.with_suffix(""))  # diamond --db wants the base name

    base = src.with_suffix("")
    out = Path(f"{base}.dmnd")
    if not (src.parent.exists() and os.access(src.parent, os.W_OK)):
        base = config.work_dir / "uniprot_sprot"
        out = Path(f"{base}.dmnd")
    if not (out.exists() and out.stat().st_mtime >= src.stat().st_mtime):
        log.info("Building diamond DB %s from %s ...", out.name, src.name)
        run_cmd(
            ["diamond", "makedb", "--in", str(src), "--db", str(base),
             "--threads", str(config.num_cpu), "--quiet"],
            cwd=config.work_dir,
        )
    return str(base)


# ---------------------------------------------------------------------------
# Step entry point
# ---------------------------------------------------------------------------


def run_strand_correction(config: AssemblyConfig) -> None:
    """Convert the Trinity transcript->genome BAMs to models and (when enabled)
    homology-correct strand."""
    wd = config.work_dir

    # 1. Always: each Trinity transcript->genome BAM -> gene>mRNA>exon GFF3 (the SL
    #    cut consumes these, replacing the conversion that used to live in run_sl_cut).
    for bam_name in _DENOVO_BAMS:
        bam = wd / bam_name
        if not bam.exists():
            continue
        stem = bam_name[: -len(_GENOME_BAM_SUFFIX)]
        n = bam_to_transcript_gff3(bam, wd / f"{stem}.genome.gff3", source=stem)
        log.info("Converted %s -> %s.genome.gff3 (%d models).", bam.name, stem, n)

    # Clear any stranded models from a prior run before deciding whether to rewrite
    # them. When correction is now a no-op (stranded library, no --uniprot, or a
    # re-run that re-clipped a track), a stale *.stranded.gff3 would otherwise
    # shadow the fresh resolve_model_source() fallback in sl_cut and keep the
    # de-fused models out of combinr. The active path rewrites them below.
    for stem in mapped_transcript_stems():
        (wd / f"{stem}.stranded.gff3").unlink(missing_ok=True)

    # 2. Gate: only for unstranded libraries with a protein DB supplied.
    if config.strand_specific is not None:
        log.info(
            "Strand-specific library (%s); skipping strand correction.",
            config.strand_specific,
        )
        return
    if config.uniprot_db is None:
        log.info("No --uniprot DB; skipping homology-based strand correction.")
        return

    # 3. The model sets present — the raw genome GFF3 of each Trinity track (built
    #    in step 1), correcting into a per-track ``<stem>.stranded.gff3``.
    sets: list[tuple[str, Path, Path]] = []
    for stem in mapped_transcript_stems():
        raw = wd / f"{stem}.gff3"
        if raw.exists():
            sets.append((stem, raw, wd / f"{stem}.stranded.gff3"))
    if not sets:
        log.warning("No transcript models to strand-correct.")
        return
    models_by_tag = {tag: _parse_transcript_models(inp) for tag, inp, _ in sets}

    db = _resolve_diamond_db(config)

    with ContigIndex(config.genome) as contigs:
        # 4. Both-strand blastx -> per-hit query frame -> the coding strand of each
        #    confirmed transcript (positive frame = label is coding, negative = its
        #    reverse complement is). This both confirms and flips labels.
        query = wd / _QUERY_FASTA
        with open(query, "w") as fh:
            for tag, models in models_by_tag.items():
                for tx in sorted(models, key=lambda t: (t.chrom, t.exons[0][0])):
                    fh.write(f">{tag}:{tx.tid}\n{_stitch(tx, contigs)}\n")
        run_cmd(
            ["diamond", "blastx",
             "--db", db,
             "--query", str(query),
             "--out", str(wd / _HITS_TSV),
             "--strand", "both",
             "--query-gencode", str(config.genetic_code_obj.ncbi_id),
             "--evalue", f"{config.strand_blastx_evalue:g}",
             "--max-target-seqs", "1",
             "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qframe",
             "--block-size", _DIAMOND_BLOCK_SIZE,
             "--index-chunks", _DIAMOND_INDEX_CHUNKS,
             "--threads", str(config.num_cpu),
             "--quiet"],
            cwd=wd,
        )
        frames = parse_hits(wd / _HITS_TSV)
        coding: dict[str, str] = {}
        for tag, models in models_by_tag.items():
            for tx in models:
                qid = f"{tag}:{tx.tid}"
                if qid in frames:
                    coding[qid] = _coding_strand(tx.strand, frames[qid])

        # 5. Learn the dominant splice consensus from confirmed introns, read on the
        #    homology-determined coding strand (not the original label).
        tally: Counter[str] = Counter()
        for tag, models in models_by_tag.items():
            for tx in models:
                cstrand = coding.get(f"{tag}:{tx.tid}")
                if cstrand is None:
                    continue
                for istart, iend in introns_of(tx.exons):
                    motif = consensus_on_strand(
                        _dinucleotide(contigs, tx.chrom, istart, iend), cstrand
                    )
                    if motif:
                        tally[motif] += 1
        total = sum(tally.values())
        dominant, rc = _pick_consensus(tally, config.min_strand_consensus)
        if total < config.min_strand_consensus:
            log.warning(
                "Only %d confirmed introns (< %d); falling back to %s consensus.",
                total, config.min_strand_consensus, dominant,
            )
        log.info(
            "Strand consensus %s (rc-twin %s); %d confirmed transcripts, %d introns.",
            dominant, rc, len(coding), total,
        )

        # 6. Correct each transcript; write stranded models + the audit TSV.
        counts: Counter[str] = Counter()
        with open(wd / _AUDIT_TSV, "w") as tsv:
            tsv.write("set\ttid\tn_introns\told_strand\tnew_strand\tdecision\thit\n")
            for tag, _inp, out in sets:
                models = models_by_tag[tag]
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
                _write_transcript_models_gff3(models, out)

    log.info(
        "Strand correction: %d flipped, %d assigned, %d kept, %d mono-exon, "
        "%d ambiguous.",
        counts["flip"], counts["assign"], counts["keep"],
        counts["mono-exon"], counts["ambiguous"],
    )
