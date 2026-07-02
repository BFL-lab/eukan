"""Homology-grounded transcript de-fusion: split chimeric transcripts.

In gene-dense / polycistronic genomes a noisy assembler can fuse two adjacent genes
into one transcript. Overlap arithmetic can't tell a fusion from a real long gene,
but protein homology
can: an ``--ultra-sensitive`` ``diamond blastx`` of each transcript against SwissProt
that finds **>=2 distinct subjects on non-overlapping query ranges** is direct
evidence of a chimera, and the gap between the hit regions tells us where to cut.

This is :mod:`eukan.assembly.jaccard` with a *homology* clip signal instead of a
*read-coverage* one: it reuses jaccard's exon-partitioning split
(:func:`jaccard._partition_exons`) and writer, feeding it spliced cut offsets derived
from the inter-hit gaps rather than read-pair troughs. Each split piece takes the
coding strand implied by the hit inside it (:func:`strand_correction._coding_strand`),
so an opposite-strand fusion is re-oriented correctly.

Runs **after** ``strand_correct`` (so it splits the homology-stranded models) and
**before** ``sl_cut`` (which prefers the ``*.defuse.gff3`` outputs). A no-op unless
``--defuse`` and ``--uniprot`` are both given; it then always writes a
``<track>.genome.defuse.gff3`` per Trinity track (copy-through when nothing is split)
plus a ``defuse.tsv`` audit table. Read coverage at each seam is logged as
**advisory** only (never gates a split).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import pysam

from eukan.assembly.jaccard import (
    _parse_transcript_models,
    _partition_exons,
    _Tx,
    _write_transcript_models_gff3,
)
from eukan.assembly.strand_correction import (
    _DIAMOND_BLOCK_SIZE,
    _DIAMOND_INDEX_CHUNKS,
    _coding_strand,
    _resolve_diamond_db,
    _stitch,
)
from eukan.assembly.tracks import mapped_transcript_stems
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_QUERY_FASTA = "defuse_query.fasta"
_HITS_TSV = "defuse_blastx.tsv"
_AUDIT_TSV = "defuse.tsv"


@dataclass
class _Hit:
    """One blastx HSP in query (spliced, 5'->3') coordinates."""

    sseqid: str
    qlo: int  # 1-based inclusive, qlo <= qhi
    qhi: int
    bitscore: float
    qframe: int


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def parse_positional_hits(path: Path) -> dict[str, list[_Hit]]:
    """Group blastx HSPs by query id.

    Tabular columns are ``qseqid sseqid qstart qend bitscore qframe`` (diamond
    ``--outfmt 6``). Query coordinates are normalised so ``qlo <= qhi`` (a minus-frame
    hit reports ``qstart > qend``); the frame sign is preserved for per-piece strand.
    """
    out: dict[str, list[_Hit]] = {}
    with path.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 6 or row[0].startswith("#"):
                continue
            try:
                qstart, qend, bitscore, qframe = (
                    int(row[2]), int(row[3]), float(row[4]), int(row[5])
                )
            except ValueError:
                continue
            lo, hi = (qstart, qend) if qstart <= qend else (qend, qstart)
            out.setdefault(row[0], []).append(_Hit(row[1], lo, hi, bitscore, qframe))
    return out


def _overlap_frac(a: _Hit, b: _Hit) -> float:
    """Query overlap of two hits as a fraction of the shorter hit's length."""
    ov = min(a.qhi, b.qhi) - max(a.qlo, b.qlo) + 1
    if ov <= 0:
        return 0.0
    return ov / min(a.qhi - a.qlo + 1, b.qhi - b.qlo + 1)


def distinct_nonoverlapping(hits: list[_Hit], tolerance: float) -> list[_Hit] | None:
    """The distinct-subject hits supporting a fusion, ordered along the transcript.

    Keeps the best-bitscore HSP per subject, then greedily selects (highest bitscore
    first) those that overlap every already-selected hit by at most *tolerance* of the
    shorter — i.e. *distinct, non-overlapping* protein evidence. Returns the selected
    set sorted by query position when >=2 survive (a chimera), else ``None``.
    """
    best_by_subj: dict[str, _Hit] = {}
    for h in hits:
        cur = best_by_subj.get(h.sseqid)
        if cur is None or h.bitscore > cur.bitscore:
            best_by_subj[h.sseqid] = h
    selected: list[_Hit] = []
    for h in sorted(best_by_subj.values(), key=lambda x: x.bitscore, reverse=True):
        if all(_overlap_frac(h, s) <= tolerance for s in selected):
            selected.append(h)
    if len(selected) < 2:
        return None
    selected.sort(key=lambda h: (h.qlo, h.qhi))
    return selected


def _piece_strand(label: str, chosen: list[_Hit], lo: int, hi: int) -> str:
    """Coding strand for a split piece spanning spliced ``[lo, hi]``.

    Uses the frame of the hit that overlaps the piece most (so each gene's piece is
    oriented by its own protein); falls back to the parent label when no hit overlaps.
    """
    best: _Hit | None = None
    best_ov = 0
    for h in chosen:
        ov = min(hi, h.qhi) - max(lo, h.qlo) + 1
        if ov > best_ov:
            best_ov, best = ov, h
    return _coding_strand(label, best.qframe) if best is not None else label


def split_fused(tx: _Tx, chosen: list[_Hit], min_segment: int) -> list[_Tx] | None:
    """Split *tx* between consecutive distinct-subject hit regions.

    Each gap midpoint ``(prev.qhi + next.qlo)//2`` is a spliced "cut after base P"
    offset (the same space jaccard's read clips use). Reuses
    :func:`jaccard._partition_exons` to cut the exon blocks; pieces shorter than
    *min_segment* are dropped; each piece is relabelled ``.d1/.d2/...`` and oriented by
    its own hit. ``None`` if no valid cut remains.
    """
    total_len = sum(e - s + 1 for s, e in tx.exons)
    clips = sorted(
        (chosen[k].qhi + chosen[k + 1].qlo) // 2 for k in range(len(chosen) - 1)
    )
    clips = [c for c in clips if 0 < c < total_len]
    if not clips:
        return None
    exons_5to3 = tx.exons if tx.strand != "-" else list(reversed(tx.exons))
    cut_points = [0, *clips, total_len]
    out: list[_Tx] = []
    for i, blocks in enumerate(_partition_exons(exons_5to3, clips, tx.strand, total_len)):
        if sum(e - s + 1 for s, e in blocks) < min_segment:
            continue
        strand = _piece_strand(tx.strand, chosen, cut_points[i] + 1, cut_points[i + 1])
        out.append(
            _Tx(f"{tx.tid}.d{len(out) + 1}", tx.chrom, strand, tx.source, sorted(blocks))
        )
    # Only a genuine multi-piece split counts; a single survivor (the rest dropped as
    # sub-min_segment) would just trim the transcript, so leave it unchanged.
    return out if len(out) >= 2 else None


def _seams(pieces: list[_Tx]) -> list[int]:
    """Genomic seam coordinates (gap midpoints) between consecutive split pieces."""
    ordered = sorted(pieces, key=lambda p: p.exons[0][0])
    return [
        (ordered[k].exons[-1][1] + ordered[k + 1].exons[0][0]) // 2
        for k in range(len(ordered) - 1)
    ]


def _seam_depth(bam: pysam.AlignmentFile | None, chrom: str, pos: int) -> int | None:
    """Advisory read depth at a 1-based genomic *pos*, or ``None`` if unavailable."""
    if bam is None:
        return None
    try:
        cov = bam.count_coverage(chrom, max(pos - 1, 0), pos)
    except (ValueError, OSError):
        return None
    return int(sum(col[0] for col in cov))


def _open_indexed_bam(bam_path: Path) -> pysam.AlignmentFile | None:
    """Open *bam_path* for region queries, building a ``.bai`` if missing. Best-effort."""
    if not bam_path.exists():
        return None
    try:
        if not Path(f"{bam_path}.bai").exists():
            pysam.index(str(bam_path))
        return pysam.AlignmentFile(str(bam_path), "rb")
    except (OSError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Step entry point
# ---------------------------------------------------------------------------


def _defuse_one_set(
    tag: str,
    in_gff: Path,
    out_gff: Path,
    hits_by_query: dict[str, list[_Hit]],
    config: AssemblyConfig,
    bam: pysam.AlignmentFile | None,
    tsv: TextIO,
) -> tuple[int, int]:
    """Split every fused transcript in one model set; write *out_gff*. Returns
    ``(n_models_in, n_fusions_split)``."""
    models = _parse_transcript_models(in_gff)
    out_models: list[_Tx] = []
    n_split = 0
    for tx in models:
        chosen = distinct_nonoverlapping(
            hits_by_query.get(f"{tag}:{tx.tid}", []), config.defuse_overlap_tolerance
        )
        pieces = (
            split_fused(tx, chosen, config.min_sl_fragment) if chosen else None
        )
        if chosen and pieces:
            out_models.extend(pieces)
            n_split += 1
            subjects = ",".join(h.sseqid for h in chosen)
            seams = _seams(pieces)
            depths = [str(_seam_depth(bam, tx.chrom, s) if bam else "NA") for s in seams]
            tsv.write(
                f"{tag}\t{tx.tid}\t{len(pieces)}\t{subjects}\t"
                f"{','.join(map(str, seams))}\t{','.join(depths)}\n"
            )
        else:
            out_models.append(tx)
    _write_transcript_models_gff3(out_models, out_gff)
    return len(models), n_split


def run_defuse(config: AssemblyConfig) -> None:
    """Split chimeric transcripts flagged by >=2 distinct non-overlapping protein hits.

    No-op (skipped by the pipeline) unless ``--defuse`` and ``--uniprot`` are set.
    Each Trinity track is rewritten to ``<track>.genome.defuse.gff3`` (copy-through
    when nothing is split).
    """
    wd = config.work_dir

    # Inputs per track: prefer the homology-stranded models from strand_correct,
    # else the raw genome GFF3 (do NOT read our own .defuse.gff3 output).
    sets: list[tuple[str, Path, Path]] = []
    for stem in mapped_transcript_stems():
        src = wd / f"{stem}.stranded.gff3"
        if not src.exists():
            src = wd / f"{stem}.gff3"
        if src.exists():
            sets.append((stem, src, wd / f"{stem}.defuse.gff3"))
    if not sets:
        log.warning("No transcript models to de-fuse.")
        return

    models_by_tag = {tag: _parse_transcript_models(inp) for tag, inp, _ in sets}
    db = _resolve_diamond_db(config)

    query = wd / _QUERY_FASTA
    with ContigIndex(config.genome) as contigs, open(query, "w") as fh:
        for tag, models in models_by_tag.items():
            for tx in sorted(models, key=lambda t: (t.chrom, t.exons[0][0])):
                fh.write(f">{tag}:{tx.tid}\n{_stitch(tx, contigs)}\n")

    run_cmd(
        ["diamond", "blastx",
         "--db", db,
         "--query", str(query),
         "--out", str(wd / _HITS_TSV),
         "--ultra-sensitive",
         "--strand", "both",
         "--query-gencode", str(config.genetic_code_obj.ncbi_id),
         "--evalue", f"{config.defuse_blastx_evalue:g}",
         "--max-target-seqs", "25",
         "--outfmt", "6", "qseqid", "sseqid", "qstart", "qend", "bitscore", "qframe",
         "--block-size", _DIAMOND_BLOCK_SIZE,
         "--index-chunks", _DIAMOND_INDEX_CHUNKS,
         "--threads", str(config.num_cpu),
         "--quiet"],
        cwd=wd,
    )

    hits_by_query = parse_positional_hits(wd / _HITS_TSV)
    bam = _open_indexed_bam(wd / config.aligner_bam)
    try:
        total_split = 0
        with open(wd / _AUDIT_TSV, "w") as tsv:
            tsv.write("set\ttid\tn_pieces\tsubjects\tseam_genomic\tseam_depth\n")
            for tag, inp, out in sets:
                n_in, n_split = _defuse_one_set(
                    tag, inp, out, hits_by_query, config, bam, tsv
                )
                total_split += n_split
                log.info(
                    "De-fuse %s -> %s: %d models (%d fused transcripts split).",
                    inp.name, out.name, n_in, n_split,
                )
    finally:
        if bam is not None:
            bam.close()

    log.info("De-fusion split %d chimeric transcript(s) total.", total_split)
