"""Jaccard clipping of fused transcripts — a standalone reimplementation of
Trinity's ``--jaccard_clip`` post-processing.

In gene-dense genomes with overlapping UTRs an assembler can string two
independently-transcribed genes into a single contig. Trinity's
``--jaccard_clip`` detects these by mapping read pairs back to the assembled
contigs and cutting where no read-pair fragment *bridges* a position; but it
only ever clips Trinity's own output. This module applies the same logic
uniformly to the de novo assembled transcripts (rnaSPAdes), so the
consolidation step (combinr) never sees an un-clipped fused contig.

Algorithm (faithful to Trinity, see ``trinityrnaseq/util/support_scripts/``):

1. Map paired reads to the transcript FASTA, ungapped (STAR ``EndToEnd`` with
   introns disabled). One transcript == one BAM reference.
2. Extract proper-pair fragment spans ``(lend, rend)`` (insert 100-500 bp).
3. Per-position jaccard over a sliding window ``W`` (two adjacent W-wide
   windows): ``(n_both + 1) / (n_single + n_both + 1)`` where ``n_both`` counts
   fragments spanning the whole window (bridging the junction) and ``n_single``
   counts fragments touching exactly one window edge. A fusion junction shows a
   trough: read pairs from the two genes don't bridge it, so ``n_both`` drops to
   ~0 while ``n_single`` stays high.
4. Call troughs (``jaccard <= 0.05``) flanked by "hills" (``jaccard >= trough +
   0.35``) within ``trough_win`` on both sides; reposition each clip to the
   local minimum-coverage position.
5. Split the contig at the clip points (segments shorter than ``min_segment``
   dropped).

The genome-anchored StringTie GTF is clipped the same way (:func:`_clip_stringtie_gtf`):
reads are mapped to StringTie's spliced transcripts, the same coverage troughs are
found, and the transcript>exon models are split at them (:func:`_split_models_at_clips`)
into ``stringtie.jaccard.gff3`` — so a StringTie transcript fusing two adjacent loci is
broken just as the de novo path keeps them separate. ``strand_correct``/``sl_cut`` then
prefer that clipped GFF3 over the raw GTF (:func:`resolve_stringtie_models`).

Jaccard clipping needs read PAIRS; with single-end reads the step is a no-op.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypedDict

import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from eukan.assembly.minimap2 import map_reads_to_transcripts
from eukan.gff import create_gff_db
from eukan.gff.io import iter_assembled_sequences
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# Trinity defaults (trinityrnaseq/util/support_scripts/). These are well-tuned
# and rarely retuned, so they are module constants rather than config knobs.
_WINDOW = 100              # jaccard sliding-window width (W)
_PSEUDO = 1                # jaccard pseudocount
_TROUGH_WIN = 200          # trough-scan window for clip detection
_MAX_TROUGH_VAL = 0.05     # a trough must dip to <= this jaccard
_MIN_JACCARD_DELTA = 0.35  # flanking "hills" must rise this far above the trough
_MIN_INSERT = 100          # proper-pair insert-size floor
_MAX_INSERT = 500          # proper-pair insert-size ceiling
_MIN_SEGMENT = 25          # shortest split segment kept (Trinity's k-mer floor)
_REPOSITION_HALF_WIN = _TROUGH_WIN // 2  # coverage-reposition search radius (100)
# Ceiling for the coverage-adapted trough gate: below the depth where a clean
# junction would need a jaccard above this, the bridging signal is too weak to
# tell a fusion from noise, so the gate stops relaxing (never over-clips a real
# single transcript that simply has a thin patch).
_MAX_ADAPTIVE_TROUGH = 0.30

# The Trinity FASTAs the map_transcripts step consumes (both de novo and
# genome-guided); the jaccard step rewrites each into a ``.jaccard.fasta`` sibling
# that ``star.map_transcripts`` prefers when present.
_TRANSCRIPT_FASTAS = (
    "trinity-denovo.fasta",
    "trinity-gg.fasta",
)

# DORMANT (StringTie path): the genome-guided StringTie GTF and the clipped sibling
# the jaccard step writes for it. Not produced by the active Trinity pipeline, so
# the run_jaccard StringTie branch and resolve_stringtie_models below are no-ops
# unless the kept-but-unwired stringtie module is re-activated.
_STRINGTIE_GTF = "stringtie.gtf"
STRINGTIE_JACCARD_GFF3 = "stringtie.jaccard.gff3"


def jaccard_output_name(fasta_name: str) -> str:
    """The ``.jaccard.fasta`` name the jaccard step writes for *fasta_name*."""
    return fasta_name.replace(".fasta", ".jaccard.fasta")


# ---------------------------------------------------------------------------
# Core metric: per-position jaccard and coverage (pure, unit-testable)
# ---------------------------------------------------------------------------


def jaccard_array(
    frags: list[tuple[int, int]],
    length: int,
    *,
    window: int = _WINDOW,
    pseudo: int = _PSEUDO,
) -> list[float]:
    """Per-position jaccard for one transcript of length *length*.

    *frags* are 1-based inclusive proper-pair fragment spans. Returns a list
    indexed 1..length (index 0 is an unused 0.0), where ``jac[mid]`` is the
    jaccard for the window whose left edge is ``mid - (window-1)//2`` — matching
    Trinity's ``mid``-indexed wig. Computed in one O(length + len(frags)) sweep
    via difference arrays:

    * ``touch[x]`` = number of fragments covering coordinate ``x``;
    * ``both[j]``  = number of fragments spanning the whole window ``[j, j+W-1]``
      (i.e. ``lend <= j`` and ``rend >= j+W-1``), which holds for ``j`` in
      ``[lend, rend-W+1]``.

    Then ``n_single(j) = (touch[j] - both[j]) + (touch[j+W-1] - both[j])`` and
    ``jaccard(j) = (both[j] + pseudo) / (n_single + both[j] + pseudo)``.
    """
    L = length
    W = window
    if L <= 0:
        return [0.0]
    size = L + W + 2
    touch = [0] * (size + 2)
    both = [0] * (size + 2)
    for lend, rend in frags:
        a = max(1, lend)
        b = min(L, rend)
        if a > b:
            continue
        touch[a] += 1
        touch[b + 1] -= 1
        # spans the full window [j, j+W-1]  <=>  a <= j  and  j+W-1 <= b
        #   => j in [a, b - W + 1]
        hi = b - W + 1
        if hi >= a:
            both[a] += 1
            both[hi + 1] -= 1
    for i in range(1, size + 1):
        touch[i] += touch[i - 1]
        both[i] += both[i - 1]

    jac = [0.0] * (L + 1)
    half = (W - 1) // 2
    for j in range(1, L + 1):
        wr = j + W - 1
        nb = both[j]
        t_left = touch[j]
        t_right = touch[wr] if wr <= size else 0
        n_single = (t_left - nb) + (t_right - nb)
        val = round((nb + pseudo) / (n_single + nb + pseudo), 4)
        mid = j + half
        if 1 <= mid <= L:
            jac[mid] = val
    return jac


def coverage_array(frags: list[tuple[int, int]], length: int) -> list[int]:
    """Per-position fragment coverage, indexed 1..length (index 0 is 0)."""
    L = length
    cov = [0] * (L + 2)
    for lend, rend in frags:
        a = max(1, lend)
        b = min(L, rend)
        if a > b:
            continue
        cov[a] += 1
        cov[b + 1] -= 1
    for i in range(1, L + 1):
        cov[i] += cov[i - 1]
    return cov[: L + 1]


# ---------------------------------------------------------------------------
# Trough -> clip detection (pure, unit-testable)
# ---------------------------------------------------------------------------


@dataclass
class Trough:
    """A candidate clip point: position, its jaccard, and dip depth."""

    pos: int
    jaccard: float
    avg_delta: float


def _candidate_troughs(
    jac: list[float],
    trough_win: int,
    max_trough: float,
    *,
    cov: list[int] | None = None,
    greed: float = 0.0,
    max_adaptive_trough: float = _MAX_ADAPTIVE_TROUGH,
) -> list[Trough]:
    """Positions whose jaccard dips to <= the (optionally coverage-adapted) trough.

    With ``greed > 0`` and a coverage array, the absolute *max_trough* floor is
    widened per position toward the deepest jaccard a clean (``n_both ≈ 0``)
    junction can physically reach at the local read-pair depth —
    ``pseudo / (left + right + pseudo)`` from the two flanking hill coverages —
    times *greed*, capped at *max_adaptive_trough*. At high depth that
    product sits below *max_trough*, so the strict floor stands; at low depth
    (where the pseudocount alone keeps a real junction above 0.05) the gate relaxes
    so the junction becomes a candidate. The flanking-hill and coverage-reposition
    tests downstream still decide whether a candidate is a true fusion, so relaxing
    the gate adds sensitivity without by itself widening what is finally cut.
    """
    L = len(jac) - 1
    half = trough_win // 2
    troughs: list[Trough] = []
    for i in range(trough_win, L + 1):
        center = i - half
        c = jac[center]
        thr = max_trough
        if cov is not None and greed > 0.0:
            # Hill coverage on both sides; clamp the left index off the unused
            # cov[0] sentinel at the 5' start so the gate isn't lopsidedly relaxed
            # at the very first scanned position.
            flank = cov[max(1, i - trough_win)] + cov[i]
            achievable = _PSEUDO / (flank + _PSEUDO)
            thr = max(max_trough, min(greed * achievable, max_adaptive_trough))
        if c <= thr:
            left = jac[i - trough_win]
            right = jac[i]
            avg_delta = ((left - c) + (right - c)) / 2
            troughs.append(Trough(center, c, avg_delta))
    return troughs


def _group_and_pick_best(troughs: list[Trough], trough_win: int) -> list[Trough]:
    """Group troughs within *trough_win* of each other; keep the deepest each.

    Candidates arrive in ascending position. Grouping is chained (gap measured
    to the last group member). Best = lowest jaccard, tie-break highest delta.
    """
    if not troughs:
        return []
    groups: list[list[Trough]] = [[troughs[0]]]
    for t in troughs[1:]:
        if t.pos - groups[-1][-1].pos <= trough_win:
            groups[-1].append(t)
        else:
            groups.append([t])
    return [min(g, key=lambda t: (t.jaccard, -t.avg_delta)) for g in groups]


def _require_hills(
    clips: list[Trough], jac: list[float], trough_win: int, min_delta: float
) -> list[Trough]:
    """Keep clips flanked by a jaccard "hill" within *trough_win* on both sides."""
    L = len(jac) - 1
    validated: list[Trough] = []
    for clip in clips:
        hill_min = clip.jaccard + min_delta
        left_lo = max(1, clip.pos - trough_win)
        right_hi = min(L, clip.pos + trough_win)
        left_hill = any(jac[i] >= hill_min for i in range(left_lo, clip.pos))
        right_hill = any(jac[i] >= hill_min for i in range(clip.pos + 1, right_hi + 1))
        if left_hill and right_hill:
            validated.append(clip)
    return validated


def _reposition_by_coverage(clips: list[Trough], cov: list[int], half_win: int) -> None:
    """Move each clip to the minimum-coverage position within +/- *half_win*."""
    L = len(cov) - 1
    for clip in clips:
        lo = max(1, clip.pos - half_win)
        hi = min(L, clip.pos + half_win)
        min_cov = cov[clip.pos]
        min_pos = clip.pos
        for i in range(lo, hi + 1):
            if cov[i] < min_cov:
                min_cov = cov[i]
                min_pos = i
        clip.pos = min_pos


def find_clip_points(
    frags: list[tuple[int, int]],
    length: int,
    *,
    window: int = _WINDOW,
    trough_win: int = _TROUGH_WIN,
    max_trough: float = _MAX_TROUGH_VAL,
    min_delta: float = _MIN_JACCARD_DELTA,
    pseudo: int = _PSEUDO,
    greed: float = 0.0,
    max_adaptive_trough: float = _MAX_ADAPTIVE_TROUGH,
) -> list[int]:
    """Clip positions (sorted, 1-based) for one transcript's fragment spans.

    *greed* (>0) makes the trough gate coverage-adaptive so low-coverage fusions
    are split too; 0 keeps Trinity's fixed *max_trough* floor (see
    :func:`_candidate_troughs`). *max_adaptive_trough* caps how far that
    adaptation may relax the floor at low depth.
    """
    jac = jaccard_array(frags, length, window=window, pseudo=pseudo)
    cov = coverage_array(frags, length)
    troughs = _candidate_troughs(
        jac, trough_win, max_trough, cov=cov, greed=greed,
        max_adaptive_trough=max_adaptive_trough,
    )
    if not troughs:
        return []
    clips = _require_hills(_group_and_pick_best(troughs, trough_win), jac, trough_win, min_delta)
    if not clips:
        return []
    _reposition_by_coverage(clips, cov, _REPOSITION_HALF_WIN)
    return sorted(c.pos for c in clips)


def split_fasta_record(seq: str, clips: list[int], min_segment: int) -> list[str]:
    """Cut *seq* at *clips* (1-based, cut after that base); drop short pieces."""
    segments: list[str] = []
    start = 1
    for clip in sorted(clips):
        piece = seq[start - 1 : clip]
        if len(piece) >= min_segment:
            segments.append(piece)
        start = clip + 1
    tail = seq[start - 1 :]
    if len(tail) >= min_segment:
        segments.append(tail)
    return segments


# ---------------------------------------------------------------------------
# Read mapping (STAR EndToEnd, ungapped) + fragment extraction
# ---------------------------------------------------------------------------


def _genome_stats(fasta: Path) -> tuple[int, int]:
    """Total sequence length and sequence count in *fasta* (one pass)."""
    total = n = 0
    for rec in SeqIO.parse(str(fasta), "fasta"):
        total += len(rec.seq)
        n += 1
    return total, n


def iter_fragment_spans(
    bam_path: Path, *, min_insert: int = _MIN_INSERT, max_insert: int = _MAX_INSERT
) -> Iterator[tuple[str, list[tuple[int, int]]]]:
    """Yield ``(reference, fragment_spans)`` per BAM reference.

    A fragment span is ``(lend, rend)`` (1-based inclusive) covering a proper
    read pair, kept only when ``min_insert <= insert <= max_insert``. Mates are
    paired by query name within each reference (so memory is bounded by one
    transcript's read depth); secondary/supplementary/improper reads are skipped.
    """
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for ref in bam.references:
            spans: list[tuple[int, int]] = []
            mates: dict[str, pysam.AlignedSegment] = {}
            for read in bam.fetch(reference=ref):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.is_proper_pair:
                    continue
                qname, start, end = read.query_name, read.reference_start, read.reference_end
                if qname is None or start is None or end is None:
                    continue
                mate = mates.pop(qname, None)
                if mate is None:
                    mates[qname] = read
                    continue
                m_start, m_end = mate.reference_start, mate.reference_end
                if m_start is None or m_end is None:
                    continue
                lend = min(start, m_start) + 1
                rend = max(end, m_end)
                insert = rend - lend + 1
                if min_insert <= insert <= max_insert:
                    spans.append((lend, rend))
            yield ref, spans
    finally:
        bam.close()


class _ClipKnobs(TypedDict):
    """The tunable :func:`find_clip_points` detection knobs (from config)."""

    max_trough: float
    min_delta: float
    greed: float
    max_adaptive_trough: float


def _clip_knobs(config: AssemblyConfig) -> _ClipKnobs:
    """Bundle the config's jaccard detection knobs as ``find_clip_points`` kwargs.

    Single source for both the de novo FASTA and StringTie GTF clip paths so the
    two stay in lockstep as knobs are added.
    """
    return {
        "max_trough": config.jaccard_max_trough,
        "min_delta": config.jaccard_min_delta,
        "greed": config.jaccard_greediness,
        "max_adaptive_trough": config.jaccard_max_adaptive_trough,
    }


def _clip_one_fasta(config: AssemblyConfig, src: Path, out: Path) -> tuple[int, int, int]:
    """Jaccard-clip every record in *src* into *out*; return (n_in, n_out, n_clipped)."""
    bam = map_reads_to_transcripts(config, src, out.name.replace(".fasta", ""))

    handle = pysam.AlignmentFile(str(bam), "rb")
    ref_len = dict(zip(handle.references, handle.lengths, strict=False))
    handle.close()

    clip_map: dict[str, list[int]] = {}
    for ref, spans in iter_fragment_spans(bam):
        if not spans:
            continue
        clips = find_clip_points(spans, ref_len[ref], **_clip_knobs(config))
        if clips:
            clip_map[ref] = clips

    n_in = n_out = n_clipped = 0
    with open(out, "w") as fh:
        for rec in SeqIO.parse(str(src), "fasta"):
            n_in += 1
            seq = str(rec.seq)
            rec_clips = clip_map.get(rec.id)
            segments = (
                split_fasta_record(seq, rec_clips, config.min_sl_fragment) if rec_clips else []
            )
            if not segments:
                # No clips, or every piece fell below the length floor: keep the
                # original so a spurious end-clip never deletes a real transcript.
                fh.write(f">{rec.id}\n{seq}\n")
                n_out += 1
                continue
            n_clipped += 1
            for k, piece in enumerate(segments, start=1):
                fh.write(f">{rec.id}.j{k}\n{piece}\n")
                n_out += 1
    return n_in, n_out, n_clipped


def run_jaccard(config: AssemblyConfig) -> None:
    """Jaccard-clip the de novo transcript FASTAs and the StringTie GTF.

    No-op on single-end reads (read-pair bridging is the clip signal).
    """
    wd = config.work_dir
    if not (config.left_reads and config.right_reads):
        log.warning(
            "Jaccard clipping needs paired reads (read-pair bridging is the "
            "signal); single-end input — skipping."
        )
        return

    for name in _TRANSCRIPT_FASTAS:
        src = wd / name
        if not src.exists() or src.stat().st_size == 0:
            continue
        out = wd / jaccard_output_name(name)
        n_in, n_out, n_clipped = _clip_one_fasta(config, src, out)
        log.info(
            "Jaccard clip %s -> %s: %d -> %d sequences (%d contigs split).",
            src.name, out.name, n_in, n_out, n_clipped,
        )

    # Genome-guided StringTie models: clip the GTF the same way (read-pair troughs
    # on the spliced transcripts) so a fused StringTie locus is split too.
    gtf = wd / _STRINGTIE_GTF
    if gtf.exists() and gtf.stat().st_size > 0:
        n_in, n_split = _clip_stringtie_gtf(config)
        log.info(
            "Jaccard clip %s -> %s: %d models (%d fused models split).",
            gtf.name, STRINGTIE_JACCARD_GFF3, n_in, n_split,
        )


# ---------------------------------------------------------------------------
# GFF3 / GTF clip path (genome-anchored transcript>exon models)
# ---------------------------------------------------------------------------


def _partition_exons(
    exons: list[tuple[int, int]],
    clips: list[int],
    strand: str,
    total_len: int,
) -> list[list[tuple[int, int]]]:
    """Partition a transcript's exon blocks at spliced clip positions.

    *exons* are genomic ``(start, end)`` blocks ordered 5'->3' along the
    transcript (ascending for ``+``, descending for ``-``). Each *clip* ``P``
    cuts after spliced base ``P``. Returns one segment per inter-clip interval;
    each segment is a list of genomic ``(start, end)`` blocks. The within-exon
    split mirrors :func:`eukan.annotation.orf._map_orf_plus_strand` /
    ``_map_orf_minus_strand`` (minus strand maps high genomic -> low spliced).
    """
    cut_points = [0, *sorted(p for p in clips if 0 < p < total_len), total_len]
    spans: list[tuple[int, int, int, int]] = []  # (tx_start0, tx_end, gstart, gend)
    consumed = 0
    for gstart, gend in exons:
        elen = gend - gstart + 1
        spans.append((consumed, consumed + elen, gstart, gend))
        consumed += elen

    segments: list[list[tuple[int, int]]] = []
    for s in range(len(cut_points) - 1):
        u, v = cut_points[s], cut_points[s + 1]  # this segment is spliced bases [u+1, v]
        blocks: list[tuple[int, int]] = []
        for tx_start, tx_end, gstart, gend in spans:
            lo = max(u + 1, tx_start + 1)
            hi = min(v, tx_end)
            if lo > hi:
                continue
            g1 = lo - tx_start  # 1-based offset within exon (5'->3')
            g2 = hi - tx_start
            if strand == "-":
                blocks.append((gend - g2 + 1, gend - g1 + 1))
            else:
                blocks.append((gstart + g1 - 1, gstart + g2 - 1))
        segments.append(blocks)
    return segments


@dataclass
class _Tx:
    """A transcript model: id, location, and genomic exon blocks (start-sorted)."""

    tid: str
    chrom: str
    strand: str
    source: str
    exons: list[tuple[int, int]] = field(default_factory=list)


def _parse_attrs(col9: str) -> dict[str, str]:
    """Parse column 9 of either GFF3 (``key=value``) or GTF (``key "value"``)."""
    attrs: dict[str, str] = {}
    for part in col9.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:  # GFF3
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
        else:  # GTF: key "value"
            bits = part.split(None, 1)
            if len(bits) == 2:
                attrs[bits[0].strip()] = bits[1].strip().strip('"')
    return attrs


def _parse_transcript_models(gff: str | Path) -> list[_Tx]:
    """Parse a transcript>exon GFF3 or GTF into transcript records.

    Groups ``exon`` features by their transcript id (``Parent=`` in GFF3,
    ``transcript_id`` in GTF), tolerating files with no explicit transcript/mRNA
    lines (StringTie emits both; combinr-style GFF3 may not). Exons are sorted by
    genomic start; transcript order follows first appearance.
    """
    by_id: dict[str, _Tx] = {}
    order: list[str] = []

    def _get(tid: str, cols: list[str]) -> _Tx:
        tx = by_id.get(tid)
        if tx is None:
            tx = _Tx(tid, cols[0], cols[6], cols[1])
            by_id[tid] = tx
            order.append(tid)
        return tx

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = _parse_attrs(cols[8])
            if cols[2] in ("transcript", "mRNA"):
                tid = attrs.get("ID") or attrs.get("transcript_id")
                if tid:
                    _get(tid, cols)
            elif cols[2] == "exon":
                tid = attrs.get("Parent") or attrs.get("transcript_id")
                if tid:
                    _get(tid, cols).exons.append((int(cols[3]), int(cols[4])))

    models = []
    for tid in order:
        tx = by_id[tid]
        if tx.exons:
            tx.exons.sort()
            models.append(tx)
    return models


def _split_transcript(tx: _Tx, clips: list[int], min_segment: int) -> list[_Tx] | None:
    """Split *tx* at spliced clip positions; ``None`` if nothing valid to cut."""
    exons_5to3 = tx.exons if tx.strand != "-" else list(reversed(tx.exons))
    total_len = sum(e - s + 1 for s, e in tx.exons)
    valid = [p for p in clips if 0 < p < total_len]
    if not valid:
        return None
    out: list[_Tx] = []
    for blocks in _partition_exons(exons_5to3, valid, tx.strand, total_len):
        if sum(e - s + 1 for s, e in blocks) < min_segment:
            continue
        out.append(_Tx(f"{tx.tid}.j{len(out) + 1}", tx.chrom, tx.strand, tx.source, sorted(blocks)))
    return out or None


def _write_transcript_models_gff3(models: list[_Tx], out_gff: str | Path) -> None:
    """Write *models* as a gene>mRNA>exon GFF3 (one synthetic gene per transcript)."""
    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for tx in sorted(models, key=lambda t: (t.chrom, t.exons[0][0])):
            gid, gstart, gend = f"{tx.tid}.gene", tx.exons[0][0], tx.exons[-1][1]
            loc = f"{tx.chrom}\t{tx.source}"
            fh.write(f"{loc}\tgene\t{gstart}\t{gend}\t.\t{tx.strand}\t.\tID={gid}\n")
            fh.write(f"{loc}\tmRNA\t{gstart}\t{gend}\t.\t{tx.strand}\t.\tID={tx.tid};Parent={gid}\n")
            for k, (s, e) in enumerate(tx.exons, start=1):
                fh.write(
                    f"{loc}\texon\t{s}\t{e}\t.\t{tx.strand}\t.\tID={tx.tid}.exon{k};Parent={tx.tid}\n"
                )


def _write_tx(fh, tx: _Tx) -> None:
    """Write one transcript as gene>mRNA>exon GFF3 (combinr-ingestible).

    The streaming, single-transcript twin of :func:`_write_transcript_models_gff3`
    (same per-model format), so output round-trips through
    :func:`_parse_transcript_models` / :func:`_iter_transcript_models`. Used by the
    max-intron split and the genomic SL cut, which stream one model at a time.
    """
    gid, gstart, gend = f"{tx.tid}.gene", tx.exons[0][0], tx.exons[-1][1]
    loc = f"{tx.chrom}\t{tx.source}"
    fh.write(f"{loc}\tgene\t{gstart}\t{gend}\t.\t{tx.strand}\t.\tID={gid}\n")
    fh.write(f"{loc}\tmRNA\t{gstart}\t{gend}\t.\t{tx.strand}\t.\tID={tx.tid};Parent={gid}\n")
    for k, (s, e) in enumerate(tx.exons, start=1):
        fh.write(
            f"{loc}\texon\t{s}\t{e}\t.\t{tx.strand}\t.\tID={tx.tid}.exon{k};Parent={tx.tid}\n"
        )


def _iter_transcript_models(gff: str | Path) -> Iterator[_Tx]:
    """Stream transcript models from a GFF3/GTF, one ``_Tx`` at a time.

    The streaming twin of :func:`_parse_transcript_models`. Groups consecutive
    ``exon`` rows by transcript id (``Parent=`` / ``transcript_id``), assuming each
    transcript's exon rows are contiguous — true for StringTie GTF and the
    BAM-derived GFF3. Streaming keeps memory bounded on a genome-wide model set.
    """
    cur: _Tx | None = None
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = _parse_attrs(cols[8])
            tid = attrs.get("Parent") or attrs.get("transcript_id")
            if not tid:
                continue
            if cur is None or cur.tid != tid:
                if cur is not None and cur.exons:
                    cur.exons.sort()
                    yield cur
                cur = _Tx(tid, cols[0], cols[6], cols[1])
            cur.exons.append((int(cols[3]), int(cols[4])))
    if cur is not None and cur.exons:
        cur.exons.sort()
        yield cur


def _split_models_at_clips(
    gff: str | Path, clip_map: dict[str, list[int]], min_segment: int
) -> tuple[list[_Tx], int]:
    """Parse transcript models from *gff* and split each at its clip positions.

    *clip_map* maps transcript id -> spliced-transcript clip positions (the same
    coordinate space as a read-to-spliced-FASTA mapping). Returns
    ``(out_models, n_split)``; transcripts absent from *clip_map* (or whose clips
    all fall in sub-*min_segment* pieces) pass through unchanged.
    """
    out_models: list[_Tx] = []
    n_split = 0
    for tx in _parse_transcript_models(gff):
        clips = clip_map.get(tx.tid)
        split = _split_transcript(tx, clips, min_segment) if clips else None
        if split:
            out_models.extend(split)
            n_split += 1
        else:
            out_models.append(tx)
    return out_models, n_split


def clip_gff3(
    gff: str | Path,
    genome: str | Path,
    clip_map: dict[str, list[int]],
    out_gff: str | Path,
    out_fasta: str | Path,
) -> None:
    """Split fused transcripts in a transcript>exon GFF3/GTF at *clip_map* positions.

    Writes the split models to *out_gff* (gene>mRNA>exon) and their re-spliced
    sequences to *out_fasta* (minus-strand reverse-complement handled by
    ``iter_assembled_sequences``). Transcripts with no clips pass through unchanged.
    """
    out_models, _ = _split_models_at_clips(gff, clip_map, _MIN_SEGMENT)
    _write_transcript_models_gff3(out_models, out_gff)
    db = create_gff_db(str(out_gff))
    with open(out_fasta, "w") as fh:
        for mrna, seq in iter_assembled_sequences(db, genome, child_featuretype="exon"):
            fh.write(f">{mrna.id}\n{seq}\n")


# ---------------------------------------------------------------------------
# Genome-guided StringTie GTF clip (read-pair troughs on spliced transcripts)
# ---------------------------------------------------------------------------


def resolve_stringtie_models(work_dir: Path) -> Path:
    """The StringTie model file downstream steps read.

    The jaccard-clipped ``stringtie.jaccard.gff3`` when the jaccard step produced it
    (so a fused StringTie locus is already split), else the raw ``stringtie.gtf``.
    Mirrors how ``map_transcripts`` prefers the de novo ``.jaccard.fasta`` sibling.
    """
    clipped = work_dir / STRINGTIE_JACCARD_GFF3
    if clipped.exists() and clipped.stat().st_size > 0:
        return clipped
    return work_dir / _STRINGTIE_GTF


def _write_spliced_fasta(models: list[_Tx], genome: str | Path, out_fasta: Path) -> None:
    """Write each model's spliced sequence (5'->3', RC on ``-``) keyed by transcript id.

    Same orientation convention as :func:`eukan.gff.io.iter_assembled_sequences`
    (ascending-genomic exon concat, reverse-complemented on the minus strand), so
    read-to-FASTA fragment spans land in the same 5'->3' spliced space the clip
    splitter (:func:`_split_transcript` / :func:`_partition_exons`) expects.
    """
    with ContigIndex(genome) as contigs, open(out_fasta, "w") as fh:
        for tx in models:
            seq = "".join(str(contigs[tx.chrom][s - 1 : e].seq) for s, e in tx.exons)
            if tx.strand == "-":
                seq = str(Seq(seq).reverse_complement())
            fh.write(f">{tx.tid}\n{seq}\n")


def _clip_stringtie_gtf(config: AssemblyConfig) -> tuple[int, int]:
    """Jaccard-clip the genome-guided StringTie GTF; return ``(n_models, n_split)``.

    Maps read pairs to StringTie's spliced transcripts, finds coverage troughs the
    same way as the de novo FASTA path (coverage-adaptive via ``jaccard_greediness``),
    and splits the genome-anchored models at those junctions into
    ``stringtie.jaccard.gff3``. A no-op (empty clip map) still rewrites the file so
    downstream consistently reads the clipped artifact.
    """
    wd = config.work_dir
    gtf = wd / _STRINGTIE_GTF
    models = _parse_transcript_models(gtf)
    out_gff = wd / STRINGTIE_JACCARD_GFF3
    if not models:
        _write_transcript_models_gff3([], out_gff)
        return 0, 0

    spliced = wd / "stringtie.spliced.fasta"
    _write_spliced_fasta(models, config.genome, spliced)
    bam = map_reads_to_transcripts(config, spliced, "stringtie")

    handle = pysam.AlignmentFile(str(bam), "rb")
    ref_len = dict(zip(handle.references, handle.lengths, strict=False))
    handle.close()

    clip_map: dict[str, list[int]] = {}
    for ref, spans in iter_fragment_spans(bam):
        if not spans:
            continue
        clips = find_clip_points(spans, ref_len[ref], **_clip_knobs(config))
        if clips:
            clip_map[ref] = clips

    out_models, n_split = _split_models_at_clips(gtf, clip_map, config.min_sl_fragment)
    _write_transcript_models_gff3(out_models, out_gff)
    spliced.unlink(missing_ok=True)
    return len(models), n_split
