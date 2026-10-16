"""Jaccard clipping of fused transcripts — a standalone reimplementation of
Trinity's ``--jaccard_clip`` post-processing.

In gene-dense genomes with overlapping UTRs an assembler can string two
independently-transcribed genes into a single contig. Trinity's
``--jaccard_clip`` detects these by mapping read pairs back to the assembled
contigs and cutting where no read-pair fragment *bridges* a position. This
module applies the same logic, as a standalone step, to both Trinity transcript
sets (de novo and genome-guided; Trinity's own ``--jaccard_clip`` is not passed),
so the consolidation step (combinr) never sees an un-clipped fused contig.

Algorithm (faithful to Trinity, see ``trinityrnaseq/util/support_scripts/``):

1. Map paired reads to the transcript FASTA, ungapped (minimap2 ``-x sr``, the
   non-spliced short-read preset). One transcript == one BAM reference.
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

DORMANT with the unwired stringtie module (Trinity replaced StringTie): the
genome-anchored StringTie GTF is clipped the same way (:func:`_clip_stringtie_gtf`):
reads are mapped to StringTie's spliced transcripts, the same coverage troughs are
found, and the transcript>exon models are split at them (:func:`_split_models_at_clips`)
into ``stringtie.jaccard.gff3`` — so a StringTie transcript fusing two adjacent loci is
broken just as the de novo path keeps them separate. :func:`resolve_stringtie_models`
prefers that clipped GFF3 over the raw GTF (no live step calls it).

Jaccard clipping needs read PAIRS; with single-end reads the step is a no-op.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

import pysam
from Bio import SeqIO

from eukan.assembly.minimap2 import map_reads_to_transcripts
from eukan.assembly.tracks import TRANSCRIPT_SETS, jaccard_output_name
from eukan.gff.transcripts import Transcript, read_transcripts, split_at, stitch, write_transcripts
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
_REPOSITION_HALF_WIN = _TROUGH_WIN // 2  # coverage-reposition search radius (100)
# Ceiling for the coverage-adapted trough gate: below the depth where a clean
# junction would need a jaccard above this, the bridging signal is too weak to
# tell a fusion from noise, so the gate stops relaxing (never over-clips a real
# single transcript that simply has a thin patch).
_MAX_ADAPTIVE_TROUGH = 0.30

# DORMANT (StringTie path): the genome-guided StringTie GTF and the clipped sibling
# the jaccard step writes for it. Not produced by the active Trinity pipeline, so
# the run_jaccard StringTie branch and resolve_stringtie_models below are no-ops
# unless the kept-but-unwired stringtie module is re-activated.
_STRINGTIE_GTF = "stringtie.gtf"
STRINGTIE_JACCARD_GFF3 = "stringtie.jaccard.gff3"


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
# Read mapping (minimap2 -x sr, ungapped) + fragment extraction
# ---------------------------------------------------------------------------


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
    """Jaccard-clip both Trinity transcript FASTAs (and a dormant StringTie GTF).

    No-op on single-end reads (read-pair bridging is the clip signal).
    """
    wd = config.work_dir
    if not (config.left_reads and config.right_reads):
        log.warning(
            "Jaccard clipping needs paired reads (read-pair bridging is the "
            "signal); single-end input — skipping."
        )
        return

    # Each track's Trinity FASTA is rewritten into a ``.jaccard.fasta`` sibling,
    # which map_transcripts prefers when present (tracks.resolve_query).
    for name, _ in TRANSCRIPT_SETS:
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


def _split_models_at_clips(
    gff: str | Path, clip_map: dict[str, list[int]], min_segment: int
) -> tuple[list[Transcript], int]:
    """Parse transcript models from *gff* and split each at its clip positions.

    *clip_map* maps transcript id -> spliced-transcript clip positions (the same
    coordinate space as a read-to-spliced-FASTA mapping). Returns
    ``(out_models, n_split)``; transcripts absent from *clip_map* (or whose clips
    all fall in sub-*min_segment* pieces) pass through unchanged.
    """
    out_models: list[Transcript] = []
    n_split = 0
    for tx in read_transcripts(gff):
        clips = clip_map.get(tx.tid)
        split = split_at(tx, clips, min_segment) if clips else None
        if split:
            out_models.extend(split)
            n_split += 1
        else:
            out_models.append(tx)
    return out_models, n_split


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


def _write_spliced_fasta(
    models: list[Transcript], genome: str | Path, out_fasta: Path
) -> None:
    """Write each model's spliced sequence (5'->3', RC on ``-``) keyed by transcript id.

    Same orientation convention as :func:`eukan.gff.io.iter_assembled_sequences`
    (ascending-genomic exon concat, reverse-complemented on the minus strand), so
    read-to-FASTA fragment spans land in the same 5'->3' spliced space the clip
    splitter (:func:`eukan.gff.transcripts.split_at`) expects.
    """
    with ContigIndex(genome) as contigs, open(out_fasta, "w") as fh:
        for tx in models:
            fh.write(f">{tx.tid}\n{stitch(tx, contigs)}\n")


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
    models = read_transcripts(gtf)
    out_gff = wd / STRINGTIE_JACCARD_GFF3
    if not models:
        write_transcripts([], out_gff)
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
    write_transcripts(out_models, out_gff)
    spliced.unlink(missing_ok=True)
    return len(models), n_split
