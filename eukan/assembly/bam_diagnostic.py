"""Diagnostic walk over a BAM for soft-clips + intron dinucleotides.

Extracts and aggregates:
- Soft-clip ends (>= ``min_clip_len``) normalized to mRNA orientation
  (reverse-complement of clip and swapped side label for reverse-mapped
  reads), keyed by an *anchored* substring -- the K bp of the clip that
  touches the aligned read base. The anchored end is what's conserved
  in a splice-leader sequence and what's locus-specific in an
  intron-spillover clip, so this is the right end to cluster by.
- Distinct genomic loci ``(chrom, anchor_pos, side)`` per cluster, so
  callers can read off the "few clusters, many loci" trans-splicing
  signature vs. the "many clusters, one locus each" non-canonical
  splice signature.
- Intron splice-site dinucleotides for every N cigar op, looked up from
  the genome FASTA and rotated to mRNA orientation for reverse reads.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from eukan.assembly.polya import PolyAStats, tally_clip

if TYPE_CHECKING:
    import pysam


_CIGAR_M = 0
_CIGAR_D = 2
_CIGAR_N = 3
_CIGAR_SOFT_CLIP = 4
_CIGAR_EQ = 7
_CIGAR_X = 8
_CIGAR_REF_CONSUMING = frozenset([_CIGAR_M, _CIGAR_D, _CIGAR_N, _CIGAR_EQ, _CIGAR_X])

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Per-locus key cap. Statistically the dominant key arrives among the
# first few clips at any locus (≥94% prevalence at depth ≥ a few), so
# 5 distinct keys captures it plus a few error variants. Past that we
# stop tracking new distinct keys but keep incrementing the ones we have.
LOCUS_KEY_CAP = 5

_ALPHABET = "ACGT"

# compute_verdict heuristic thresholds (subject to revision as more datasets
# come in). Trans-splicing is called from the top non-trivial soft-clip
# cluster's locus + read counts; non-canonical splice from the canonical
# dinucleotide percentage.
_TS_STRONG_MIN_LOCI = 1000
_TS_STRONG_MIN_READS = 10_000
_TS_MODERATE_MIN_LOCI = 100
_TS_MODERATE_MIN_READS = 1000
_NC_CANONICAL_PCT_EXTENSIVE = 80.0  # canonical% below this -> EXTENSIVE
_NC_CANONICAL_PCT_MODERATE = 95.0   # canonical% below this -> MODERATE, else ABSENT
_NC_DINUC_REPORT_MIN_PCT = 1.0      # min % for a non-canonical dinuc to be named


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _hamming_distance(a: str, b: str) -> int:
    """Hamming distance between equal-length strings.

    Returns ``max(len(a), len(b))`` (i.e. "definitely far") if the
    strings differ in length — keeps callers from special-casing the
    sub-K short keys.
    """
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b, strict=True) if x != y)


def _hamming1_neighbors(key: str) -> Iterable[str]:
    """All strings at exactly Hamming distance 1 from ``key`` over ACGT."""
    for i, ch in enumerate(key):
        prefix = key[:i]
        suffix = key[i + 1:]
        for b in _ALPHABET:
            if b != ch:
                yield prefix + b + suffix


def _hamming2_neighbors(key: str) -> Iterable[str]:
    """All strings at exactly Hamming distance 2 from ``key`` over ACGT."""
    n = len(key)
    for i in range(n):
        ci = key[i]
        for j in range(i + 1, n):
            cj = key[j]
            for bi in _ALPHABET:
                if bi == ci:
                    continue
                for bj in _ALPHABET:
                    if bj == cj:
                        continue
                    yield key[:i] + bi + key[i + 1:j] + bj + key[j + 1:]


def _canonicalize_keys(
    keys_with_counts: list[tuple[str, int]],
    hamming_tolerance: int,
) -> dict[str, str]:
    """Online clustering: each raw key → canonical seed within Hamming H.

    Process keys in descending order of read support so that high-coverage
    motifs become seeds and their low-coverage near-neighbors (mostly
    sequencer-error variants) get absorbed. Greedy + order-dependent but
    deterministic given the input sort.

    With ``hamming_tolerance == 0`` this is the identity mapping.
    """
    canonical: dict[str, str] = {}
    seeds: set[str] = set()

    if hamming_tolerance <= 0:
        for key, _ in keys_with_counts:
            canonical[key] = key
        return canonical

    for key, _ in keys_with_counts:
        if key in seeds:
            canonical[key] = key
            continue
        chosen: str | None = None
        for nb in _hamming1_neighbors(key):
            if nb in seeds:
                chosen = nb
                break
        if chosen is None and hamming_tolerance >= 2:
            for nb in _hamming2_neighbors(key):
                if nb in seeds:
                    chosen = nb
                    break
        if chosen is not None:
            canonical[key] = chosen
        else:
            seeds.add(key)
            canonical[key] = key
    return canonical


def _cluster_consensus(
    clip_seqs: list[str],
    side: str,
    *,
    min_coverage: int = 5,
    min_majority_fraction: float = 0.6,
    max_extension: int = 60,
) -> str:
    """Majority-base consensus across clips, anchored at the alignment edge.

    For 5p clips the anchor is the rightmost base of the clip (the bp
    touching the alignment); the consensus extends leftward (upstream).
    For 3p clips the anchor is the leftmost base; extension is rightward
    (downstream).

    Walks columns outward from the anchor. At each column, takes the
    majority base if coverage >= ``min_coverage`` AND
    majority/coverage >= ``min_majority_fraction``. Stops at the first
    failing column or at ``max_extension``. Non-ACGT bases (``N`` etc.)
    don't contribute to the count for that column. Returns the consensus
    in 5'->3' order (the anchor base is at the right end for 5p, left
    end for 3p). Empty string if column 0 fails.
    """
    if not clip_seqs:
        return ""
    upper = [s.upper() for s in clip_seqs]

    bases: list[str] = []
    for col in range(max_extension):
        counts: Counter[str] = Counter()
        for s in upper:
            if side == "5p":
                idx = len(s) - 1 - col
                if idx < 0:
                    continue
            else:
                idx = col
                if idx >= len(s):
                    continue
            ch = s[idx]
            if ch in _ALPHABET:
                counts[ch] += 1
        coverage = sum(counts.values())
        if coverage < min_coverage:
            break
        best, n_best = counts.most_common(1)[0]
        if n_best / coverage < min_majority_fraction:
            break
        bases.append(best)

    if not bases:
        return ""
    # For 5p, columns were collected right-to-left (anchor outward).
    # Reverse to get 5'->3' order. For 3p, the order is already 5'->3'.
    if side == "5p":
        return "".join(reversed(bases))
    return "".join(bases)


@dataclass
class SoftClipStats:
    """Aggregated soft-clip stats from one BAM walk."""

    n_reads_scanned: int
    n_reads_with_clip: int
    n_clips_total: int
    n_clips_by_side: dict[str, int]
    n_loci: int
    n_clusters: int
    cluster_to_loci: dict[str, int] = field(default_factory=dict)
    cluster_to_reads: dict[str, int] = field(default_factory=dict)
    cluster_examples: dict[str, str] = field(default_factory=dict)
    cluster_consensus: dict[str, str] = field(default_factory=dict)
    top_clusters: list[tuple[str, int, int]] = field(default_factory=list)


@dataclass
class IntronStats:
    """Aggregated intron dinucleotide stats from one BAM walk."""

    n_introns_total: int
    by_dinucleotide: dict[str, int]
    canonical_pct: float


@dataclass
class LocusData:
    """Per-locus state accumulated during the BAM walk.

    ``long_key_counts`` is capped at ``LOCUS_KEY_CAP`` distinct keys, but
    once a key is tracked its count keeps incrementing — so the dominant
    key (which arrives early and often) remains accurate.
    """

    n_clips: int = 0
    longest_clip: str = ""
    long_key_counts: Counter[str] = field(default_factory=Counter)
    short_keys: set[str] = field(default_factory=set)


@dataclass
class LocusRow:
    """One emitted row of per-locus data (non-singleton loci only)."""

    chrom: str
    pos: int
    side: str
    n_clips: int
    status: str          # "consistent" | "inconsistent" | "short_only"
    motif_key: str       # "" for inconsistent / short_only
    motif_share: int     # 0 for inconsistent / short_only
    longest_clip: str


@dataclass
class LocusConsistencyStats:
    """Whether clips at each locus consolidate to a single motif.

    The user's biological hypothesis: at a non-canonical splice site,
    reads of varying clip lengths at the same locus all consolidate to
    one (locus-specific) intron-sequence motif. Trans-splicing looks
    the same locally (one motif per locus) but the motif is shared
    across many loci. So we classify each locus by consolidation
    status, and for consistent loci bucket by how widely the motif is
    shared across the genome.
    """

    n_loci_total: int
    n_loci_singleton: int
    n_loci_consistent: int
    n_loci_inconsistent: int
    n_loci_short_only: int
    motif_share_histogram: dict[str, int] = field(default_factory=dict)
    deepest_loci: list[LocusRow] = field(default_factory=list)
    all_rows: list[LocusRow] = field(default_factory=list)


@dataclass
class DiagnosticReport:
    """All views of one BAM walk."""

    softclip: SoftClipStats
    intron: IntronStats
    locus_consistency: LocusConsistencyStats
    # Poly-A / poly-T soft-clip tallies, accumulated in the same pass but kept
    # separate from the SL/trans-splice verdict (see eukan.assembly.polya).
    polya: PolyAStats = field(default_factory=PolyAStats)


@dataclass
class TransSplicingCall:
    """Categorical verdict on whether trans-splicing is present."""

    call: str  # "STRONG" | "MODERATE" | "ABSENT"
    top_non_trivial_cluster_key: str
    top_non_trivial_cluster_consensus: str
    top_non_trivial_cluster_n_loci: int
    top_non_trivial_cluster_n_reads: int
    sl_bucket_pct_of_consistent: float


@dataclass
class NonCanonicalSpliceCall:
    """Categorical verdict on non-canonical splice prevalence."""

    call: str  # "EXTENSIVE" | "MODERATE" | "ABSENT"
    canonical_pct: float
    top_non_canonical_dinuc: str


@dataclass
class Verdict:
    """Empirical-verdict labels + supporting numbers from a BAM walk."""

    trans_splicing: TransSplicingCall
    non_canonical_splice: NonCanonicalSpliceCall


def _iter_primary_alignments(
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


def _extract_clips(
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
    if op0 == _CIGAR_SOFT_CLIP and len0 >= min_clip_len:
        clip_seq = seq[:len0]
        anchor_pos = read.reference_start
        if read.is_reverse:
            yield "3p", _reverse_complement(clip_seq), anchor_pos
        else:
            yield "5p", clip_seq, anchor_pos

    op_last, len_last = cigar[-1]
    if op_last == _CIGAR_SOFT_CLIP and len_last >= min_clip_len:
        clip_seq = seq[-len_last:]
        anchor_pos = ref_end - 1
        if read.is_reverse:
            yield "5p", _reverse_complement(clip_seq), anchor_pos
        else:
            yield "3p", clip_seq, anchor_pos


def _cluster_key(side: str, seq: str, k: int) -> str:
    """Return the K-bp substring anchored at the alignment boundary."""
    seq = seq.upper()
    if len(seq) <= k:
        return seq
    if side == "5p":
        return seq[-k:]
    return seq[:k]


def _walk_introns(
    read: pysam.AlignedSegment,
) -> Iterable[tuple[str, int, int]]:
    """Yield ``(chrom, intron_start, intron_end_exclusive)`` for each N op."""
    chrom = read.reference_name
    cigar = read.cigartuples
    if chrom is None or cigar is None:
        return
    pos = read.reference_start
    for op, length in cigar:
        if op == _CIGAR_N:
            yield chrom, pos, pos + length
        if op in _CIGAR_REF_CONSUMING:
            pos += length


def _short_fits_long(short: str, long: str, side: str) -> bool:
    """Is ``short`` consistent with ``long`` for the given clip side?

    A length-<K clip at a locus consolidates with the longer K-bp
    anchored key if its full sequence is a suffix (5p) / prefix (3p) of
    the K-bp key — same alignment boundary, just less intron read into.
    """
    if side == "5p":
        return long.endswith(short)
    return long.startswith(short)


def _classify_locus(
    ld: LocusData, side: str, hamming_tolerance: int,
    min_consistency_fraction: float,
) -> tuple[str, str]:
    """Classify a locus's consolidation status.

    Returns ``(status, motif_key)`` where status is one of
    ``"singleton"``, ``"consistent"``, ``"inconsistent"``,
    ``"short_only"``. ``motif_key`` is the locus's K-bp anchored motif
    (the dominant key by read count) when status == "consistent", else
    an empty string.

    A locus is consistent if (a) every short clip fits the dominant as
    suffix/prefix, and (b) the fraction of tracked long-clip reads
    whose K-bp key is within ``hamming_tolerance`` of the dominant key
    meets ``min_consistency_fraction``. Fraction-based is more robust
    than "all-must-match" at high depth, where a small minority of
    multi-error reads can otherwise flip a locus to inconsistent.
    """
    if ld.n_clips == 1:
        return "singleton", ""

    if ld.long_key_counts:
        # Dominant = most-common tracked long key. Ties broken
        # deterministically by Counter.most_common (insertion order).
        dominant_key, _ = ld.long_key_counts.most_common(1)[0]
        for sk in ld.short_keys:
            if not _short_fits_long(sk, dominant_key, side):
                return "inconsistent", ""

        matched = 0
        total = 0
        for k, c in ld.long_key_counts.items():
            total += c
            if _hamming_distance(k, dominant_key) <= hamming_tolerance:
                matched += c

        if total == 0 or matched / total < min_consistency_fraction:
            return "inconsistent", ""
        return "consistent", dominant_key

    # No long keys at all — only sub-K clips at this locus.
    if len(ld.short_keys) <= 1:
        return "short_only", ""

    longest_short = max(ld.short_keys, key=len)
    for sk in ld.short_keys:
        if sk == longest_short:
            continue
        if not _short_fits_long(sk, longest_short, side):
            return "inconsistent", ""
    return "short_only", ""


def _share_bucket(share: int) -> str:
    """Map a motif's locus-share count into a coarse histogram bucket."""
    if share <= 1:
        return "1"
    if share <= 10:
        return "2-10"
    if share <= 100:
        return "11-100"
    if share <= 1000:
        return "101-1000"
    return ">1000"


def _dinucleotide(
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


def diagnose_bam(
    bam_path: Path,
    genome_path: Path,
    *,
    min_clip_len: int = 8,
    cluster_key_len: int = 12,
    min_mapq: int = 20,
    top_n: int = 15,
    hamming_tolerance: int = 2,
    cluster_hamming_tolerance: int = 1,
    min_consistency_fraction: float = 0.95,
    consensus_min_majority_fraction: float = 0.6,
) -> DiagnosticReport:
    """Stream a BAM and return aggregated soft-clip + intron + locus stats.

    One sequential pass over the BAM; no index required. Three knobs
    tune motif consolidation over the same K-bp anchored window:

    - ``hamming_tolerance`` — per-locus consolidation. Tracked clip
      keys at a locus are considered to match the dominant key if
      within this distance.
    - ``min_consistency_fraction`` — per-locus fraction-based check.
      A locus is consistent if at least this fraction of tracked
      long-clip reads have a key within ``hamming_tolerance`` of the
      dominant. Tolerates a minority of multi-error reads at deep
      loci that would otherwise flip the locus to "inconsistent".
    - ``cluster_hamming_tolerance`` — cross-locus canonicalization.
      Raw cluster keys within this distance of an earlier-seen-by-count
      seed get merged. Larger values absorb sequencer-error
      neighbors and biological paralogs but raise random-collision
      risk in cross-locus share counts.

    Set the Hamming knobs to 0 to disable; set the fraction to 1.0 to
    require all-tracked-keys-must-match.

    ``consensus_min_majority_fraction`` controls how strict the
    per-column majority vote is when building a cluster's consensus
    sequence. Lower values let the consensus extend further into noisier
    columns; higher values terminate sooner. Default 0.6.
    """
    import pysam

    from eukan.infra.genome import ContigIndex

    n_reads_scanned = 0
    n_reads_with_clip = 0
    clips_by_side: Counter[str] = Counter()
    cluster_loci: dict[str, set[tuple[str, int, str]]] = {}
    cluster_reads: Counter[str] = Counter()
    cluster_examples: dict[str, str] = {}
    locus_data: dict[tuple[str, int, str], LocusData] = {}

    dinuc: Counter[str] = Counter()
    n_introns = 0

    polya = PolyAStats()

    with ContigIndex(genome_path) as contigs:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
        try:
            for read in _iter_primary_alignments(bam, min_mapq=min_mapq):
                n_reads_scanned += 1
                had_clip = False
                chrom_name = read.reference_name
                if chrom_name is None:
                    continue
                for side, seq, anchor_pos in _extract_clips(read, min_clip_len):
                    had_clip = True
                    clips_by_side[side] += 1
                    tally_clip(polya, side, seq, chrom_name)
                    locus = (chrom_name, anchor_pos, side)
                    key = _cluster_key(side, seq, cluster_key_len)
                    cluster_loci.setdefault(key, set()).add(locus)
                    cluster_reads[key] += 1
                    cluster_examples.setdefault(key, seq)

                    ld = locus_data.setdefault(locus, LocusData())
                    ld.n_clips += 1
                    if len(seq) > len(ld.longest_clip):
                        ld.longest_clip = seq
                    if len(seq) >= cluster_key_len:
                        if key in ld.long_key_counts or len(ld.long_key_counts) < LOCUS_KEY_CAP:
                            ld.long_key_counts[key] += 1
                    else:
                        if len(ld.short_keys) < LOCUS_KEY_CAP:
                            ld.short_keys.add(seq.upper())
                if had_clip:
                    n_reads_with_clip += 1

                for chrom, intron_start, intron_end in _walk_introns(read):
                    n_introns += 1
                    pair = _dinucleotide(contigs, chrom, intron_start, intron_end)
                    if pair is not None:
                        dinuc[pair] += 1
        finally:
            bam.close()

    # Canonicalize cluster keys: high-coverage motifs become seeds first
    # so their sequencer-error neighbors get absorbed into them.
    keys_sorted = sorted(
        cluster_reads.items(), key=lambda kv: (-kv[1], kv[0]),
    )
    canonical_map = _canonicalize_keys(keys_sorted, cluster_hamming_tolerance)

    canon_cluster_loci: dict[str, set[tuple[str, int, str]]] = {}
    canon_cluster_reads: Counter[str] = Counter()
    canon_cluster_examples: dict[str, str] = {}
    # Track which raw key contributed each example, so a later high-coverage
    # raw key can replace a lower-coverage one's example.
    canon_example_src_reads: dict[str, int] = {}
    for raw_key, loci in cluster_loci.items():
        seed = canonical_map[raw_key]
        canon_cluster_loci.setdefault(seed, set()).update(loci)
        n_reads_here = cluster_reads[raw_key]
        canon_cluster_reads[seed] += n_reads_here
        if n_reads_here > canon_example_src_reads.get(seed, -1):
            canon_cluster_examples[seed] = cluster_examples[raw_key]
            canon_example_src_reads[seed] = n_reads_here

    canon_cluster_to_loci_count = {k: len(v) for k, v in canon_cluster_loci.items()}
    top_keys = sorted(
        canon_cluster_to_loci_count.items(),
        key=lambda kv: (-kv[1], -canon_cluster_reads[kv[0]], kv[0]),
    )[:top_n]
    top_clusters = [(k, n_loci, canon_cluster_reads[k]) for k, n_loci in top_keys]

    # Build a consensus sequence for each top-N cluster from its member
    # longest-clip sequences on the dominant side. Only top-N to bound
    # work on long-tail clusters whose consensus is never reported.
    # ``LocusData.longest_clip`` is the longest soft-clip of any read at
    # the locus, regardless of its K-bp anchor. At a deep SL locus, an
    # unrelated long contamination read can win that slot and pollute the
    # consensus. Filter to clips whose own K-bp anchor is one of the
    # cluster's raw keys — keeps only reads that ARE the SL/motif, not
    # co-located off-anchor contamination.
    canon_to_raw: dict[str, set[str]] = {}
    for raw_key, seed in canonical_map.items():
        canon_to_raw.setdefault(seed, set()).add(raw_key)

    canon_cluster_consensus: dict[str, str] = {}
    for seed, _n_loci, _n_reads in top_clusters:
        loci = canon_cluster_loci.get(seed, set())
        raw_keys = canon_to_raw.get(seed, {seed})
        loci_by_side: dict[str, list[tuple[str, int, str]]] = {"5p": [], "3p": []}
        for loc in loci:
            loci_by_side[loc[2]].append(loc)
        side_counts = Counter({s: len(v) for s, v in loci_by_side.items()})
        if not any(side_counts.values()):
            continue
        dominant_side = side_counts.most_common(1)[0][0]
        clips: list[str] = []
        for loc in loci_by_side[dominant_side]:
            lc = locus_data[loc].longest_clip
            if len(lc) < cluster_key_len:
                continue
            if _cluster_key(dominant_side, lc, cluster_key_len) not in raw_keys:
                continue
            clips.append(lc)
        if not clips:
            continue
        consensus = _cluster_consensus(
            clips, dominant_side,
            min_majority_fraction=consensus_min_majority_fraction,
        )
        if consensus:
            canon_cluster_consensus[seed] = consensus

    soft = SoftClipStats(
        n_reads_scanned=n_reads_scanned,
        n_reads_with_clip=n_reads_with_clip,
        n_clips_total=sum(clips_by_side.values()),
        n_clips_by_side=dict(clips_by_side),
        n_loci=len(locus_data),
        n_clusters=len(canon_cluster_to_loci_count),
        cluster_to_loci=canon_cluster_to_loci_count,
        cluster_to_reads=dict(canon_cluster_reads),
        cluster_examples=canon_cluster_examples,
        cluster_consensus=canon_cluster_consensus,
        top_clusters=top_clusters,
    )

    canonical = dinuc.get("GT-AG", 0) + dinuc.get("CT-AC", 0)
    canonical_pct = 100.0 * canonical / n_introns if n_introns else 0.0
    intron = IntronStats(
        n_introns_total=n_introns,
        by_dinucleotide=dict(dinuc.most_common()),
        canonical_pct=canonical_pct,
    )

    locus_consistency = _build_locus_consistency_stats(
        locus_data, canon_cluster_to_loci_count, canonical_map,
        hamming_tolerance=hamming_tolerance,
        min_consistency_fraction=min_consistency_fraction,
        top_n=top_n,
    )

    return DiagnosticReport(
        softclip=soft, intron=intron, locus_consistency=locus_consistency, polya=polya,
    )


def _build_locus_consistency_stats(
    locus_data: dict[tuple[str, int, str], LocusData],
    canon_cluster_to_loci_count: dict[str, int],
    canonical_map: dict[str, str],
    *,
    hamming_tolerance: int,
    min_consistency_fraction: float,
    top_n: int,
) -> LocusConsistencyStats:
    """Classify every locus and build histograms + per-row TSV input.

    For each consistent locus, the emitted motif_key is the *canonical
    seed* of the dominant raw key, so its motif_share lookup uses the
    merged cross-locus view.
    """
    n_singleton = n_consistent = n_inconsistent = n_short_only = 0
    histogram: Counter[str] = Counter()
    all_rows: list[LocusRow] = []

    for locus, ld in locus_data.items():
        chrom, pos, side = locus
        status, dominant_key = _classify_locus(
            ld, side, hamming_tolerance, min_consistency_fraction,
        )

        if status == "singleton":
            n_singleton += 1
            continue

        canonical_motif = canonical_map.get(dominant_key, dominant_key) if dominant_key else ""
        motif_share = (
            canon_cluster_to_loci_count.get(canonical_motif, 0) if canonical_motif else 0
        )
        all_rows.append(
            LocusRow(
                chrom=chrom,
                pos=pos,
                side=side,
                n_clips=ld.n_clips,
                status=status,
                motif_key=canonical_motif,
                motif_share=motif_share,
                longest_clip=ld.longest_clip,
            )
        )

        if status == "consistent":
            n_consistent += 1
            histogram[_share_bucket(motif_share)] += 1
        elif status == "inconsistent":
            n_inconsistent += 1
        else:
            n_short_only += 1

    deepest = sorted(all_rows, key=lambda r: -r.n_clips)[:top_n]

    return LocusConsistencyStats(
        n_loci_total=len(locus_data),
        n_loci_singleton=n_singleton,
        n_loci_consistent=n_consistent,
        n_loci_inconsistent=n_inconsistent,
        n_loci_short_only=n_short_only,
        motif_share_histogram=dict(histogram),
        deepest_loci=deepest,
        all_rows=all_rows,
    )


# ---------------------------------------------------------------------------
# Empirical verdict — categorical calls on top of the raw stats
# ---------------------------------------------------------------------------


def _is_low_complexity(key: str, threshold: float = 0.7) -> bool:
    """True if any single base accounts for more than ``threshold`` of ``key``.

    Used to skip poly-A / poly-T / similar near-mononucleotide cluster keys
    when picking the "top non-trivial cluster" for the verdict. The 0.7
    threshold passes SL motifs like ``CTGTACTTTATT`` (T=5/12=0.42) and
    real intron sequences but rejects pure runs.
    """
    if not key:
        return False
    counts = Counter(key)
    return max(counts.values()) / len(key) > threshold


def _first_non_trivial_top_cluster(
    soft: SoftClipStats,
) -> tuple[str, int, int] | None:
    """Walk ``soft.top_clusters`` and return the first non-low-complexity row."""
    for key, n_loci, n_reads in soft.top_clusters:
        if not _is_low_complexity(key):
            return key, n_loci, n_reads
    return None


def compute_verdict(report: DiagnosticReport) -> Verdict:
    """Derive categorical empirical-verdict calls + supporting numbers.

    Heuristic thresholds (subject to revision as more datasets come in):

    Trans-splicing:
      - STRONG    if top non-trivial cluster ≥ 1000 loci AND ≥ 10,000 reads
      - MODERATE  if top non-trivial cluster ≥  100 loci AND ≥  1,000 reads
      - ABSENT    otherwise

    Non-canonical splice:
      - EXTENSIVE if canonical_pct < 80%
      - MODERATE  if 80% ≤ canonical_pct < 95%
      - ABSENT    if canonical_pct ≥ 95%

    The supporting numbers (top-cluster key/loci/reads, sl-bucket %,
    canonical_pct, top non-canonical dinuc) are kept alongside the
    label so callers can override the judgement without re-reading
    the spec.
    """
    soft = report.softclip
    intron = report.intron
    lc = report.locus_consistency

    top = _first_non_trivial_top_cluster(soft)
    if top is not None:
        key, n_loci, n_reads = top
        if n_loci >= _TS_STRONG_MIN_LOCI and n_reads >= _TS_STRONG_MIN_READS:
            ts_call = "STRONG"
        elif n_loci >= _TS_MODERATE_MIN_LOCI and n_reads >= _TS_MODERATE_MIN_READS:
            ts_call = "MODERATE"
        else:
            ts_call = "ABSENT"
    else:
        key, n_loci, n_reads = "", 0, 0
        ts_call = "ABSENT"

    sl_bucket_pct = 0.0
    if lc.n_loci_consistent:
        sl_bucket_pct = 100.0 * lc.motif_share_histogram.get(">1000", 0) / lc.n_loci_consistent

    if intron.canonical_pct < _NC_CANONICAL_PCT_EXTENSIVE:
        nc_call = "EXTENSIVE"
    elif intron.canonical_pct < _NC_CANONICAL_PCT_MODERATE:
        nc_call = "MODERATE"
    else:
        nc_call = "ABSENT"

    nc_top_label = "none above 1%"
    for dinuc, n in intron.by_dinucleotide.items():
        if dinuc in ("GT-AG", "CT-AC"):
            continue
        pct = 100.0 * n / intron.n_introns_total if intron.n_introns_total else 0.0
        if pct >= _NC_DINUC_REPORT_MIN_PCT:
            nc_top_label = f"{dinuc} {pct:.2f}%"
        break

    return Verdict(
        trans_splicing=TransSplicingCall(
            call=ts_call,
            top_non_trivial_cluster_key=key,
            top_non_trivial_cluster_consensus=soft.cluster_consensus.get(key, ""),
            top_non_trivial_cluster_n_loci=n_loci,
            top_non_trivial_cluster_n_reads=n_reads,
            sl_bucket_pct_of_consistent=sl_bucket_pct,
        ),
        non_canonical_splice=NonCanonicalSpliceCall(
            call=nc_call,
            canonical_pct=intron.canonical_pct,
            top_non_canonical_dinuc=nc_top_label,
        ),
    )


def to_summary_dict(report: DiagnosticReport, verdict: Verdict) -> dict:
    """Pack a slim summary of the diagnostic + verdict for JSON output.

    Includes scalars + histograms + top-15 clusters + top-15 deepest
    loci + the full verdict — but NOT the per-cluster ``cluster_to_*``
    maps or the per-locus ``all_rows`` (those are too heavy for a
    routine pipeline summary; the experiment driver still dumps them
    as TSV side files).
    """
    soft = report.softclip
    intron = report.intron
    lc = report.locus_consistency
    return {
        "softclip": {
            "n_reads_scanned": soft.n_reads_scanned,
            "n_reads_with_clip": soft.n_reads_with_clip,
            "n_clips_total": soft.n_clips_total,
            "n_clips_by_side": soft.n_clips_by_side,
            "n_loci": soft.n_loci,
            "n_clusters": soft.n_clusters,
            "top_clusters": [
                {
                    "key": k, "n_loci": n_loci, "n_reads": n_reads,
                    "example": soft.cluster_examples.get(k, ""),
                    "consensus": soft.cluster_consensus.get(k, ""),
                }
                for k, n_loci, n_reads in soft.top_clusters
            ],
        },
        "intron": {
            "n_introns_total": intron.n_introns_total,
            "canonical_pct": intron.canonical_pct,
            "by_dinucleotide": intron.by_dinucleotide,
        },
        "locus_consistency": {
            "n_loci_total": lc.n_loci_total,
            "n_loci_singleton": lc.n_loci_singleton,
            "n_loci_consistent": lc.n_loci_consistent,
            "n_loci_inconsistent": lc.n_loci_inconsistent,
            "n_loci_short_only": lc.n_loci_short_only,
            "motif_share_histogram": lc.motif_share_histogram,
            "deepest_loci": [
                {
                    "chrom": r.chrom, "pos": r.pos, "side": r.side,
                    "n_clips": r.n_clips, "status": r.status,
                    "motif_key": r.motif_key, "motif_share": r.motif_share,
                    "longest_clip": r.longest_clip,
                }
                for r in lc.deepest_loci
            ],
        },
        "verdict": {
            "trans_splicing": {
                "call": verdict.trans_splicing.call,
                "top_non_trivial_cluster_key":
                    verdict.trans_splicing.top_non_trivial_cluster_key,
                "top_non_trivial_cluster_consensus":
                    verdict.trans_splicing.top_non_trivial_cluster_consensus,
                "top_non_trivial_cluster_n_loci":
                    verdict.trans_splicing.top_non_trivial_cluster_n_loci,
                "top_non_trivial_cluster_n_reads":
                    verdict.trans_splicing.top_non_trivial_cluster_n_reads,
                "sl_bucket_pct_of_consistent":
                    verdict.trans_splicing.sl_bucket_pct_of_consistent,
            },
            "non_canonical_splice": {
                "call": verdict.non_canonical_splice.call,
                "canonical_pct": verdict.non_canonical_splice.canonical_pct,
                "top_non_canonical_dinuc":
                    verdict.non_canonical_splice.top_non_canonical_dinuc,
            },
        },
    }
