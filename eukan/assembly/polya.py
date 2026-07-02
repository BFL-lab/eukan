"""Poly-A / poly-T characterization of soft-clips and assembled transcripts.

Standalone, SL-independent statistics on where poly-A tails surface in the
assembly pipeline. Tools that align transcripts **pairwise** (gmap/blat) recommend
trimming poly-A tails first, where an untrimmed tail forces terminal mismatches/insertions
and degrades the alignment. eukan instead maps with STAR/segemehl in ``Local``
(soft-clip) mode, so a poly-A tail simply *soft-clips* rather than degrading the
body alignment. This module quantifies that, so the choice can be revisited from
data rather than assumed.

A poly-A tail sits at the mRNA 3' end, so it is an **A-rich ``3p`` soft-clip**; a
5' poly-T (antisense poly-A) is a **T-rich ``5p`` clip**. The clip sequences fed in
here are expected in mRNA 5'->3' orientation, exactly as
:func:`eukan.assembly.bam_diagnostic._extract_clips` yields them.

The output is a per-section ``polyA_diagnostic.json`` (read-BAM "reads", de novo
transcript-BAM "transcripts", and "unmapped_transcripts" tails) plus an INFO log —
deliberately decoupled from the SL trans-splice verdict.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from eukan.infra.logging import get_logger

log = get_logger(__name__)

# Well-tuned defaults (module constants, not config knobs — mirrors jaccard).
POLYA_MIN_LEN = 8     # shortest soft-clip considered (matches diagnose_bam min_clip_len)
POLYA_MIN_FRAC = 0.8  # fraction of the clip that must be the homopolymer base

# Diagnostic JSON filename (also registered as Artifact.POLYA_DIAGNOSTIC).
POLYA_DIAGNOSTIC = "polyA_diagnostic.json"


def classify_clip(
    side: str, seq: str, *, min_len: int = POLYA_MIN_LEN, min_frac: float = POLYA_MIN_FRAC
) -> str | None:
    """``"polyA"`` for an A-rich 3' clip, ``"polyT"`` for a T-rich 5' clip, else ``None``.

    *seq* must be in mRNA 5'->3' orientation (as ``_extract_clips`` yields): a poly-A
    tail is at the mRNA 3' end (``side == "3p"``) and reads as A's; a 5' poly-T
    (antisense poly-A) is at ``side == "5p"`` and reads as T's. Shorter than
    *min_len* or below *min_frac* homopolymer content → not a poly-tail.
    """
    if len(seq) < min_len:
        return None
    s = seq.upper()
    if side == "3p" and s.count("A") / len(s) >= min_frac:
        return "polyA"
    if side == "5p" and s.count("T") / len(s) >= min_frac:
        return "polyT"
    return None


@dataclass
class PolyAStats:
    """Poly-A / poly-T soft-clip tallies accumulated over one BAM walk."""

    label: str = ""
    n_clips_examined: int = 0  # soft-clips of len >= min_len seen
    n_polya: int = 0           # A-rich 3' clips (poly-A tails)
    n_polyt: int = 0           # T-rich 5' clips (antisense poly-A)
    polya_len_sum: int = 0
    polya_len_max: int = 0
    contigs_with_polya: set[str] = field(default_factory=set)

    @property
    def polya_pct_of_clips(self) -> float:
        return 100.0 * self.n_polya / self.n_clips_examined if self.n_clips_examined else 0.0

    @property
    def polya_mean_len(self) -> float:
        return self.polya_len_sum / self.n_polya if self.n_polya else 0.0


def tally_clip(
    stats: PolyAStats,
    side: str,
    seq: str,
    contig: str = "",
    *,
    min_len: int = POLYA_MIN_LEN,
    min_frac: float = POLYA_MIN_FRAC,
) -> None:
    """Fold one (mRNA-oriented) soft-clip into *stats* in place.

    Designed to be called from inside an existing BAM clip loop (e.g.
    :func:`bam_diagnostic.diagnose_bam`) so the read BAM is walked only once.
    """
    if len(seq) < min_len:
        return
    stats.n_clips_examined += 1
    kind = classify_clip(side, seq, min_len=min_len, min_frac=min_frac)
    if kind == "polyA":
        stats.n_polya += 1
        stats.polya_len_sum += len(seq)
        stats.polya_len_max = max(stats.polya_len_max, len(seq))
        if contig:
            stats.contigs_with_polya.add(contig)
    elif kind == "polyT":
        stats.n_polyt += 1


def characterize_polya_bam(
    bam_path: Path,
    label: str,
    *,
    min_clip_len: int = POLYA_MIN_LEN,
    min_mapq: int = 0,
) -> PolyAStats:
    """Own-pass poly-A characterization of *bam_path* (for the transcript BAM).

    Uses :func:`bam_diagnostic._extract_clips` (mRNA-oriented clips) over primary
    alignments. ``min_mapq`` defaults to 0 so multi-mapping transcripts (segemehl
    ``-H 1`` can assign low MAPQ) are still characterized.
    """
    import pysam

    # Local import: bam_diagnostic imports this module at top level, so importing
    # it back here only inside the function avoids a circular import.
    from eukan.assembly.bam_diagnostic import _extract_clips, _iter_primary_alignments

    stats = PolyAStats(label=label)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in _iter_primary_alignments(bam, min_mapq=min_mapq):
            contig = read.reference_name or ""
            for side, seq, _anchor in _extract_clips(read, min_clip_len):
                tally_clip(stats, side, seq, contig)
    finally:
        bam.close()
    return stats


def scan_fasta_polya(
    fasta_path: Path, *, min_len: int = POLYA_MIN_LEN, min_frac: float = POLYA_MIN_FRAC
) -> tuple[int, int]:
    """Return ``(n_seqs, n_with_polyA_tail)`` for a FASTA (e.g. the unmapped set).

    A poly-A tail is the trailing *min_len* bases being >= *min_frac* A (sense) **or**
    the leading *min_len* bases being >= *min_frac* T (antisense poly-A). Both are
    checked because the de novo transcript library is unstranded, so a transcript may
    be assembled in either orientation and carry its tail as a 3' poly-A or a 5'
    poly-T. A proxy (it does not measure tail length) — enough to flag whether
    failures to map are poly-A-laden.
    """
    from Bio import SeqIO

    n = n_polya = 0
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        n += 1
        seq = str(rec.seq).upper()
        if len(seq) < min_len:
            continue
        tail_polya = seq[-min_len:].count("A") / min_len >= min_frac
        head_polyt = seq[:min_len].count("T") / min_len >= min_frac
        if tail_polya or head_polyt:
            n_polya += 1
    return n, n_polya


def has_section(wd: Path, section: str) -> bool:
    """True if ``polyA_diagnostic.json`` in *wd* already carries *section*.

    Lets a producer skip a redundant backfill pass (e.g. the read-BAM "reads" section
    on a resumed run) without clobbering work another step already wrote.
    """
    path = wd / POLYA_DIAGNOSTIC
    if not path.exists():
        return False
    try:
        return section in json.loads(path.read_text())
    except (ValueError, OSError):
        return False


def stats_to_dict(stats: PolyAStats) -> dict:
    """JSON-serialisable summary of one :class:`PolyAStats` section."""
    return {
        "n_softclips_examined": stats.n_clips_examined,
        "n_polyA_3p": stats.n_polya,
        "n_polyT_5p": stats.n_polyt,
        "polyA_pct_of_softclips": round(stats.polya_pct_of_clips, 4),
        "polyA_mean_len": round(stats.polya_mean_len, 2),
        "polyA_max_len": stats.polya_len_max,
        "n_contigs_with_polyA": len(stats.contigs_with_polya),
    }


def write_polya_section(wd: Path, section: str, payload: dict) -> Path:
    """Merge *payload* under key *section* into ``wd/polyA_diagnostic.json``.

    Multiple steps (read mapping, then transcript mapping) contribute different
    sections to the one file; each call loads, updates its own key, and rewrites,
    so ordering between the producing steps does not matter.
    """
    path = wd / POLYA_DIAGNOSTIC
    data: dict = {}
    if path.exists():
        try:
            data = json.loads(path.read_text())
        except (ValueError, OSError):
            data = {}
    data[section] = payload
    path.write_text(json.dumps(data, indent=2))
    return path
