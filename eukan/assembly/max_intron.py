"""Max-intron split: break transcript models at over-long introns.

A strand-agnostic sanitization pass — independent of spliced-leader detection —
that splits any transcript model carrying an intron longer than
``max_intron_len`` into separate genes. It exists because the de novo **segemehl**
transcript→genome path has no native intron bound (STAR/STARlong cap introns at
mapping time via ``--alignIntronMax``), so a segemehl-mapped contig can bridge two
distant loci across one implausibly long gap. Splitting that gap keeps the loci
separate before consensus.

Runs before the genomic SL cut (:mod:`eukan.assembly.sl_cut`): for each mapped
Trinity track it reads the latest model variant (defuse > stranded > raw, via
:func:`eukan.assembly.tracks.resolve_model_source`) and writes
``{stem}.maxintron.gff3`` — **always**, copying the model through unchanged when
nothing needs cutting, so the file is a stable input the SL cut reads directly.
Reuses the shared exon-split machinery (:mod:`eukan.assembly.jaccard`) and streams
one model at a time so a genome-wide model set stays bounded. ``max_intron_len <= 0``
disables the split (a pure copy-through).
"""

from __future__ import annotations

from pathlib import Path

from eukan.assembly.jaccard import (
    _iter_transcript_models,
    _split_transcript,
    _Tx,
    _write_tx,
)
from eukan.assembly.tracks import mapped_transcript_stems, resolve_model_source
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def _long_intron_cut_offsets(
    exons_5to3: list[tuple[int, int]], strand: str, max_intron_len: int
) -> set[int]:
    """Spliced "cut after base P" offsets severing every intron > *max_intron_len*.

    *exons_5to3* are 1-based inclusive genomic blocks in 5'->3' order (ascending for
    ``+``/``.``, descending for ``-``) — the order :func:`jaccard._split_transcript`
    and :func:`jaccard._partition_exons` expect. Each gap between consecutive blocks
    is an intron; when its genomic length exceeds the limit, cut at the cumulative
    spliced end of the 5'-side exon, which severs the model at that intron. Returns
    offsets in the same spliced space as SL clips so the two unite cleanly.
    ``max_intron_len <= 0`` disables the cut.
    """
    if max_intron_len <= 0 or len(exons_5to3) < 2:
        return set()
    cuts: set[int] = set()
    consumed = 0
    for i, (gstart, gend) in enumerate(exons_5to3):
        consumed += gend - gstart + 1
        if i + 1 < len(exons_5to3):
            ns, ne = exons_5to3[i + 1]
            intron_len = (ns - gend - 1) if strand != "-" else (gstart - ne - 1)
            if intron_len > max_intron_len:
                cuts.add(consumed)
    return cuts


def _count_long_introns(tx: _Tx, max_intron_len: int) -> int:
    """Number of introns in *tx* longer than *max_intron_len* (genomic gaps)."""
    if max_intron_len <= 0:
        return 0
    ex = tx.exons  # start-sorted ascending
    return sum(
        1 for i in range(len(ex) - 1) if ex[i + 1][0] - ex[i][1] - 1 > max_intron_len
    )


def _cut_one(tx: _Tx, min_segment: int, *, max_intron_len: int) -> list[_Tx] | None:
    """Split *tx* at every over-long intron; ``None`` if nothing to cut.

    Strand-agnostic: a genomic gap is severed regardless of strand (a ``.``-strand
    transcript is treated as ``+``/ascending). The cut points are computed in the
    same spliced space :func:`jaccard._split_transcript` uses, which also drops
    sub-*min_segment* fragments and returns ``None`` when nothing valid remains.
    """
    exons_5to3 = tx.exons if tx.strand != "-" else list(reversed(tx.exons))
    clips = _long_intron_cut_offsets(exons_5to3, tx.strand, max_intron_len)
    if not clips:
        return None
    return _split_transcript(tx, sorted(clips), min_segment)


def cut_models_at_max_intron(
    gff_or_gtf: str | Path,
    out_gff: str | Path,
    *,
    min_segment: int,
    max_intron_len: int,
) -> tuple[int, int]:
    """Split every transcript in *gff_or_gtf* at over-long introns into *out_gff*.

    Streams one transcript at a time. Transcripts with no over-long intron pass
    through unchanged; *out_gff* is always written (a copy-through when nothing is
    cut) so it can serve as the SL cut's stable input. Returns
    ``(transcripts_cut, over_long_introns_severed)``.
    """
    n_cut = 0
    n_long = 0
    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for tx in _iter_transcript_models(gff_or_gtf):
            long_here = _count_long_introns(tx, max_intron_len)
            result = _cut_one(tx, min_segment, max_intron_len=max_intron_len)
            if result is None:
                _write_tx(fh, tx)
            else:
                for piece in result:
                    _write_tx(fh, piece)
                n_cut += 1
                n_long += long_here
    return n_cut, n_long


def run_max_intron_split(config: AssemblyConfig) -> None:
    """Split each Trinity track's models at over-long introns → ``{stem}.maxintron.gff3``.

    Reads the latest model variant per track (defuse > stranded > raw) and writes
    one ``{stem}.maxintron.gff3`` per track that produced models. This is the single
    place the ``max_intron_len`` limit is hard-imposed on transcript models (the de
    novo segemehl path has no native intron bound); the resume fingerprint folds in
    ``max_intron_len`` (see :mod:`eukan.assembly.pipeline`), so tightening ``-M``
    re-runs this step and cascades to the SL cut and combinr.
    """
    wd = config.work_dir
    min_segment = config.min_sl_fragment
    for stem in mapped_transcript_stems():
        src = resolve_model_source(wd, stem)
        if src is None:
            continue
        out = wd / f"{stem}.maxintron.gff3"
        n_cut, n_long = cut_models_at_max_intron(
            src, out, min_segment=min_segment, max_intron_len=config.max_intron_len,
        )
        log.info(
            "Max-intron split %s -> %s (%d transcripts cut, %d over-long introns severed).",
            src.name, out.name, n_cut, n_long,
        )
