"""Transcript models as plain records: read, write, split and stitch.

The assembly steps hand transcript models to each other as gene>mRNA>exon GFF3
(and read StringTie's transcript>exon GTF). This module holds the one record type
and the one reader, writer, splitter and sequence stitcher for them. It reads rows
with :mod:`eukan.gff.rows` instead of building a gffutils database, which is about
30x slower on a genome-wide model set.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from Bio.Seq import Seq

from eukan.gff.rows import iter_gff3_rows, parse_attrs

if TYPE_CHECKING:
    from eukan.infra.genome import ContigIndex


@dataclass
class Transcript:
    """A transcript model: id, location, and genomic exon blocks (start-sorted)."""

    tid: str
    chrom: str
    strand: str
    source: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def length(self) -> int:
        """Spliced length: the summed exon lengths."""
        return sum(e - s + 1 for s, e in self.exons)

    def exons_5to3(self) -> list[tuple[int, int]]:
        """Exon blocks in transcript order (descending genomic on ``-``)."""
        return self.exons if self.strand != "-" else list(reversed(self.exons))


def iter_transcripts(path: str | Path) -> Iterator[Transcript]:
    """Stream the transcript models of a GFF3 or GTF file, one at a time.

    A transcript is an ``mRNA``/``transcript`` row (id from ``ID`` or
    ``transcript_id``) plus its ``exon`` rows (``Parent`` or ``transcript_id``);
    the transcript row is optional. Each transcript's rows must be contiguous, as
    eukan, StringTie and combinr write them: rows interleaved with another
    transcript's come out as separate records. Contig, strand and source are taken
    from the transcript's first row, exons are sorted by start, and transcripts
    without exons are dropped.
    """
    cur: Transcript | None = None
    for cols in iter_gff3_rows(path):
        if cols[2] in ("transcript", "mRNA"):
            attrs = parse_attrs(cols[8])
            tid = attrs.get("ID") or attrs.get("transcript_id")
        elif cols[2] == "exon":
            attrs = parse_attrs(cols[8])
            tid = attrs.get("Parent") or attrs.get("transcript_id")
        else:
            continue
        if not tid:
            continue
        if cur is None or cur.tid != tid:
            if cur is not None and cur.exons:
                cur.exons.sort()
                yield cur
            cur = Transcript(tid, cols[0], cols[6], cols[1])
        if cols[2] == "exon":
            cur.exons.append((int(cols[3]), int(cols[4])))
    if cur is not None and cur.exons:
        cur.exons.sort()
        yield cur


def read_transcripts(path: str | Path) -> list[Transcript]:
    """All transcript models of *path*, in file order (see :func:`iter_transcripts`)."""
    return list(iter_transcripts(path))


def write_transcript(fh: TextIO, tx: Transcript) -> None:
    """Write *tx* as gene>mRNA>exon GFF3 rows, under a synthetic ``<tid>.gene``."""
    gid = f"{tx.tid}.gene"
    loc = f"{tx.chrom}\t{tx.source}"
    fh.write(f"{loc}\tgene\t{tx.start}\t{tx.end}\t.\t{tx.strand}\t.\tID={gid}\n")
    fh.write(f"{loc}\tmRNA\t{tx.start}\t{tx.end}\t.\t{tx.strand}\t.\tID={tx.tid};Parent={gid}\n")
    for k, (s, e) in enumerate(tx.exons, start=1):
        fh.write(
            f"{loc}\texon\t{s}\t{e}\t.\t{tx.strand}\t.\tID={tx.tid}.exon{k};Parent={tx.tid}\n"
        )


def write_transcripts(models: Iterable[Transcript], out_gff: str | Path) -> None:
    """Write *models* to *out_gff* as gene>mRNA>exon GFF3, sorted by contig and start."""
    with open(out_gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for tx in sorted(models, key=lambda t: (t.chrom, t.start)):
            write_transcript(fh, tx)


def partition_exons(
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


def split_at(
    tx: Transcript,
    clips: list[int],
    min_segment: int,
    *,
    suffix: str = "j",
    strand_of: Callable[[int, int], str] | None = None,
) -> list[Transcript] | None:
    """Split *tx* after each spliced position in *clips*; ``None`` if nothing to cut.

    Clips outside ``(0, length)`` are ignored. Pieces shorter than *min_segment*
    are dropped; the rest are named ``<tid>.<suffix>1``, ``<tid>.<suffix>2``, ...
    A piece keeps *tx*'s strand unless *strand_of* is given, which is called with
    the piece's 1-based spliced span ``(first, last)`` and returns its strand.
    """
    total_len = tx.length
    valid = [p for p in clips if 0 < p < total_len]
    if not valid:
        return None
    cut_points = [0, *sorted(valid), total_len]
    out: list[Transcript] = []
    for i, blocks in enumerate(partition_exons(tx.exons_5to3(), valid, tx.strand, total_len)):
        if sum(e - s + 1 for s, e in blocks) < min_segment:
            continue
        strand = strand_of(cut_points[i] + 1, cut_points[i + 1]) if strand_of else tx.strand
        out.append(
            Transcript(f"{tx.tid}.{suffix}{len(out) + 1}", tx.chrom, strand, tx.source,
                       sorted(blocks))
        )
    return out or None


def stitch(tx: Transcript, contigs: ContigIndex) -> str:
    """Spliced sequence of *tx*, 5'->3': exons joined, reverse-complemented on ``-``."""
    seq = "".join(str(contigs[tx.chrom][s - 1 : e].seq) for s, e in tx.exons)
    return str(Seq(seq).reverse_complement()) if tx.strand == "-" else seq
