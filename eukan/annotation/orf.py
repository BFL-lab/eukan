"""ORF identification in transcript assemblies.

Finds start/stop codons, selects the longest ORF per transcript,
and maps transcript-space coordinates back to genome coordinates.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

import gffutils
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq

from eukan.gff import parser as gffparser
from eukan.gff.intersecter import merge_fully_overlapping_transcript_genes

# ---------------------------------------------------------------------------
# Sequence extraction from FeatureDB
# ---------------------------------------------------------------------------


def fetch_aligned_sequences(
    gff3db: gffutils.FeatureDB,
    fasta: str | Path,
    featuretype: str = "exon",
) -> list[tuple[str, str, str]]:
    """Extract concatenated exon/CDS sequences for each mRNA.

    Returns:
        List of (mRNA_id, strand, sequence) tuples.
    """
    from eukan.infra.genome import ContigIndex

    result: list[tuple[str, str, str]] = []
    # Sort mRNAs by chromosome so the ContigIndex single-record cache stays warm.
    mrnas = sorted(gff3db.features_of_type("mRNA"), key=lambda m: (m.chrom, m.start))
    with ContigIndex(fasta) as contigs:
        for mRNA in mrnas:
            contig = contigs[mRNA.chrom].seq
            assert contig is not None  # FASTA records always carry a sequence
            seq_parts = [
                str(contig[child.start - 1 : child.end]).upper()
                for child in gff3db.children(mRNA, featuretype=featuretype)
            ]
            seq = "".join(seq_parts)
            if mRNA.strand == "-":
                seq = str(Seq(seq).reverse_complement())
            result.append((mRNA.id, mRNA.strand, seq))

    return result


# ---------------------------------------------------------------------------
# Start/stop codon scanning
# ---------------------------------------------------------------------------


def find_starts_stops(
    seqs: list[tuple[str, str, str]],
    genetic_code: int = 1,
) -> pd.DataFrame:
    """Find all start and stop codon positions in transcript sequences.

    Returns:
        DataFrame with columns: id, seqlen, pos, strand, codon, frame
    """
    code = CodonTable.unambiguous_dna_by_id[genetic_code]
    start_pattern = "(" + "|".join(code.start_codons) + ")"
    stop_pattern = "(" + "|".join(code.stop_codons) + ")"

    rows: list[tuple[str, int, int, str, str, int]] = []

    for name, strand, sequence in seqs:
        seqlen = len(sequence)
        for codon_type, pattern in [("start", start_pattern), ("stop", stop_pattern)]:
            for match in re.finditer(pattern, sequence):
                pos = match.start() + 1
                frame = (pos - 1) % 3
                rows.append((name, seqlen, pos, strand, codon_type, frame))

    return pd.DataFrame(
        rows, columns=["id", "seqlen", "pos", "strand", "codon", "frame"]
    )


# ---------------------------------------------------------------------------
# Longest ORF selection
# ---------------------------------------------------------------------------


_ORF_COLUMNS = ["id", "seqlen", "start", "stop", "strand", "frame", "orflen"]


def fetch_longest_orf(
    codon_df: pd.DataFrame,
    min_orf_len: int = 90,
    orf_len_frac: float = 0.50,
) -> pd.DataFrame:
    """Select the longest ORF per transcript from start/stop codon positions.

    For each (transcript, strand, frame), pairs each start with the
    nearest downstream stop and tracks the longest such pair.  Across
    frames, the longest ORF per (transcript, strand) wins.

    The previous implementation built an O(starts * stops) merged
    DataFrame per frame; this version walks the codons once and
    binary-searches for the nearest stop, which is O((S+T) log T) per
    frame and avoids materialising the cross product.

    Args:
        codon_df: Output of find_starts_stops.
        min_orf_len: Minimum ORF length in nucleotides (strict >).
        orf_len_frac: Minimum ORF length as fraction of transcript length (strict >).

    Returns:
        DataFrame with columns: id, seqlen, start, stop, strand, frame, orflen
    """
    if codon_df.empty:
        return pd.DataFrame(columns=_ORF_COLUMNS)

    from bisect import bisect_right
    from collections import defaultdict
    from dataclasses import dataclass, field

    @dataclass
    class _FrameData:
        seqlen: int = 0
        starts: list[int] = field(default_factory=list)
        stops: list[int] = field(default_factory=list)

    # Group codons by (id, strand, frame). seqlen is invariant per id.
    groups: dict[tuple[str, str, int], _FrameData] = defaultdict(_FrameData)
    for row in codon_df.itertuples(index=False):
        slot = groups[(row.id, row.strand, row.frame)]
        slot.seqlen = row.seqlen
        (slot.starts if row.codon == "start" else slot.stops).append(row.pos)

    # For each frame, find the longest in-frame ORF; keep the per-transcript winner.
    best_per_transcript: dict[tuple[str, str], dict] = {}
    for (tid, strand, frame), data in groups.items():
        stops_sorted = sorted(data.stops)
        if not stops_sorted:
            continue
        seqlen = data.seqlen

        best_orflen = 0
        best_start = best_stop = 0
        for start in data.starts:
            idx = bisect_right(stops_sorted, start)
            if idx == len(stops_sorted):
                continue
            stop = stops_sorted[idx]
            orflen = stop - start
            if orflen > best_orflen:
                best_orflen = orflen
                best_start, best_stop = start, stop

        if best_orflen <= min_orf_len:
            continue
        if seqlen <= 0 or best_orflen / seqlen <= orf_len_frac:
            continue

        key = (tid, strand)
        prev = best_per_transcript.get(key)
        if prev is None or prev["orflen"] < best_orflen:
            best_per_transcript[key] = {
                "id": tid,
                "seqlen": seqlen,
                "start": best_start,
                "stop": best_stop + 2,  # include the full stop codon
                "strand": strand,
                "frame": frame,
                "orflen": best_orflen,
            }

    if not best_per_transcript:
        return pd.DataFrame(columns=_ORF_COLUMNS)
    return pd.DataFrame(list(best_per_transcript.values()), columns=_ORF_COLUMNS)


# ---------------------------------------------------------------------------
# ORF coordinate mapping (transcript → genome)
# ---------------------------------------------------------------------------


def orf_to_genome_coords(
    longest_orfs: pd.DataFrame,
    gff3: gffutils.FeatureDB,
) -> Iterator[gffutils.Feature]:
    """Map transcript-space ORF coordinates to genome CDS features.

    For each mRNA with an identified ORF, walks through its exons and
    creates CDS features at the corresponding genome positions.
    """
    orf_lookup = longest_orfs.set_index("id")
    for mRNA in gff3.features_of_type("mRNA"):
        if mRNA.id not in orf_lookup.index:
            continue

        orf = orf_lookup.loc[mRNA.id]
        orf_start, orf_stop = int(orf["start"]), int(orf["stop"])

        if mRNA.strand == "+":
            yield from _map_orf_plus_strand(gff3, mRNA, orf_start, orf_stop)
        else:
            yield from _map_orf_minus_strand(gff3, mRNA, orf_start, orf_stop)


def _make_cds_from_exon(
    exon: gffutils.Feature, start: int, end: int, phase: int
) -> gffutils.Feature:
    """Create a CDS feature derived from an exon."""
    attrs = dict(exon.attributes)
    exon_id = exon.attributes["ID"][0]
    cds_id = re.sub("exon", "CDS", exon_id)
    if cds_id == exon_id:
        cds_id = f"{exon_id}:CDS"
    attrs["ID"] = [cds_id]
    return gffutils.Feature(
        seqid=exon.chrom,
        source=exon.source,
        featuretype="CDS",
        start=start,
        end=end,
        strand=exon.strand,
        frame=phase,
        attributes=attrs,
    )


def _map_orf_plus_strand(
    gff3: gffutils.FeatureDB,
    mRNA: gffutils.Feature,
    orf_start: int,
    orf_stop: int,
) -> Iterator[gffutils.Feature]:
    """Map ORF to genome coords on the plus strand."""
    exons = list(gff3.children(mRNA, featuretype="exon", order_by="start"))
    cds_count = 0
    next_phase = 0
    consumed = 0  # transcript bases consumed so far

    for exon in exons:
        exon_len = exon.end - exon.start + 1
        exon_transcript_start = consumed
        exon_transcript_end = consumed + exon_len
        consumed = exon_transcript_end

        # Does this exon overlap the ORF in transcript space?
        orf_rel_start = max(orf_start - 1, exon_transcript_start) - exon_transcript_start
        orf_rel_end = min(orf_stop, exon_transcript_end) - exon_transcript_start

        if orf_rel_start >= orf_rel_end or orf_start - 1 >= exon_transcript_end or orf_stop <= exon_transcript_start:
            continue

        cds_count += 1
        genome_start = exon.start + orf_rel_start
        genome_end = exon.start + orf_rel_end - 1

        phase = 0 if cds_count == 1 else next_phase
        cds_len = genome_end - genome_start + 1
        next_phase = (3 - ((cds_len - phase) % 3)) % 3

        yield _make_cds_from_exon(exon, genome_start, genome_end, phase)


def _map_orf_minus_strand(
    gff3: gffutils.FeatureDB,
    mRNA: gffutils.Feature,
    orf_start: int,
    orf_stop: int,
) -> Iterator[gffutils.Feature]:
    """Map ORF to genome coords on the minus strand."""
    exons = list(
        gff3.children(mRNA, featuretype="exon", order_by="start", reverse=True)
    )
    cds_count = 0
    next_phase = 0
    consumed = 0

    # On minus strand, transcript coords go from 3'->5' (high->low genome coord)
    orf_start_adj = orf_start - 1
    orf_stop_adj = orf_stop - 1

    for exon in exons:
        exon_len = exon.end - exon.start + 1
        exon_transcript_start = consumed
        exon_transcript_end = consumed + exon_len
        consumed = exon_transcript_end

        orf_rel_start = max(orf_start_adj, exon_transcript_start) - exon_transcript_start
        orf_rel_end = min(orf_stop_adj, exon_transcript_end) - exon_transcript_start

        if orf_rel_start >= orf_rel_end or orf_start_adj >= exon_transcript_end or orf_stop_adj <= exon_transcript_start:
            continue

        cds_count += 1
        # Minus strand: genome coords are reversed
        genome_end = exon.end - orf_rel_start
        genome_start = exon.end - orf_rel_end + 1

        phase = 0 if cds_count == 1 else next_phase
        cds_len = genome_end - genome_start + 1
        next_phase = (3 - ((cds_len - phase) % 3)) % 3

        yield _make_cds_from_exon(exon, genome_start, genome_end, phase)


# ---------------------------------------------------------------------------
# Transcriptome ORF database creation
# ---------------------------------------------------------------------------


def create_transcriptome_orf_db(
    gff: str | Path, genome: str | Path, genetic_code: int = 1,
) -> gffutils.FeatureDB:
    """Build a FeatureDB of transcript ORFs mapped to genome coordinates.

    1. Parse transcript GFF into a normalized FeatureDB.
    2. Extract sequences, find start/stop codons.
    3. Select longest ORFs.
    4. Map ORF coords back to the genome.
    5. Merge overlapping transcript genes.
    """
    dialect = gffutils.DataIterator(str(gff)).dialect
    dialect["fmt"] = "gtf"
    xcripts = gffutils.create_db(
        str(gff), ":memory:", dialect=dialect,
        gtf_transcript_key="Parent", gtf_gene_key="Parent",
    )
    dialect["fmt"] = "gff3"
    xcripts.dialect = dialect
    xcripts = gffutils.create_db(
        xcripts, ":memory:", dialect=dialect, transform=gffparser.gff3_it,
    )

    seqs = fetch_aligned_sequences(xcripts, genome)
    starts_stops = find_starts_stops(seqs, genetic_code=genetic_code)
    longest_orfs = fetch_longest_orf(starts_stops)

    xcripts.update(
        orf_to_genome_coords(longest_orfs, xcripts),
        merge_strategy="create_unique",
    )
    return merge_fully_overlapping_transcript_genes(xcripts)
