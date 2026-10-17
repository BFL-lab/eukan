"""ORF identification in transcript assemblies.

Finds start/stop codons, selects the longest ORF per transcript,
and maps transcript-space coordinates back to genome coordinates.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

import gffutils
from Bio.Data import CodonTable

from eukan.gff import create_gff_db
from eukan.gff.concordance import merge_fully_overlapping_transcript_genes
from eukan.gff.io import iter_assembled_sequences
from eukan.gff.transforms import derived_feature, gff3_it, swap_id_kind

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
        List of (mRNA_id, strand, sequence) tuples. Sequences are uppercased
        so soft-masked input doesn't break codon regex matching downstream.
    """
    return [
        (mrna.id, mrna.strand, seq.upper())
        for mrna, seq in iter_assembled_sequences(
            gff3db, fasta, child_featuretype=featuretype,
        )
    ]


# ---------------------------------------------------------------------------
# Longest ORF selection
# ---------------------------------------------------------------------------


def find_longest_orfs(
    seqs: list[tuple[str, str, str]],
    genetic_code: int = 1,
    min_orf_len: int = 90,
    orf_len_frac: float = 0.50,
) -> dict[str, tuple[int, int]]:
    """Longest ORF of each transcript, as 1-based ``(start, end)``.

    Scans the three frames of each sequence (already on the transcript's
    strand). In a frame, the ORF ending at a stop codon starts at the first
    start codon after the previous stop; the longest ORF wins, the upstream
    one on a tie. Across frames a tie goes to the frame whose first start
    codon comes first. The ORF length, start codon to stop codon with the
    stop excluded, must exceed *min_orf_len* and *orf_len_frac* of the
    sequence length. ``end`` is the last base of the stop codon.

    Args:
        seqs: ``(transcript_id, strand, sequence)`` from fetch_aligned_sequences.
        genetic_code: NCBI genetic code for the start and stop codon sets.
        min_orf_len: Minimum ORF length in nucleotides (strict >).
        orf_len_frac: Minimum ORF length as fraction of transcript length (strict >).
    """
    code = CodonTable.unambiguous_dna_by_id[genetic_code]
    # Zero-width lookahead so overlapping codons are all found: a plain
    # re.finditer resumes after each match, and several codes' start (or
    # stop) sets overlap, e.g. ATT/TTG in code 11 (ATTG), TAG/AGA in code 2.
    start_re = re.compile("(?=(" + "|".join(code.start_codons) + "))")
    stop_re = re.compile("(?=(" + "|".join(code.stop_codons) + "))")

    orfs: dict[str, tuple[int, int]] = {}
    for name, _strand, seq in seqs:
        starts = [m.start() for m in start_re.finditer(seq)]
        stops = [m.start() for m in stop_re.finditer(seq)]
        # best = ((length, -first start of its frame), start, stop), 0-based
        best: tuple[tuple[int, int], int, int] | None = None
        for frame in range(3):
            frame_starts = [p for p in starts if p % 3 == frame]
            i = 0  # first start codon after the previous stop
            for stop in (p for p in stops if p % 3 == frame):
                if i == len(frame_starts):
                    break
                if frame_starts[i] < stop:
                    rank = (stop - frame_starts[i], -frame_starts[0])
                    if best is None or rank > best[0]:
                        best = (rank, frame_starts[i], stop)
                    while i < len(frame_starts) and frame_starts[i] < stop:
                        i += 1
        if best is None:
            continue
        length = best[0][0]
        if length > min_orf_len and length / len(seq) > orf_len_frac:
            orfs[name] = (best[1] + 1, best[2] + 3)
    return orfs


# ---------------------------------------------------------------------------
# ORF coordinate mapping (transcript → genome)
# ---------------------------------------------------------------------------


def orf_to_genome_coords(
    longest_orfs: dict[str, tuple[int, int]],
    gff3: gffutils.FeatureDB,
) -> Iterator[gffutils.Feature]:
    """Map transcript-space ORF coordinates to genome CDS features.

    For each mRNA with an identified ORF, walks through its exons and
    creates CDS features at the corresponding genome positions.
    """
    for mRNA in gff3.features_of_type("mRNA"):
        if mRNA.id not in longest_orfs:
            continue

        orf_start, orf_stop = longest_orfs[mRNA.id]

        if mRNA.strand == "+":
            yield from _map_orf_plus_strand(gff3, mRNA, orf_start, orf_stop)
        else:
            yield from _map_orf_minus_strand(gff3, mRNA, orf_start, orf_stop)


def _make_cds_from_exon(
    exon: gffutils.Feature, start: int, end: int, phase: int
) -> gffutils.Feature:
    """Create a CDS feature derived from an exon."""
    attrs = dict(exon.attributes)
    attrs["ID"] = [swap_id_kind(exon.attributes["ID"][0], "exon", "CDS")]
    return derived_feature(
        exon, "CDS", attrs, frame=phase, start=start, end=end,
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
    2. Extract sequences and find each one's longest ORF.
    3. Map ORF coords back to the genome.
    4. Merge overlapping transcript genes.
    """
    dialect = gffutils.DataIterator(str(gff)).dialect
    dialect["fmt"] = "gtf"
    xcripts = create_gff_db(
        gff, dialect=dialect,
        gtf_transcript_key="Parent", gtf_gene_key="Parent",
    )
    dialect["fmt"] = "gff3"
    xcripts.dialect = dialect
    xcripts = create_gff_db(xcripts, dialect=dialect, transform=gff3_it)

    seqs = fetch_aligned_sequences(xcripts, genome)
    longest_orfs = find_longest_orfs(seqs, genetic_code=genetic_code)

    xcripts.update(
        orf_to_genome_coords(longest_orfs, xcripts),
        merge_strategy="create_unique",
    )
    return merge_fully_overlapping_transcript_genes(xcripts)
