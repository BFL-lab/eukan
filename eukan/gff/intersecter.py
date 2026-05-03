"""Genomic interval operations using gffutils.

Handles merging overlapping genes, finding concordant models between
prediction sources, extracting training sets, and non-overlapping gene
detection — all via direct gffutils FeatureDB queries without pybedtools.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import gffutils

from eukan.gff import create_gff_db
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger


def _empty_db() -> gffutils.FeatureDB:
    """Create an empty FeatureDB.

    Workaround for gffutils raising EmptyInputError on empty input.
    Creates a DB with a dummy feature, then deletes it.
    """
    db = gffutils.create_db(
        "chr1\t.\tgene\t1\t1\t.\t+\t.\tID=_empty_placeholder",
        ":memory:", from_string=True,
    )
    db.delete(db["_empty_placeholder"])
    return db

log = get_logger(__name__)

BOUNDARY_TOLERANCE = 3  # bp tolerance for intron-exon boundary matching
MRNA_TOLERANCE = 0.05   # 5% tolerance for mRNA span comparison


# ---------------------------------------------------------------------------
# Overlap detection (pure gffutils)
# ---------------------------------------------------------------------------


def _features_overlap(a: gffutils.Feature, b: gffutils.Feature) -> bool:
    """Check if two features overlap on the same strand."""
    return (
        a.chrom == b.chrom
        and a.strand == b.strand
        and a.start <= b.end
        and b.start <= a.end
    )


def _find_overlapping_genes(
    db: gffutils.FeatureDB,
) -> list[list[gffutils.Feature]]:
    """Find clusters of overlapping genes on the same strand.

    Uses a sweep-line approach: sort genes by (chrom, strand, start),
    then merge overlapping intervals.
    """
    genes = sorted(
        db.features_of_type("gene"),
        key=lambda g: (g.chrom, g.strand, g.start),
    )

    clusters: list[list[gffutils.Feature]] = []
    if not genes:
        return clusters

    current_cluster = [genes[0]]
    cluster_end = genes[0].end

    for gene in genes[1:]:
        if (
            gene.chrom == current_cluster[0].chrom
            and gene.strand == current_cluster[0].strand
            and gene.start <= cluster_end
        ):
            current_cluster.append(gene)
            cluster_end = max(cluster_end, gene.end)
        else:
            if len(current_cluster) > 1:
                clusters.append(current_cluster)
            current_cluster = [gene]
            cluster_end = gene.end

    if len(current_cluster) > 1:
        clusters.append(current_cluster)

    return clusters


# ---------------------------------------------------------------------------
# Overlapping gene consolidation
# ---------------------------------------------------------------------------


def merge_fully_overlapping_transcript_genes(
    gff3db: gffutils.FeatureDB,
) -> gffutils.FeatureDB:
    """Merge genes that fully overlap on the same strand.

    For each cluster of overlapping genes, keeps the first (longest span
    after sorting) as canonical and re-parents all mRNAs from shorter
    genes to the canonical one.
    """
    clusters = _find_overlapping_genes(gff3db)

    # Build replacement map: duplicate gene ID → canonical gene ID
    replacement_ids: dict[str, str] = {}
    for cluster in clusters:
        canonical = cluster[0]
        for dup in cluster[1:]:
            replacement_ids[dup.id] = canonical.id

    if not replacement_ids:
        return gff3db

    log.debug("Merging %d duplicate genes into canonical parents", len(replacement_ids))

    def _consolidate() -> Iterator[gffutils.Feature]:
        for f in gff3db.all_features():
            if f.featuretype == "gene" and f.id in replacement_ids:
                continue  # drop duplicate gene features
            if f.featuretype == "mRNA":
                parent = f.attributes["Parent"][0]
                if parent in replacement_ids:
                    f.attributes["Parent"] = [replacement_ids[parent]]
            yield f

    return create_gff_db(_consolidate())


# ---------------------------------------------------------------------------
# Non-overlapping gene detection
# ---------------------------------------------------------------------------


def find_nonoverlapping_genes(
    db_source: gffutils.FeatureDB,
    db_target: gffutils.FeatureDB,
) -> list[gffutils.Feature]:
    """Find ORF-containing genes in db_source that don't overlap any gene in db_target.

    Returns all features (gene + children) for non-overlapping, ORF-containing genes.
    """
    # Build a per-(chrom, strand) index: intervals sorted by start, with a
    # parallel monotone prefix-max of ends. The prefix-max lets us bisect
    # for the left bound of overlap candidates: any target with
    # max_end_to[i] < query.start cannot overlap.
    target_starts: dict[tuple[str, str], list[int]] = defaultdict(list)
    target_ends: dict[tuple[str, str], list[int]] = defaultdict(list)
    target_max_end_to: dict[tuple[str, str], list[int]] = {}
    for gene in db_target.features_of_type("gene"):
        key = (gene.chrom, gene.strand)
        target_starts[key].append(gene.start)
        target_ends[key].append(gene.end)

    for key in target_starts:
        paired = sorted(zip(target_starts[key], target_ends[key], strict=True))
        target_starts[key] = [s for s, _ in paired]
        ends_sorted = [e for _, e in paired]
        target_ends[key] = ends_sorted
        max_end_to: list[int] = []
        running_max = 0
        for end in ends_sorted:
            if end > running_max:
                running_max = end
            max_end_to.append(running_max)
        target_max_end_to[key] = max_end_to

    result: list[gffutils.Feature] = []

    for gene in db_source.features_of_type("gene"):
        key = (gene.chrom, gene.strand)
        starts = target_starts.get(key)
        if not starts:
            overlaps = False
        else:
            ends = target_ends[key]
            max_end_to = target_max_end_to[key]
            idx_right = bisect.bisect_right(starts, gene.end)
            if idx_right == 0:
                overlaps = False
            else:
                idx_left = bisect.bisect_left(max_end_to, gene.start)
                overlaps = any(
                    ends[i] >= gene.start for i in range(idx_left, idx_right)
                )
        if overlaps:
            continue

        # Check that this gene has CDS children (i.e., contains an ORF)
        child_types = {f.featuretype for f in db_source.children(gene)}
        if "CDS" not in child_types and not any(
            "CDS" in {gc.featuretype for gc in db_source.children(c)}
            for c in db_source.children(gene, featuretype="mRNA")
        ):
            continue

        result.append(gene)
        result.extend(db_source.children(gene, order_by="featuretype", reverse=True))

    log.debug("Found %d non-overlapping ORF-containing genes", len([f for f in result if f.featuretype == "gene"]))
    return result


# ---------------------------------------------------------------------------
# Concordant model detection
# ---------------------------------------------------------------------------


def _get_cds_list(
    db: gffutils.FeatureDB, mrna: gffutils.Feature
) -> list[gffutils.Feature]:
    """Get sorted CDS features for an mRNA."""
    return list(db.children(mrna, featuretype="CDS", order_by="start"))


def _cds_boundaries_match(
    cds_a: list[gffutils.Feature],
    cds_b: list[gffutils.Feature],
) -> bool:
    """Check if two CDS lists have concordant intron-exon boundaries.

    For single-CDS genes: always match (no introns to compare).
    For multi-CDS genes: internal boundaries must match within BOUNDARY_TOLERANCE bp.
    Terminal CDS boundaries have a looser check (only the intron-facing end).
    """
    n = len(cds_a)
    if n == 1:
        return True

    for i in range(n):
        a, b = cds_a[i], cds_b[i]

        if i == 0:
            # First CDS: check 3' boundary (end)
            if abs(a.end - b.end) > BOUNDARY_TOLERANCE:
                return False
        elif i == n - 1:
            # Last CDS: check 5' boundary (start)
            if abs(a.start - b.start) > BOUNDARY_TOLERANCE:
                return False
        else:
            # Internal CDS: both boundaries must match
            if abs(a.start - b.start) > BOUNDARY_TOLERANCE:
                return False
            if abs(a.end - b.end) > BOUNDARY_TOLERANCE:
                return False

    return True


def _mrna_spans_match(
    mrna_a: gffutils.Feature, mrna_b: gffutils.Feature
) -> bool:
    """Check if two mRNA features have similar spans (within 5% tolerance).

    Tolerance is computed from the longer gene span so that genes near
    coordinate zero are not penalised by a near-zero absolute tolerance.
    """
    span = max(mrna_a.end - mrna_a.start, mrna_b.end - mrna_b.start, 1)
    tol = span * MRNA_TOLERANCE
    return (
        abs(mrna_a.start - mrna_b.start) <= tol
        and abs(mrna_a.end - mrna_b.end) <= tol
    )


def find_concordant_models(
    gff3_1: str | Path, gff3_2: str | Path
) -> gffutils.FeatureDB:
    """Find structurally concordant gene models between two GFF3 files.

    Two models are concordant if:
    1. Their mRNAs overlap on the same strand with similar spans
    2. They have the same number of CDS features
    3. Their intron-exon boundaries match within BOUNDARY_TOLERANCE bp

    Returns a FeatureDB containing only the concordant models from gff3_1.
    """
    db1 = create_gff_db(gff3_1)
    db2 = create_gff_db(gff3_2)

    # Index db2 mRNAs by (chrom, strand) -- sorted by start with a parallel
    # monotone prefix-max-end for bounded overlap lookup.
    db2_buckets: dict[tuple[str, str], list[gffutils.Feature]] = defaultdict(list)
    for mrna in db2.features_of_type("mRNA"):
        db2_buckets[(mrna.chrom, mrna.strand)].append(mrna)
    db2_index: dict[
        tuple[str, str],
        tuple[list[gffutils.Feature], list[int], list[int]],
    ] = {}
    for key, mrnas in db2_buckets.items():
        mrnas.sort(key=lambda m: m.start)
        starts = [m.start for m in mrnas]
        max_end_to: list[int] = []
        running_max = 0
        for m in mrnas:
            if m.end > running_max:
                running_max = m.end
            max_end_to.append(running_max)
        db2_index[key] = (mrnas, starts, max_end_to)

    concordant_gene_ids: set[str] = set()

    for mrna1 in db1.features_of_type("mRNA"):
        cds1 = _get_cds_list(db1, mrna1)
        if not cds1:
            continue

        key = (mrna1.chrom, mrna1.strand)
        entry = db2_index.get(key)
        if entry is None:
            continue
        mrnas2, starts2, max_end_to2 = entry
        idx_right = bisect.bisect_right(starts2, mrna1.end)
        if idx_right == 0:
            continue
        idx_left = bisect.bisect_left(max_end_to2, mrna1.start)

        for i in range(idx_left, idx_right):
            mrna2 = mrnas2[i]
            if mrna2.end < mrna1.start:
                continue

            if not _mrna_spans_match(mrna1, mrna2):
                continue

            cds2 = _get_cds_list(db2, mrna2)
            if len(cds1) != len(cds2):
                continue

            if _cds_boundaries_match(cds1, cds2):
                parent_id = mrna1.attributes["Parent"][0]
                concordant_gene_ids.add(parent_id)
                break  # found a match for this mRNA, move on

    log.debug("Found %d concordant gene models", len(concordant_gene_ids))

    # Collect all features for concordant genes
    if not concordant_gene_ids:
        return _empty_db()

    features: list[gffutils.Feature] = []
    for gene_id in concordant_gene_ids:
        try:
            features.append(db1[gene_id])
            features.extend(db1.children(gene_id))
        except gffutils.FeatureNotFoundError:
            continue

    return create_gff_db(features) if features else _empty_db()


# ---------------------------------------------------------------------------
# Non-redundant model combination
# ---------------------------------------------------------------------------


def combine_nonredundant_models(
    *feature_dbs: gffutils.FeatureDB,
) -> gffutils.FeatureDB:
    """Combine gene models from multiple FeatureDBs, removing redundancy.

    Collects all unique gene IDs across databases. For IDs that appear in
    multiple databases, the first database takes precedence.
    """
    if len(feature_dbs) not in (2, 3):
        raise ValueError(f"Expected 2 or 3 FeatureDBs, got {len(feature_dbs)}")

    seen_ids: set[str] = set()
    all_features: list[gffutils.Feature] = []

    for db in feature_dbs:
        for gene in db.features_of_type("gene"):
            if gene.id in seen_ids:
                continue
            seen_ids.add(gene.id)
            all_features.append(gene)
            all_features.extend(db.children(gene))

    log.debug("Combined %d non-redundant gene models from %d sources", len(seen_ids), len(feature_dbs))
    return create_gff_db(all_features) if all_features else _empty_db()


# ---------------------------------------------------------------------------
# Training set extraction
# ---------------------------------------------------------------------------


def extract_supported_models(
    *gff3_paths: str | Path, output_dir: Path | None = None,
) -> gffutils.FeatureDB:
    """Extract concordant gene models supported by multiple evidence sources.

    Args:
        gff3_paths: 2 or 3 GFF3 file paths to compare.
        output_dir: Directory to write training_set.gff3. Defaults to cwd.

    Returns:
        FeatureDB of concordant training set models.
    """
    paths = list(gff3_paths)
    if len(paths) == 2:
        training_set = find_concordant_models(paths[0], paths[1])
    elif len(paths) == 3:
        pair1 = find_concordant_models(paths[0], paths[1])
        pair2 = find_concordant_models(paths[0], paths[2])
        pair3 = find_concordant_models(paths[2], paths[1])
        training_set = combine_nonredundant_models(pair1, pair2, pair3)
    else:
        raise ValueError(f"Expected 2 or 3 paths, got {len(paths)}")

    out_path = (output_dir / "training_set.gff3") if output_dir else Path("training_set.gff3")
    featuredb2gff3_file(training_set, out_path)
    return create_gff_db(out_path)
