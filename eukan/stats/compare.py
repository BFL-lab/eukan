"""Annotation quality assessment: compare predicted GFF3 against a
reference/previous annotation.

Computes strand-aware overlaps at gene, mRNA, CDS, and intron levels.

Gene-level classification:
  exact       — prediction coordinates identical to reference
  inexact     — prediction overlaps a single reference with boundary differences
  missing     — reference gene with no overlapping prediction (FN)
  merged      — prediction spans 2+ references (>50% of each ref covered)
  fragmented  — 2+ predictions each overlap >=50% of a single reference

Transcript (mRNA) level:
  Matched gene pairs -> maximum pairwise nucleotide overlap matching.
  Unmatched reference mRNAs -> missing; surplus predictions -> FP.

CDS and intron levels:
  Within matched mRNA pairs -> pairwise overlap matching.
  Unmatched reference features -> missing; surplus predictions -> FP.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from itertools import permutations
from pathlib import Path

import gffutils

from eukan.stats.models import (
    FRAG_THRESHOLD,
    MERGE_THRESHOLD,
    ComparisonResult,
    FeatureRecord,
    Interval,
    gene_stats_from_records,
    subfeature_stats_from_records,
)


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------


def _load_gff3(path: Path) -> gffutils.FeatureDB:
    return gffutils.create_db(
        str(path), ":memory:", merge_strategy="create_unique",
        sort_attribute_values=True,
    )


def _parse_features(
    db: gffutils.FeatureDB, feat_type: str, parent_attr: bool = False,
) -> list[Interval]:
    intervals = []
    for f in db.features_of_type(feat_type):
        parent_id = None
        if parent_attr and "Parent" in f.attributes:
            parent_id = f.attributes["Parent"][0]
        intervals.append(Interval(
            chrom=f.chrom, strand=f.strand,
            start=f.start, end=f.end,
            feat_id=f.id, parent_id=parent_id,
        ))
    return intervals


def _derive_introns(cds_list: list[Interval]) -> list[Interval]:
    """Derive intron intervals from sorted CDS features of one mRNA."""
    if len(cds_list) < 2:
        return []
    sorted_cds = sorted(cds_list, key=lambda c: c.start)
    return [
        Interval(
            chrom=sorted_cds[i].chrom,
            strand=sorted_cds[i].strand,
            start=sorted_cds[i].end + 1,
            end=sorted_cds[i + 1].start - 1,
            feat_id=f"{sorted_cds[i].parent_id}:intron:{i + 1}",
            parent_id=sorted_cds[i].parent_id,
        )
        for i in range(len(sorted_cds) - 1)
    ]


def _group_by_parent(
    intervals: list[Interval],
) -> dict[str, list[Interval]]:
    """Group intervals by their parent_id."""
    groups: dict[str, list[Interval]] = defaultdict(list)
    for iv in intervals:
        groups[iv.parent_id].append(iv)
    return groups


# ---------------------------------------------------------------------------
# Overlap computation
# ---------------------------------------------------------------------------


# Spatial index: (chrom, strand) -> (sorted intervals, starts, prefix-max ends)
_Index = dict[tuple[str, str], tuple[list[Interval], list[int], list[int]]]


def _build_index(intervals: list[Interval]) -> _Index:
    """Build a strand-aware interval index with monotone prefix-max ends.

    ``max_end_to[i]`` is ``max(end[0..i])``.  Because intervals are sorted
    by start, this array is monotone non-decreasing, so a binary search on
    it gives the leftmost candidate whose end could reach a query.start.
    """
    buckets: dict[tuple[str, str], list[Interval]] = defaultdict(list)
    for iv in intervals:
        buckets[(iv.chrom, iv.strand)].append(iv)
    idx: _Index = {}
    for key, ivs in buckets.items():
        ivs.sort(key=lambda x: x.start)
        starts = [iv.start for iv in ivs]
        max_end_to: list[int] = []
        running_max = 0
        for iv in ivs:
            if iv.end > running_max:
                running_max = iv.end
            max_end_to.append(running_max)
        idx[key] = (ivs, starts, max_end_to)
    return idx


def _find_overlaps(
    query: Interval,
    targets_idx: _Index,
) -> list[tuple[Interval, int]]:
    """Find all overlapping targets with overlap bp. Strand-aware.

    Bounds the candidate range with two binary searches:
    * right bound: first target with ``start > query.end`` (no future overlap)
    * left bound: first target where ``max_end_to >= query.start``
      (all earlier targets strictly end before the query)
    """
    entry = targets_idx.get((query.chrom, query.strand))
    if entry is None:
        return []
    ivs, starts, max_end_to = entry
    idx_right = bisect_right(starts, query.end)
    if idx_right == 0:
        return []
    idx_left = bisect_left(max_end_to, query.start)
    results = []
    for i in range(idx_left, idx_right):
        t = ivs[i]
        if t.end < query.start:
            continue
        ovl = min(query.end, t.end) - max(query.start, t.start) + 1
        results.append((t, ovl))
    return results


def _overlap_bp(a: Interval, b: Interval) -> int:
    if a.chrom != b.chrom or a.strand != b.strand:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


# ---------------------------------------------------------------------------
# Maximum pairwise matching (for mRNA/CDS/intron within matched parents)
# ---------------------------------------------------------------------------


def _max_pairwise_matching(
    ref_feats: list[Interval],
    pred_feats: list[Interval],
) -> list[tuple[Interval, Interval, int]]:
    """Match ref features to pred features by maximum pairwise overlap sum.

    For small sets (<=8), tries all permutations for optimal matching.
    For larger sets, uses greedy matching by descending overlap.
    """
    if not ref_feats or not pred_feats:
        return []

    # Compute overlap matrix
    overlaps: list[tuple[int, int, int]] = []  # (ref_idx, pred_idx, ovl_bp)
    for ri, r in enumerate(ref_feats):
        for pi, p in enumerate(pred_feats):
            ovl = _overlap_bp(r, p)
            if ovl > 0:
                overlaps.append((ri, pi, ovl))

    if not overlaps:
        return []

    n_ref, n_pred = len(ref_feats), len(pred_feats)

    # For small problems, try optimal matching via permutation
    if n_ref <= 8 and n_pred <= 8:
        best_matches: list[tuple[int, int, int]] = []
        best_sum = 0

        ovl_lookup = {(ri, pi): ovl for ri, pi, ovl in overlaps}

        k = min(n_ref, n_pred)
        ref_indices = list(range(n_ref))
        pred_indices = list(range(n_pred))

        if n_ref <= n_pred:
            for perm in permutations(pred_indices, k):
                total = 0
                matches = []
                for ri, pi in zip(ref_indices, perm):
                    ovl = ovl_lookup.get((ri, pi), 0)
                    total += ovl
                    if ovl > 0:
                        matches.append((ri, pi, ovl))
                if total > best_sum:
                    best_sum = total
                    best_matches = matches
        else:
            for perm in permutations(ref_indices, k):
                total = 0
                matches = []
                for ri, pi in zip(perm, pred_indices):
                    ovl = ovl_lookup.get((ri, pi), 0)
                    total += ovl
                    if ovl > 0:
                        matches.append((ri, pi, ovl))
                if total > best_sum:
                    best_sum = total
                    best_matches = matches

        return [
            (ref_feats[ri], pred_feats[pi], ovl)
            for ri, pi, ovl in best_matches
        ]

    # Greedy fallback for larger sets
    overlaps.sort(key=lambda x: x[2], reverse=True)
    used_refs: set[int] = set()
    used_preds: set[int] = set()
    matches = []
    for ri, pi, ovl in overlaps:
        if ri not in used_refs and pi not in used_preds:
            matches.append((ref_feats[ri], pred_feats[pi], ovl))
            used_refs.add(ri)
            used_preds.add(pi)
    return matches


# ---------------------------------------------------------------------------
# Gene-level classification
# ---------------------------------------------------------------------------


def _classify_genes(
    ref_genes: list[Interval],
    pred_genes: list[Interval],
) -> tuple[dict[str, str], list[FeatureRecord]]:
    """Classify genes and return match map + per-feature records."""
    pred_idx = _build_index(pred_genes)
    ref_idx = _build_index(ref_genes)
    records: list[FeatureRecord] = []

    # --- Step 1: detect merged predictions ---
    # A prediction that covers >50% of 2+ reference genes
    merged_ref_ids: set[str] = set()
    merging_pred_ids: set[str] = set()

    for pred in pred_genes:
        ref_overlaps = _find_overlaps(pred, ref_idx)
        if len(ref_overlaps) < 2:
            continue
        covered_refs = [
            ref_iv.feat_id
            for ref_iv, ovl_bp in ref_overlaps
            if ovl_bp / (ref_iv.end - ref_iv.start + 1) > MERGE_THRESHOLD
        ]
        if len(covered_refs) >= 2:
            merged_ref_ids.update(covered_refs)
            merging_pred_ids.add(pred.feat_id)

    # --- Step 2: detect fragmented references ---
    fragmented_ref_ids: set[str] = set()
    fragmenting_pred_ids: set[str] = set()

    for ref in ref_genes:
        if ref.feat_id in merged_ref_ids:
            continue
        pred_overlaps = _find_overlaps(ref, pred_idx)
        if len(pred_overlaps) < 2:
            continue
        qualifying = [
            pred_iv.feat_id
            for pred_iv, ovl_bp in pred_overlaps
            if ovl_bp / (pred_iv.end - pred_iv.start + 1) >= FRAG_THRESHOLD
        ]
        if len(qualifying) >= 2:
            fragmented_ref_ids.add(ref.feat_id)
            fragmenting_pred_ids.update(qualifying)

    # --- Step 3: classify remaining refs and build match map ---
    gene_match_map: dict[str, str] = {}
    excluded_pred_ids = merging_pred_ids | fragmenting_pred_ids

    for ref in ref_genes:
        if ref.feat_id in merged_ref_ids:
            records.append(FeatureRecord.from_ref("gene", "merged", ref))
            continue

        if ref.feat_id in fragmented_ref_ids:
            records.append(FeatureRecord.from_ref("gene", "fragmented", ref))
            continue

        pred_overlaps = [
            (p, o) for p, o in _find_overlaps(ref, pred_idx)
            if p.feat_id not in excluded_pred_ids
        ]

        if not pred_overlaps:
            records.append(FeatureRecord.from_ref("gene", "missing", ref))
            continue

        best_pred, best_ovl = max(pred_overlaps, key=lambda x: x[1])

        if ref.start == best_pred.start and ref.end == best_pred.end:
            records.append(FeatureRecord.from_match(
                "gene", "exact", ref, best_pred, best_ovl,
            ))
        else:
            if ref.strand == "+":
                b5p = best_pred.start - ref.start
                b3p = best_pred.end - ref.end
            else:
                b5p = ref.end - best_pred.end
                b3p = ref.start - best_pred.start
            records.append(FeatureRecord.from_match(
                "gene", "inexact", ref, best_pred, best_ovl,
                boundary_5p=b5p, boundary_3p=b3p,
            ))

        gene_match_map[ref.feat_id] = best_pred.feat_id

    # --- Step 4: novel predictions ---
    matched_pred_ids = set(gene_match_map.values()) | excluded_pred_ids
    for pred in pred_genes:
        if pred.feat_id not in matched_pred_ids:
            if not _find_overlaps(pred, ref_idx):
                records.append(FeatureRecord.from_pred("gene", "novel", pred))

    return gene_match_map, records


# ---------------------------------------------------------------------------
# Subfeature-level classification (mRNA, CDS, intron)
# ---------------------------------------------------------------------------


def _classify_subfeatures(
    level_name: str,
    ref_feats_by_parent: dict[str, list[Interval]],
    pred_feats_by_parent: dict[str, list[Interval]],
    parent_match_map: dict[str, str],
) -> tuple[dict[str, str], list[FeatureRecord], int, int]:
    """Classify subfeatures within matched parent pairs.

    Returns (match_map, records, ref_total, pred_total).
    """
    level_key = level_name.lower()
    match_map: dict[str, str] = {}
    records: list[FeatureRecord] = []
    ref_total = 0
    pred_total = 0

    for ref_parent_id, pred_parent_id in parent_match_map.items():
        ref_feats = ref_feats_by_parent.get(ref_parent_id, [])
        pred_feats = pred_feats_by_parent.get(pred_parent_id, [])
        ref_total += len(ref_feats)
        pred_total += len(pred_feats)

        matches = _max_pairwise_matching(ref_feats, pred_feats)

        matched_ref_ids: set[str] = set()
        matched_pred_ids: set[str] = set()

        for ref_iv, pred_iv, ovl_bp in matches:
            matched_ref_ids.add(ref_iv.feat_id)
            matched_pred_ids.add(pred_iv.feat_id)
            match_map[ref_iv.feat_id] = pred_iv.feat_id
            records.append(FeatureRecord.from_match(
                level_key, "match", ref_iv, pred_iv, ovl_bp,
            ))

        for ref_iv in ref_feats:
            if ref_iv.feat_id not in matched_ref_ids:
                records.append(FeatureRecord.from_ref(
                    level_key, "missing", ref_iv,
                ))

        for pred_iv in pred_feats:
            if pred_iv.feat_id not in matched_pred_ids:
                records.append(FeatureRecord.from_pred(
                    level_key, "fp", pred_iv,
                ))

    return match_map, records, ref_total, pred_total


# ---------------------------------------------------------------------------
# Main comparison driver
# ---------------------------------------------------------------------------


def compare_annotations(ref_path: Path, pred_path: Path) -> ComparisonResult:
    """Compare predicted annotations against a reference GFF3."""
    ref_db = _load_gff3(ref_path)
    pred_db = _load_gff3(pred_path)

    # --- Parse features ---
    ref_genes = _parse_features(ref_db, "gene")
    pred_genes = _parse_features(pred_db, "gene")
    ref_mrnas = _parse_features(ref_db, "mRNA", parent_attr=True)
    pred_mrnas = _parse_features(pred_db, "mRNA", parent_attr=True)
    ref_cds = _parse_features(ref_db, "CDS", parent_attr=True)
    pred_cds = _parse_features(pred_db, "CDS", parent_attr=True)

    # Group by parent
    ref_mrnas_by_gene = _group_by_parent(ref_mrnas)
    pred_mrnas_by_gene = _group_by_parent(pred_mrnas)
    ref_cds_by_mrna = _group_by_parent(ref_cds)
    pred_cds_by_mrna = _group_by_parent(pred_cds)

    # Derive introns from CDS gaps
    ref_introns_by_mrna = {
        mid: introns
        for mid, cds in ref_cds_by_mrna.items()
        if (introns := _derive_introns(cds))
    }
    pred_introns_by_mrna = {
        mid: introns
        for mid, cds in pred_cds_by_mrna.items()
        if (introns := _derive_introns(cds))
    }

    # --- Gene-level classification ---
    gene_match_map, gene_records = _classify_genes(ref_genes, pred_genes)
    gene_stats = gene_stats_from_records(
        gene_records, len(ref_genes), len(pred_genes),
    )

    # --- mRNA-level classification (within matched gene pairs) ---
    mrna_match_map, mrna_records, mrna_ref_total, mrna_pred_total = (
        _classify_subfeatures(
            "mRNA", ref_mrnas_by_gene, pred_mrnas_by_gene, gene_match_map,
        )
    )
    mrna_stats = subfeature_stats_from_records(
        mrna_records, "mRNA", mrna_ref_total, mrna_pred_total,
    )

    # --- CDS-level classification (within matched mRNA pairs) ---
    _, cds_records, cds_ref_total, cds_pred_total = _classify_subfeatures(
        "CDS", ref_cds_by_mrna, pred_cds_by_mrna, mrna_match_map,
    )
    cds_stats = subfeature_stats_from_records(
        cds_records, "CDS", cds_ref_total, cds_pred_total,
    )

    # --- Intron-level classification (within matched mRNA pairs) ---
    _, intron_records, intron_ref_total, intron_pred_total = (
        _classify_subfeatures(
            "Intron", ref_introns_by_mrna, pred_introns_by_mrna,
            mrna_match_map,
        )
    )
    intron_stats = subfeature_stats_from_records(
        intron_records, "Intron", intron_ref_total, intron_pred_total,
    )

    all_records = gene_records + mrna_records + cds_records + intron_records

    return ComparisonResult(
        gene_stats=gene_stats,
        mrna_stats=mrna_stats,
        cds_stats=cds_stats,
        intron_stats=intron_stats,
        ref_path=str(ref_path),
        pred_path=str(pred_path),
        records=all_records,
    )
