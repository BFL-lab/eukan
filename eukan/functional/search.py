"""Homology search via pyhmmer (phmmer + hmmscan) and result annotation."""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TypedDict

import gffutils
import pyhmmer
from Bio import SeqIO

from eukan.gff import create_gff_db
from eukan.gff.io import count_gff3_features, featuredb2gff3_file
from eukan.infra.logging import get_logger

log = get_logger(__name__)


class HitInfo(TypedDict, total=False):
    description: str
    evalue: float
    accession: str  # optional — backfilled for legacy caches via _pfam_accession
    # --- KOfam-only fields (present when search came from run_kofam_search) ---
    score: float
    threshold: float       # 0.0 when ko_list had no curated threshold ("-")
    score_type: str        # "full" or "domain"
    above_threshold: bool  # score >= per-KO threshold
    ec_numbers: list[str]  # EC codes parsed out of the KO definition


# E-value at or below which an hmmscan/phmmer hit is treated as a confident
# ("good") functional assignment in the annotated GFF3; above it the hit is
# marginal (e.g. product falls back to "hypothetical protein").
_EVALUE_GOOD_HIT = 1e-2

# The annotated FASTA headers are stricter: from this E-value up a UniProt hit
# reads "hypothetical protein [<description>]" and a Pfam hit
# "<description> [marginal domain hit]".
_EVALUE_MARGINAL_IN_FASTA = 1e-3

_PFAM_VERSION_RE = re.compile(r"\.\d+$")


def _pfam_accession(info: HitInfo, hit_name: str) -> str:
    """Return the Pfam accession (e.g. ``PF00710``) for an hmmscan hit.

    Strips the trailing version suffix (``PF00710.15`` → ``PF00710``).
    Falls back to the hit name when the cached result predates accession
    capture or the HMM has no accession metadata.
    """
    acc = info.get("accession", "")
    if not acc:
        return hit_name
    return _PFAM_VERSION_RE.sub("", acc)


# query_id -> {hit_id -> HitInfo}
HitResults = dict[str, dict[str, HitInfo]]

# How the hits of one homology DB are written out: one hit as a FASTA header
# field, and a query's hits onto its mRNA/CDS feature's GFF3 attributes.
HitFormatter = Callable[[str, HitInfo], str]
HitApplier = Callable[[gffutils.Feature, dict[str, HitInfo]], None]


# ---------------------------------------------------------------------------
# Sequence loading
# ---------------------------------------------------------------------------


def _load_digital_sequences(fasta_path: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """Load protein sequences from a FASTA file into pyhmmer digital format."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(str(fasta_path), digital=True, alphabet=alphabet) as sf:
        seqs: list[pyhmmer.easel.DigitalSequence] = list(sf)  # type: ignore[arg-type]
    return seqs


@contextmanager
def _open_hmm_db(
    hmm_path: Path,
) -> Iterator[Iterable[pyhmmer.plan7.HMM | pyhmmer.plan7.OptimizedProfile]]:
    """Open an HMM database as the profiles to hand to hmmscan.

    A pressed database (``eukan db-fetch`` presses the ones it fetches) gives
    its optimized profiles, read straight into one block; an unpressed one
    gives the open file, which hmmscan configures into a block itself.
    """
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        if hf.is_pressed():
            yield pyhmmer.plan7.OptimizedProfileBlock(
                pyhmmer.easel.Alphabet.amino(), hf.optimized_profiles(),
            )
        else:
            yield hf


# ---------------------------------------------------------------------------
# Homology search via pyhmmer
# ---------------------------------------------------------------------------


def run_phmmer_search(
    queries: list[pyhmmer.easel.DigitalSequence],
    targets: list[pyhmmer.easel.DigitalSequence],
    num_cpu: int,
    evalue_threshold: float,
) -> HitResults:
    """Run phmmer (sequence vs sequence) and return best hit per query."""
    results: HitResults = {}

    for top_hits in pyhmmer.hmmer.phmmer(queries, targets, cpus=num_cpu, E=evalue_threshold):
        query_name = top_hits.query.name
        # Keep only the best hit (lowest e-value). E= gates hits in-scan
        # against the targets seen so far; only the *reported* ones pass it
        # against the whole database.
        best = next(top_hits.reported, None)
        if best is None:
            continue
        hit_name = best.name
        hit_desc = best.description or ""
        results[query_name] = {
            hit_name: {"description": hit_desc, "evalue": best.evalue}
        }

    log.info("phmmer: %d queries with hits", len(results))
    return results


def run_hmmscan_search(
    queries: list[pyhmmer.easel.DigitalSequence],
    profiles: Iterable[pyhmmer.plan7.HMM | pyhmmer.plan7.OptimizedProfile],
    num_cpu: int,
    evalue_threshold: float,
) -> HitResults:
    """Run hmmscan (sequence vs HMM profiles) and return non-overlapping domain hits.

    For each query, hits are taken best E-value first; each is placed by its
    best-scoring domain and kept unless that domain overlaps, on the query
    sequence, the domain of a hit already kept.
    """
    results: HitResults = {}

    for top_hits in pyhmmer.hmmer.hmmscan(queries, profiles, cpus=num_cpu, E=evalue_threshold):
        query_name = top_hits.query.name
        if not top_hits:
            continue

        query_hits: dict[str, HitInfo] = {}
        kept_spans: list[tuple[int, int]] = []

        # Only *reported* hits pass E= against the whole database (see
        # run_phmmer_search).
        for hit in top_hits.reported:
            if not hit.domains:
                continue
            aln = hit.best_domain.alignment
            start, end = aln.target_from, aln.target_to
            if any(start <= k_end and k_start <= end for k_start, k_end in kept_spans):
                continue
            kept_spans.append((start, end))
            hit_name = hit.name
            hit_desc = hit.description or ""
            query_hits[hit_name] = {
                "description": hit_desc,
                "evalue": hit.evalue,
                "accession": hit.accession or "",
            }

        if query_hits:
            results[query_name] = query_hits

    log.info("hmmscan: %d queries with domain hits", len(results))
    return results


# ---------------------------------------------------------------------------
# FASTA annotation
# ---------------------------------------------------------------------------


def annotate_fasta(
    proteins: Path,
    homology_res: HitResults,
    hmmscan_res: HitResults,
    format_hit: HitFormatter,
) -> Path:
    """Annotate protein FASTA headers with functional information.

    ``format_hit`` writes one hit of ``homology_res``: :func:`format_uniprot_hit`
    or :func:`format_kofam_hit`, after the homology DB searched. Returns the
    path to the annotated output file.
    """
    homology_fmt = _format_hits(homology_res, format_hit)
    hmmscan_fmt = _format_hits(hmmscan_res, _format_pfam_hit)

    output_path = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"

    count = 0
    with open(output_path, "w") as f:
        for rec in SeqIO.parse(str(proteins), "fasta"):
            parts = [homology_fmt.get(rec.id, "hypothetical protein"), f"length={len(rec.seq)}"]
            if rec.id in hmmscan_fmt:
                parts.append(hmmscan_fmt[rec.id])
            rec.description = " ;; ".join(parts)
            SeqIO.write(rec, f, "fasta")
            count += 1

    log.info("Functional FASTA: annotated %d sequences -> %s", count, output_path.name)
    return output_path


def _format_hits(results: HitResults, format_hit: HitFormatter) -> dict[str, str]:
    """Each query's hits as one FASTA header field, one ``format_hit`` per hit."""
    return {
        query_id: " ;; ".join(format_hit(hit_id, info) for hit_id, info in hits.items())
        for query_id, hits in results.items()
    }


def _format_evalue_hit(hit_id: str, info: HitInfo, marginal_label: str) -> str:
    """A UniProt or Pfam hit; a marginal one's description becomes ``marginal_label``."""
    desc, ev = info["description"], info["evalue"]
    if ev >= _EVALUE_MARGINAL_IN_FASTA:
        desc = marginal_label.format(desc=desc)
    return f"{hit_id}: {desc} ({ev:.2e})"


def format_uniprot_hit(hit_id: str, info: HitInfo) -> str:
    """A UniProt hit for a FASTA header; a marginal one reads as a hypothetical protein."""
    return _format_evalue_hit(hit_id, info, "hypothetical protein [{desc}]")


def _format_pfam_hit(hit_id: str, info: HitInfo) -> str:
    return _format_evalue_hit(hit_id, info, "{desc} [marginal domain hit]")


def format_kofam_hit(k_number: str, info: HitInfo) -> str:
    """A KOfam hit for a FASTA header.

    Above-threshold hits are emitted verbatim; below-threshold hits are
    tagged ``[marginal KO hit]`` so the user can tell at a glance whether
    a record's KO assignment passed the per-KO cutoff.
    """
    tag = "" if info.get("above_threshold") else " [marginal KO hit]"
    return f"{k_number}: {info['description']}{tag} ({info['evalue']:.2e})"


# ---------------------------------------------------------------------------
# GFF3 annotation
# ---------------------------------------------------------------------------


def annotate_gff3(
    gff3_path: Path,
    homology_res: HitResults,
    hmmscan_res: HitResults,
    apply_hits: HitApplier,
    output_dir: Path | None = None,
) -> Path:
    """Annotate GFF3 features with functional information.

    ``apply_hits`` sets a feature's attributes from its ``homology_res`` hits:
    :func:`apply_uniprot_hits` or :func:`apply_kofam_hits`, after the homology
    DB searched.

    Output filename is ``<stem>.mod<suffix>``. When *output_dir* is given,
    the file lands there; otherwise it lands next to *gff3_path*. The
    func-annot pipeline passes its work_dir so the result matches the
    ``Artifact.FINAL_FUNC_GFF3`` convention (``func-annot/final.mod.gff3``)
    that prep-submission auto-discovery and downstream tooling rely on.

    Returns path to the annotated output file.
    """
    hmmscan_strict = {
        qid: {hid: info for hid, info in hits.items() if info["evalue"] <= _EVALUE_GOOD_HIT}
        for qid, hits in hmmscan_res.items()
        if any(info["evalue"] <= _EVALUE_GOOD_HIT for info in hits.values())
    }

    gff3 = create_gff_db(gff3_path)
    gff3.update(
        _add_func_info(gff3, homology_res, hmmscan_strict, apply_hits),
        merge_strategy="replace",
    )

    out_dir = output_dir or gff3_path.parent
    output_path = out_dir / f"{gff3_path.stem}.mod{gff3_path.suffix}"
    featuredb2gff3_file(gff3, output_path)

    log.info(
        "Functional GFF3: %d genes -> %s",
        count_gff3_features(output_path), output_path.name,
    )
    return output_path


def apply_uniprot_hits(f: gffutils.Feature, hits: dict[str, HitInfo]) -> None:
    """Set product/inference from UniProt phmmer hits (E-value gated)."""
    f.attributes["product"] = [
        v["description"] if v["evalue"] <= _EVALUE_GOOD_HIT else "hypothetical protein"
        for v in hits.values()
    ]
    f.attributes["inference"] = [
        f"similar to AA sequence:UniProtKB:{k}"
        for k, v in hits.items()
        if v["evalue"] <= _EVALUE_GOOD_HIT
    ]


def apply_kofam_hits(f: gffutils.Feature, hits: dict[str, HitInfo]) -> None:
    """Set product/ec_number/Dbxref/inference from KOfam hits (per-KO threshold gated).

    Only above-threshold KOs contribute to ``product``/``Dbxref``/
    ``inference``. The top-scoring above-threshold hit drives ``product``;
    all above-threshold hits emit ``Dbxref`` and ``inference``. EC numbers
    de-duplicate across hits before emission so multi-functional KOs that
    share an EC don't produce repeats.
    """
    above = [(k, info) for k, info in hits.items() if info.get("above_threshold")]
    if not above:
        f.attributes["product"] = ["hypothetical protein"]
        return

    above.sort(key=lambda kv: kv[1].get("score", 0.0), reverse=True)
    _, top_info = above[0]
    f.attributes["product"] = [top_info["description"] or "hypothetical protein"]

    ec_numbers: list[str] = []
    seen_ec: set[str] = set()
    for _k, info in above:
        for ec in info.get("ec_numbers", []) or []:
            if ec not in seen_ec:
                seen_ec.add(ec)
                ec_numbers.append(ec)
    if ec_numbers:
        f.attributes["ec_number"] = ec_numbers

    f.attributes["Dbxref"] = [f"KEGG:{k}" for k, _ in above]
    f.attributes["inference"] = [f"protein motif:KOFAM:{k}" for k, _ in above]


def _add_func_info(
    gff3: gffutils.FeatureDB,
    homology_res: HitResults,
    hmmscan_res: HitResults,
    apply_hits: HitApplier,
) -> Iterator[gffutils.Feature]:
    """Yield mRNA/CDS features annotated with functional information."""
    for f in gff3.features_of_type(["mRNA", "CDS"]):
        f.attributes["locus_tag"] = list(f.attributes["ID"])
        feature_id = f.attributes["ID"][0]
        parent_id = f.attributes.get("Parent", [None])[0]

        lookup_id = feature_id if f.featuretype == "mRNA" else parent_id

        if lookup_id and lookup_id in homology_res:
            apply_hits(f, homology_res[lookup_id])
        else:
            f.attributes["product"] = ["hypothetical protein"]

        if lookup_id and lookup_id in hmmscan_res:
            pfam_inferences = [
                f"protein motif:PFAM:{_pfam_accession(info, hit_name)}"
                for hit_name, info in hmmscan_res[lookup_id].items()
            ]
            if "inference" in f.attributes:
                f.attributes["inference"] = f.attributes["inference"] + pfam_inferences
            else:
                f.attributes["inference"] = pfam_inferences

        yield f
