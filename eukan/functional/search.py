"""Homology search via pyhmmer (phmmer + hmmscan) and result annotation."""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path
from typing import TypedDict

import gffutils
import pyhmmer
from Bio import SeqIO

from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger

log = get_logger(__name__)


def _decode(value: str | bytes | None) -> str:
    """Decode bytes to str, pass through if already str (pyhmmer compat)."""
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else value


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


# ---------------------------------------------------------------------------
# Sequence loading
# ---------------------------------------------------------------------------


def _load_digital_sequences(fasta_path: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """Load protein sequences from a FASTA file into pyhmmer digital format."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(str(fasta_path), digital=True, alphabet=alphabet) as sf:
        seqs: list[pyhmmer.easel.DigitalSequence] = list(sf)  # type: ignore[arg-type]
    return seqs


def _load_hmm_db(hmm_path: Path) -> list[pyhmmer.plan7.HMM]:
    """Load HMM profiles from a pressed Pfam database."""
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms: list[pyhmmer.plan7.HMM] = list(hf)
    return hmms


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
        query_name = _decode(top_hits.query.name)
        if not top_hits:
            continue
        # Keep only the best hit (lowest e-value)
        best = top_hits[0]
        hit_name = _decode(best.name)
        hit_desc = _decode(best.description) if best.description else ""
        results[query_name] = {
            hit_name: {"description": hit_desc, "evalue": best.evalue}
        }

    log.info("phmmer: %d queries with hits", len(results))
    return results


def run_hmmscan_search(
    queries: list[pyhmmer.easel.DigitalSequence],
    hmms: list[pyhmmer.plan7.HMM],
    num_cpu: int,
    evalue_threshold: float,
) -> HitResults:
    """Run hmmscan (sequence vs HMM profiles) and return non-overlapping domain hits.

    For each query, keeps hits whose domains don't overlap on the query sequence.
    """
    results: HitResults = {}

    for top_hits in pyhmmer.hmmer.hmmscan(queries, hmms, cpus=num_cpu, E=evalue_threshold):
        query_name = _decode(top_hits.query.name)
        if not top_hits:
            continue

        query_hits: dict[str, HitInfo] = {}
        prev_end = -1

        for hit in top_hits:
            if not hit.domains:
                continue
            # Use the best domain for range checking
            best_domain = hit.domains[0]
            dom_start = best_domain.alignment.target_from
            dom_end = best_domain.alignment.target_to

            # Only keep non-overlapping domains
            if dom_start > prev_end:
                hit_name = _decode(hit.name)
                hit_desc = _decode(hit.description) if hit.description else ""
                query_hits[hit_name] = {
                    "description": hit_desc,
                    "evalue": hit.evalue,
                    "accession": _decode(hit.accession),
                }
                prev_end = dom_end

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
    homology_db: str = "uniprot",
) -> Path:
    """Annotate protein FASTA headers with functional information.

    ``homology_db`` selects the source of the homology hits cached in
    ``homology_res`` and only affects the wording of the marginal-hit
    fallback label. Returns the path to the annotated output file.
    """
    if homology_db == "kofam":
        homology_fmt = _format_kofam_results(homology_res)
    else:
        homology_fmt = _format_results(
            homology_res, marginal_label="hypothetical protein [{desc}]",
        )
    hmmscan_fmt = _format_results(hmmscan_res, marginal_label="{desc} [marginal domain hit]")

    output_path = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"

    with open(output_path, "w") as f:
        for rec in SeqIO.parse(str(proteins), "fasta"):
            parts = [homology_fmt.get(rec.id, "hypothetical protein"), f"length={len(rec.seq)}"]
            if rec.id in hmmscan_fmt:
                parts.append(hmmscan_fmt[rec.id])
            rec.description = " ;; ".join(parts)
            SeqIO.write(rec, f, "fasta")

    return output_path


def _format_results(results: HitResults, marginal_label: str = "{desc} [marginal hit]") -> dict[str, str]:
    """Format UniProt/Pfam-style search results with marginal hit annotations."""
    formatted = {}
    for query_id, hits in results.items():
        parts = []
        for hit_id, info in hits.items():
            desc = info["description"]
            ev = info["evalue"]
            if ev >= 1e-3:
                desc = marginal_label.format(desc=desc)
            parts.append(f"{hit_id}: {desc} ({ev:.2e})")
        formatted[query_id] = " ;; ".join(parts)
    return formatted


def _format_kofam_results(results: HitResults) -> dict[str, str]:
    """Format KOfam hits for FASTA headers.

    Above-threshold hits are emitted verbatim; below-threshold hits are
    tagged ``[marginal KO hit]`` so the user can tell at a glance whether
    a record's KO assignment passed the per-KO cutoff.
    """
    formatted = {}
    for query_id, hits in results.items():
        parts = []
        for k_number, info in hits.items():
            desc = info["description"]
            ev = info["evalue"]
            tag = "" if info.get("above_threshold") else " [marginal KO hit]"
            parts.append(f"{k_number}: {desc}{tag} ({ev:.2e})")
        formatted[query_id] = " ;; ".join(parts)
    return formatted


# ---------------------------------------------------------------------------
# GFF3 annotation
# ---------------------------------------------------------------------------


def annotate_gff3(
    gff3_path: Path,
    homology_res: HitResults,
    hmmscan_res: HitResults,
    output_dir: Path | None = None,
    homology_db: str = "uniprot",
) -> Path:
    """Annotate GFF3 features with functional information.

    Output filename is ``<stem>.mod<suffix>``. When *output_dir* is given,
    the file lands there; otherwise it lands next to *gff3_path*. The
    func-annot pipeline passes its work_dir so the result matches the
    ``Artifact.FINAL_FUNC_GFF3`` convention (``func-annot/final.mod.gff3``)
    that prep-submission auto-discovery and downstream tooling rely on.

    Returns path to the annotated output file.
    """
    hmmscan_strict = {
        qid: {hid: info for hid, info in hits.items() if info["evalue"] <= 1e-2}
        for qid, hits in hmmscan_res.items()
        if any(info["evalue"] <= 1e-2 for info in hits.values())
    }

    from eukan.gff import create_gff_db
    gff3 = create_gff_db(gff3_path)
    gff3.update(
        _add_func_info(gff3, homology_res, hmmscan_strict, homology_db),
        merge_strategy="replace",
    )

    out_dir = output_dir or gff3_path.parent
    output_path = out_dir / f"{gff3_path.stem}.mod{gff3_path.suffix}"
    featuredb2gff3_file(gff3, output_path)

    return output_path


def _apply_uniprot_hits(f: gffutils.Feature, hits: dict[str, HitInfo]) -> None:
    """Set product/inference from UniProt phmmer hits (E-value gated)."""
    f.attributes["product"] = [
        v["description"] if v["evalue"] <= 1e-2 else "hypothetical protein"
        for v in hits.values()
    ]
    f.attributes["inference"] = [
        f"similar to AA sequence:UniProtKB:{k}"
        for k, v in hits.items()
        if v["evalue"] <= 1e-2
    ]


def _apply_kofam_hits(f: gffutils.Feature, hits: dict[str, HitInfo]) -> None:
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
    homology_db: str = "uniprot",
) -> Iterator[gffutils.Feature]:
    """Yield mRNA/CDS features annotated with functional information."""
    for f in gff3.features_of_type(["mRNA", "CDS"]):
        f.attributes["locus_tag"] = list(f.attributes["ID"])
        feature_id = f.attributes["ID"][0]
        parent_id = f.attributes.get("Parent", [None])[0]

        lookup_id = feature_id if f.featuretype == "mRNA" else parent_id

        if lookup_id and lookup_id in homology_res:
            hits = homology_res[lookup_id]
            if homology_db == "kofam":
                _apply_kofam_hits(f, hits)
            else:
                _apply_uniprot_hits(f, hits)
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
