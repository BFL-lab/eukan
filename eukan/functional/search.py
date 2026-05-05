"""Homology search via pyhmmer (phmmer + hmmscan) and result annotation."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import TypedDict

import gffutils
import pyhmmer
from Bio import SeqIO

from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger

log = get_logger(__name__)


def _decode(value: str | bytes) -> str:
    """Decode bytes to str, pass through if already str (pyhmmer compat)."""
    return value.decode() if isinstance(value, bytes) else value


class HitInfo(TypedDict):
    description: str
    evalue: float


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
                }
                prev_end = dom_end

        if query_hits:
            results[query_name] = query_hits

    log.info("hmmscan: %d queries with domain hits", len(results))
    return results


def run_homology_search(
    proteins: Path,
    uniprot_db: Path,
    pfam_db: Path,
    num_cpu: int,
    evalue: str,
) -> tuple[HitResults, HitResults]:
    """Run phmmer and hmmscan using pyhmmer.

    Stages are run sequentially with the target database released
    between them: holding both UniProt-SwissProt and Pfam in memory
    while two searches ran concurrently was OOM-prone on container
    runtimes.  Each stage uses the full CPU budget.

    Returns:
        Tuple of (phmmer_results, hmmscan_results) dicts keyed by query ID.
    """
    import gc

    evalue_f = float(evalue)

    log.info("Loading query sequences from %s", proteins)
    queries = _load_digital_sequences(proteins)

    log.info("Loading UniProt database from %s", uniprot_db)
    targets = _load_digital_sequences(uniprot_db)
    log.info(
        "Running phmmer (%d queries vs %d targets, %d CPUs)...",
        len(queries), len(targets), num_cpu,
    )
    phmmer_res = run_phmmer_search(queries, targets, num_cpu, evalue_f)
    del targets
    gc.collect()

    log.info("Loading Pfam HMMs from %s", pfam_db)
    hmms = _load_hmm_db(pfam_db)
    log.info(
        "Running hmmscan (%d queries vs %d profiles, %d CPUs)...",
        len(queries), len(hmms), num_cpu,
    )
    hmmscan_res = run_hmmscan_search(queries, hmms, num_cpu, evalue_f)
    del hmms
    gc.collect()

    return phmmer_res, hmmscan_res


# ---------------------------------------------------------------------------
# FASTA annotation
# ---------------------------------------------------------------------------


def annotate_fasta(
    proteins: Path,
    phmmer_res: HitResults,
    hmmscan_res: HitResults,
) -> Path:
    """Annotate protein FASTA headers with functional information.

    Returns path to the annotated output file.
    """
    phmmer_fmt = _format_results(phmmer_res, marginal_label="hypothetical protein [{desc}]")
    hmmscan_fmt = _format_results(hmmscan_res, marginal_label="{desc} [marginal domain hit]")

    output_path = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"

    with open(output_path, "w") as f:
        for rec in SeqIO.parse(str(proteins), "fasta"):
            parts = [phmmer_fmt.get(rec.id, "hypothetical protein"), f"length={len(rec.seq)}"]
            if rec.id in hmmscan_fmt:
                parts.append(hmmscan_fmt[rec.id])
            rec.description = " ;; ".join(parts)
            SeqIO.write(rec, f, "fasta")

    return output_path


def _format_results(results: HitResults, marginal_label: str = "{desc} [marginal hit]") -> dict[str, str]:
    """Format search results with marginal hit annotations."""
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


# ---------------------------------------------------------------------------
# GFF3 annotation
# ---------------------------------------------------------------------------


def annotate_gff3(
    gff3_path: Path,
    phmmer_res: HitResults,
    hmmscan_res: HitResults,
) -> Path:
    """Annotate GFF3 features with functional information.

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
        _add_func_info(gff3, phmmer_res, hmmscan_strict),
        merge_strategy="replace",
    )

    stem = gff3_path.stem
    suffix = gff3_path.suffix
    output_path = gff3_path.parent / f"{stem}.mod{suffix}"
    featuredb2gff3_file(gff3, output_path)

    return output_path


def _add_func_info(
    gff3: gffutils.FeatureDB,
    phmmer_res: HitResults,
    hmmscan_res: HitResults,
) -> Iterator[gffutils.Feature]:
    """Yield mRNA/CDS features annotated with functional information."""
    for f in gff3.features_of_type(["mRNA", "CDS"]):
        f.attributes["locus_tag"] = list(f.attributes["ID"])
        feature_id = f.attributes["ID"][0]
        parent_id = f.attributes.get("Parent", [None])[0]

        lookup_id = feature_id if f.featuretype == "mRNA" else parent_id

        if lookup_id and lookup_id in phmmer_res:
            hits = phmmer_res[lookup_id]
            f.attributes["product"] = [
                v["description"] if v["evalue"] <= 1e-2 else "hypothetical protein"
                for v in hits.values()
            ]
            f.attributes["inference"] = [
                f"similar to AA sequence:UniProtKB:{k}"
                for k, v in hits.items()
                if v["evalue"] <= 1e-2
            ]
        else:
            f.attributes["product"] = ["hypothetical protein"]

        if lookup_id and lookup_id in hmmscan_res:
            pfam_inferences = [
                f"protein motif:PFAM:{k}" for k in hmmscan_res[lookup_id]
            ]
            if "inference" in f.attributes:
                f.attributes["inference"] = f.attributes["inference"] + pfam_inferences
            else:
                f.attributes["inference"] = pfam_inferences

        yield f
