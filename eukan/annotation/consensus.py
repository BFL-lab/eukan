"""Final consensus model building: combinr consensus + prettification."""

from __future__ import annotations

from pathlib import Path

import gffutils

from eukan.annotation.combinr_consensus import run_combinr_consensus
from eukan.gff import concordance, create_gff_db
from eukan.gff.hierarchy import ORF_PATCH_SOURCE, prettify_gff3
from eukan.gff.intervals import IntervalIndex
from eukan.gff.io import count_gff3_features
from eukan.gff.rows import iter_gff3_rows
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def _patch_in_unmatched_orfs(consdb, orf_path: Path, protein_alignments: Path | None = None):
    """Add transcript ORFs that don't overlap any consensus gene, on either strand.

    With *protein_alignments* (``--orf-patch supported``), only those whose CDS
    span overlaps one of its CDS rows on the same strand.

    Every row added says so in column 2 (:data:`ORF_PATCH_SOURCE`, which
    :func:`prettify_gff3` keeps) and every mRNA added carries
    ``support=transcript_orf``, beside combinr's ``support=consensus``.

    Returns ``consdb``, updated in place when there are ORFs to add.
    """
    if not orf_path.exists():
        return consdb
    orf_db = create_gff_db(orf_path)
    support = None
    if protein_alignments is not None:
        support = IntervalIndex(
            gffutils.Feature(seqid=c[0], strand=c[6], start=int(c[3]), end=int(c[4]))
            for c in iter_gff3_rows(protein_alignments) if c[2] == "CDS"
        )
    missing = concordance.find_nonoverlapping_genes(
        orf_db, consdb, ignore_strand=True, support=support,
    )
    if not missing:
        return consdb
    for f in missing:
        f.source = ORF_PATCH_SOURCE
        if f.featuretype == "mRNA":
            f.attributes["support"] = ["transcript_orf"]
    log.info(
        "Reintroduced %d transcript ORFs not overlapping the consensus",
        sum(f.featuretype == "gene" for f in missing),
    )
    return consdb.update(missing, merge_strategy="merge")


def _write_prettified_gff3(consdb, shortname: str, out_path: Path) -> None:
    """Fix CDS phases, assign locus tags, and write the final GFF3."""
    consdb.dialect["order"].append("locus_tag")
    with open(out_path, "w") as outfile:
        for feature in prettify_gff3(consdb, shortname):
            outfile.write(f"{feature}\n")


def build_consensus_models(
    config: PipelineConfig,
    *evidence: Path,
    transcripts: Path | None = None,
) -> Path:
    """Build final consensus models from all predictions.

    Phases:
      1. Consensus building into consensus_models.gff3 via combinr consensus
         (which folds UTRs and isoforms in itself via --alt-splice).
      2. Under --orf-patch, patch in transcript ORFs that overlap no consensus
         gene (``supported``: only those on a protein alignment's CDS)
      3. Recompute CDS phases, assign locus tags, write final.gff3

    ``transcripts`` is the transcript-alignments evidence: the assembled
    ``nr_transcripts.gff3`` when transcripts were assembled, else ``None``.
    """
    sdir = step_dir(config.work_dir, "consensus_models")
    log.info("Building consensus gene models...")

    run_combinr_consensus(config, sdir, list(evidence), transcripts=transcripts)
    log.info(
        "combinr consensus: %d gene predictions",
        count_gff3_features(sdir / "consensus_models.gff3"),
    )
    consdb = create_gff_db(sdir / "consensus_models.gff3")

    orfs = config.work_dir / "orf_finder" / "transcript_orfs.gff3"
    if config.orf_patch is None:
        if orfs.exists():  # a protein-only run has none to reintroduce
            log.info("ORF patch off; pass --orf-patch supported|all to reintroduce transcript ORFs")
    elif config.orf_patch == "all":
        consdb = _patch_in_unmatched_orfs(consdb, orfs)
    else:
        # The protein alignment among the evidence, picked as combinr's staging picks it.
        prot = next(ev for ev in evidence if ev.name == "prot.gff3")
        consdb = _patch_in_unmatched_orfs(consdb, orfs, protein_alignments=prot)

    final_path = config.work_dir / Artifact.FINAL_GFF3.value
    _write_prettified_gff3(consdb, config.shortname, final_path)
    log.info(
        "Final consensus: %d gene predictions", count_gff3_features(final_path),
    )
    return final_path
