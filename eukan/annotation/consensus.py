"""Final consensus model building: combinr consensus + prettification."""

from __future__ import annotations

from pathlib import Path

from eukan.annotation.combinr_consensus import run_combinr_consensus
from eukan.gff import concordance, create_gff_db
from eukan.gff.hierarchy import fix_CDS_phases, prettify_gff3
from eukan.gff.io import count_gff3_features
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def _patch_in_unmatched_orfs(consdb, orf_path: Path):
    """Add transcript ORFs that don't overlap any consensus gene.

    Returns either the original ``consdb`` (no changes) or a new FeatureDB
    with the missing ORFs merged in.
    """
    if not orf_path.exists():
        return consdb
    orf_db = create_gff_db(orf_path)
    missing = concordance.find_nonoverlapping_genes(orf_db, consdb)
    if not missing:
        return consdb
    log.info("Reintroduced %d transcript ORFs not overlapping the consensus", len(missing))
    all_features = [*consdb.all_features(), *missing]
    return create_gff_db(all_features, merge_strategy="merge")


def _write_prettified_gff3(consdb, shortname: str, out_path: Path) -> None:
    """Fix CDS phases, assign locus tags, and write the final GFF3."""
    consdb.dialect["order"].append("locus_tag")
    consdb = create_gff_db(fix_CDS_phases(consdb), merge_strategy="merge")
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
      2. Patch in transcript ORFs that don't overlap any consensus gene
      3. Recompute CDS phases, assign locus tags, write final.gff3

    ``transcripts`` is the transcript-alignments evidence. It's the assembled
    ``nr_transcripts.gff3`` when transcripts were assembled, the GeneMark
    predictions standing in for the consensus engine in the protein-only
    non-fungus/protist branch, or ``None``. (combinr ignores the GeneMark
    stand-in — it gates transcript evidence on ``config.has_transcripts``.)
    """
    sdir = step_dir(config.work_dir, "consensus_models")
    log.info("Building consensus gene models...")

    run_combinr_consensus(config, sdir, list(evidence), transcripts=transcripts)
    log.info(
        "combinr consensus: %d gene predictions",
        count_gff3_features(sdir / "consensus_models.gff3"),
    )
    consdb = create_gff_db(sdir / "consensus_models.gff3")

    consdb = _patch_in_unmatched_orfs(
        consdb, config.work_dir / "orf_finder" / "transcript_orfs.gff3",
    )

    final_path = config.work_dir / Artifact.FINAL_GFF3.value
    _write_prettified_gff3(consdb, config.shortname, final_path)
    log.info(
        "Final consensus: %d gene predictions", count_gff3_features(final_path),
    )
    return final_path
