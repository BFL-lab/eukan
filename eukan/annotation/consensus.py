"""Final consensus model building: EVM + optional PASA UTRs + prettification."""

from __future__ import annotations

from pathlib import Path

from eukan.annotation.evm import run_evm
from eukan.assembly.pasa import write_pasa_configs
from eukan.gff import create_gff_db
from eukan.gff import intersecter as gffintersecter
from eukan.gff import parser as gffparser
from eukan.infra.logging import count_gff3_features, get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def add_utrs_from_pasa(config: PipelineConfig, sdir: Path, pasa_db: Path) -> None:
    """Add UTRs and model alternative splicing from a PASA database."""
    write_pasa_configs(sdir, pasa_db)

    run_cmd(
        [
            "Load_Current_Gene_Annotations.dbi",
            "-c", "alignAssembly.config",
            "-g", str(config.genome),
            "-P", "consensus_models.gff3",
        ],
        cwd=sdir,
    )

    gc_args = config.genetic_code_obj.pasa_flag

    run_cmd(
        [
            "Launch_PASA_pipeline.pl",
            "-c", "annotCompare.config",
            "-A", "-g", str(config.genome),
            "--CPU", str(config.num_cpu),
            "-t", str(config.transcripts_fasta),
            *gc_args,
        ],
        cwd=sdir,
    )


def _resolve_consensus_path(sdir: Path) -> Path:
    """Pick the latest PASA-updated GFF3 if present, else raw EVM consensus."""
    pasa_outputs = sorted(sdir.glob("*gene_structures_post_PASA_updates.*.gff3"))
    return pasa_outputs[0] if pasa_outputs else sdir / "consensus_models.gff3"


def _patch_in_unmatched_orfs(consdb, orf_path: Path):
    """Add transcript ORFs that don't overlap any consensus gene.

    Returns either the original ``consdb`` (no changes) or a new FeatureDB
    with the missing ORFs merged in.
    """
    if not orf_path.exists():
        return consdb
    orf_db = create_gff_db(orf_path)
    missing = gffintersecter.find_nonoverlapping_genes(orf_db, consdb)
    if not missing:
        return consdb
    log.info("Reintroduced %d transcript ORFs not overlapping EVM consensus", len(missing))
    all_features = [*consdb.all_features(), *missing]
    return create_gff_db(all_features, merge_strategy="merge")


def _write_prettified_gff3(consdb, shortname: str, out_path: Path) -> None:
    """Fix CDS phases, assign locus tags, and write the final GFF3."""
    consdb.dialect["order"].append("locus_tag")
    consdb = create_gff_db(gffparser.fix_CDS_phases(consdb), merge_strategy="merge")
    with open(out_path, "w") as outfile:
        for feature in gffparser.prettify_gff3(consdb, shortname):
            outfile.write(f"{feature}\n")


def build_consensus_models(config: PipelineConfig, *evidence: Path) -> Path:
    """Build final consensus models from all predictions.

    Phases:
      1. EVM merges all evidence into consensus_models.gff3
      2. (optional) PASA UTR/altsplice update if config.utrs_db is set
      3. Patch in transcript ORFs that don't overlap any consensus gene
      4. Recompute CDS phases, assign locus tags, write final.gff3
    """
    sdir = step_dir(config.work_dir, "evm_consensus_models")
    log.info("Building consensus gene models...")

    run_evm(config, list(evidence))
    log.info(
        "EVM consensus: %d genes",
        count_gff3_features(sdir / "consensus_models.gff3"),
    )

    if config.utrs_db:
        add_utrs_from_pasa(config, sdir, config.utrs_db)

    consdb = create_gff_db(_resolve_consensus_path(sdir))
    consdb = _patch_in_unmatched_orfs(
        consdb, config.work_dir / "orf_finder" / "transcript_orfs.gff3",
    )

    final_path = config.work_dir / "final.gff3"
    _write_prettified_gff3(consdb, config.shortname, final_path)
    return final_path
