"""Final consensus model building: EVM + optional PASA UTRs + prettification."""

from __future__ import annotations

from pathlib import Path

import gffutils

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


def build_consensus_models(config: PipelineConfig, *evidence: Path) -> Path:
    """Build final consensus models from all predictions using EVM."""
    sdir = step_dir(config.work_dir, "evm_consensus_models")
    log.info("Building consensus gene models...")

    run_evm(config, list(evidence))

    evm_genes = count_gff3_features(sdir / "consensus_models.gff3")
    log.info("EVM consensus: %d genes", evm_genes)

    if config.utrs_db:
        add_utrs_from_pasa(config, sdir, config.utrs_db)

    pasa_outputs = sorted(sdir.glob("*gene_structures_post_PASA_updates.*.gff3"))
    consensus_path = pasa_outputs[0] if pasa_outputs else sdir / "consensus_models.gff3"

    # Format the final GFF3
    consdb = create_gff_db(consensus_path)

    # Patch in transcript ORFs that don't overlap consensus models
    orf_path = config.work_dir / "orf_finder" / "transcript_orfs.gff3"
    if orf_path.exists():
        orf_db = create_gff_db(orf_path)
        missing = gffintersecter.find_nonoverlapping_genes(orf_db, consdb)
        if missing:
            log.info("Reintroduced %d transcript ORFs not overlapping EVM consensus", len(missing))
        all_features = list(consdb.all_features()) + missing
        consdb = gffutils.create_db(
            all_features, ":memory:",
            merge_strategy="merge", from_string=True,
        )

    consdb.dialect["order"].append("locus_tag")
    consdb = create_gff_db(
        gffparser.fix_CDS_phases(consdb), merge_strategy="merge",
    )
    prettified = gffparser.prettify_gff3(consdb, config.shortname)

    final_path = config.work_dir / "final.gff3"
    with open(final_path, "w") as outfile:
        for feature in prettified:
            outfile.write(f"{feature}\n")

    return final_path
