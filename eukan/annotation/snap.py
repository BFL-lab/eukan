"""SNAP and CodingQuarry gene prediction."""

from __future__ import annotations

from pathlib import Path

import gffutils

from eukan.annotation.augustus import build_training_set
from eukan.gff import GFF3_DIALECT, create_gff_db, transform_db
from eukan.gff import parser as gffparser
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.infra.utils import symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def run_snap(config: PipelineConfig, *evidence: Path) -> Path:
    """Train and run SNAP gene prediction."""
    output = "snap.gff3"
    sdir = step_dir(config.work_dir, "snap")
    log.info("Running SNAP training and prediction...")
    symlink(config.genome, sdir / "genome.dna")

    build_training_set(config, evidence, sdir)
    gffparser.gff2zff(sdir / "training_set.gff3", str(config.genome), output_dir=sdir)

    run_cmd(["fathom", "-categorize", "1000", "genome.ann", "genome.dna"], cwd=sdir)
    run_cmd(["fathom", "-export", "1000", "-plus", "uni.ann", "uni.dna"], cwd=sdir)
    run_cmd(["forge", "export.ann", "export.dna"], cwd=sdir)
    run_cmd(["hmm-assembler.pl", "snap", "."], cwd=sdir, out_file="snap.hmm")
    run_cmd(
        ["snap", "-gff", "snap.hmm", str(config.genome)],
        cwd=sdir, out_file="snap.gff",
    )

    snap = create_gff_db(
        sdir / "snap.gff", transform=gffparser.Snap.make_featuretype_transform(),
    )
    dialect = gffutils.DataIterator(str(sdir / "snap.gff")).dialect
    dialect["fmt"] = "gtf"
    snap = gffutils.create_db(
        snap, ":memory:", dialect=dialect,
        gtf_transcript_key="Parent", gtf_gene_key="Parent",
    )
    snap.dialect = GFF3_DIALECT
    snap = transform_db(snap, gffparser.gff3_it)
    snap.update(
        gffparser.add_missing_feats_to_gff3(snap),
        merge_strategy="create_unique",
    )
    snap = transform_db(snap, gffparser.Snap.homogenize_source)
    featuredb2gff3_file(snap, sdir / output)
    return sdir / output


def run_codingquarry(config: PipelineConfig, evidence: Path) -> Path:
    """Run CodingQuarry gene prediction (fungi/protists)."""
    output = "codingquarry.gff3"
    sdir = step_dir(config.work_dir, "codingquarry")
    log.info("Running CodingQuarry gene prediction...")

    flag = "-t" if config.transcripts_gff else "-a"
    run_cmd(
        [
            "CodingQuarry",
            "-p", str(config.num_cpu),
            "-f", str(config.genome),
            flag, str(evidence),
        ],
        cwd=sdir,
    )

    cq = create_gff_db(sdir / "out" / "PredictedPass.gff3", transform=gffparser.CodingQuarry.add_mRNA)
    cq.update(
        gffparser.add_missing_feats_to_gff3(cq),
        merge_strategy="create_unique",
    )
    cq = transform_db(cq, gffparser.CodingQuarry.fix_dup_ids)
    featuredb2gff3_file(cq, sdir / output)
    return sdir / output
