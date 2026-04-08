"""GeneMark-ES/ET gene prediction."""

from __future__ import annotations

from pathlib import Path

import gffutils
import pandas as pd

from eukan.annotation.validation import validate_gff
from eukan.gff import GFF3_DIALECT, parser as gffparser
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_complete, step_dir
from eukan.infra.logging import get_logger
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def run_genemark(config: PipelineConfig, hints: Path | None = None) -> Path:
    """Run GeneMark-ES/ET gene prediction."""
    output = "genemark.gff3"
    existing = step_complete(config.work_dir, "genemark", output)
    if existing:
        return existing

    sdir = step_dir(config.work_dir, "genemark")
    log.info("Running GeneMark gene prediction...")

    if config.rnaseq_hints is not None:
        validate_gff(config.rnaseq_hints)

    hgd_flag = ["--fungus"] if config.is_fungus else []

    if hints is None:
        training_type = ["--ES"]
    else:
        # Extract intron hints for guided training
        hints_df = pd.read_csv(hints, sep="\t", header=None, low_memory=False)
        introns = hints_df[hints_df[2] == "intron"]
        introns.to_csv(sdir / "introns.gff", sep="\t", header=None, index=False)
        training_type = (
            ["--ET=introns.gff", "--et_score=3"]
            if len(introns) >= 150
            else ["--ES"]
        )

    if not (sdir / "genemark.gtf").exists():
        run_cmd(
            ["gmes_petap.pl", "--soft", "1000"]
            + training_type
            + [f"--cores={config.num_cpu}", f"--sequence={config.genome}"]
            + hgd_flag,
            cwd=sdir,
        )

    # Convert GeneMark GTF to GFF3
    gmgtf = gffutils.create_db(
        str(sdir / "genemark.gtf"), ":memory:", verbose=False,
        disable_infer_genes=True, disable_infer_transcripts=True,
        merge_strategy="create_unique",
        id_spec={"gene": "gene_id", "mRNA": "transcript_id", "transcript": "transcript_id",
                 "CDS": ["gene_id", "transcript_id"], "exon": ["gene_id", "transcript_id"],
                 "start_codon": ["gene_id", "transcript_id"], "stop_codon": ["gene_id", "transcript_id"]},
    )
    gmgff3 = gffutils.create_db(
        gmgtf, ":memory:", dialect=GFF3_DIALECT,
        transform=gffparser.gtf2gff3, verbose=False,
    )
    gmgff3.update(
        gffparser.add_missing_feats_to_gff3(gmgff3),
        merge_strategy="create_unique",
    )
    gmgff3 = gffutils.create_db(
        gmgff3, ":memory:", dialect=GFF3_DIALECT,
        transform=gffparser.fix_contig_names, verbose=False,
    )
    featuredb2gff3_file(gmgff3, sdir / output)
    return sdir / output
