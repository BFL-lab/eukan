"""GeneMark-ES/ET gene prediction."""

from __future__ import annotations

from pathlib import Path

import gffutils

from eukan.gff import parser as gffparser
from eukan.gff import transform_db
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger, validate_gff
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def run_genemark(config: PipelineConfig, hints: Path | None = None) -> Path:
    """Run GeneMark-ES/ET gene prediction."""
    output = "genemark.gff3"
    sdir = step_dir(config.work_dir, "genemark")
    log.info("Running GeneMark gene prediction...")

    if config.rnaseq_hints is not None:
        validate_gff(config.rnaseq_hints)

    hgd_flag = ["--fungus"] if config.is_fungus else []

    # Determine training mode: ES (self-training) or ET (with RNA-seq intron hints)
    has_intron_hints = False
    if hints is not None:
        intron_count = 0
        with open(hints) as fin, open(sdir / "introns.gff", "w") as fout:
            for line in fin:
                cols = line.split("\t")
                if len(cols) >= 3 and cols[2] == "intron":
                    fout.write(line)
                    intron_count += 1
        has_intron_hints = intron_count >= 150

    training_type = (
        ["--ET=introns.gff", "--et_score=3"] if has_intron_hints else ["--ES"]
    )
    gcode_flag = config.genetic_code_obj.genemark_flag

    if not (sdir / "genemark.gtf").exists():
        run_cmd(
            [
                "gmes_petap.pl", "--soft", "1000",
                *training_type,
                f"--cores={config.num_cpu}", f"--sequence={config.genome}",
                *hgd_flag,
                *gcode_flag,
            ],
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
    gmgff3 = transform_db(gmgtf, gffparser.gtf2gff3)
    gmgff3.update(
        gffparser.add_missing_feats_to_gff3(gmgff3),
        merge_strategy="create_unique",
    )
    gmgff3 = transform_db(gmgff3, gffparser.fix_contig_names)
    featuredb2gff3_file(gmgff3, sdir / output)
    return sdir / output
