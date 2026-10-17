"""GeneMark-ES/ET gene prediction."""

from __future__ import annotations

from pathlib import Path

import gffutils

from eukan.gff import create_gff_db, transform_db
from eukan.gff import transforms as gff_transforms
from eukan.gff.normalize import normalize_to_gff3
from eukan.gff.rows import iter_gff3_rows
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig
from eukan.validation import validate_gff

log = get_logger(__name__)


def _genemark_homogenize_source(f: gffutils.Feature) -> gffutils.Feature:
    """Set ``source`` to ``genemark`` for all features.

    GeneMark stamps column 2 as ``GeneMark.hmm`` (or ``GeneMark.hmm3`` in
    newer releases); genemark.gff3 says ``genemark`` whatever the release, as
    the AUGUSTUS, SNAP and CodingQuarry outputs carry their own names. Nothing
    downstream reads the column: the GeneMark models choose between spaln and
    gth, supply spaln's gene lengths (and its introns when there are no intron
    hints), and train AUGUSTUS, but never reach the consensus engine.
    """
    f.source = "genemark"
    return f


def _read_genemark_gtf(path: Path) -> gffutils.FeatureDB:
    """Parse GeneMark's GTF output into a normalised GFF3 FeatureDB.

    GeneMark emits GTF with ``gene_id`` / ``transcript_id`` keys on every
    feature (including start/stop codons). The id_spec walks each feature
    type to its identifying attributes; the second pass converts GTF to
    GFF3 conventions via ``gff_transforms.gtf2gff3``.
    """
    gmgtf = create_gff_db(
        path,
        disable_infer_genes=True, disable_infer_transcripts=True,
        id_spec={
            "gene": "gene_id",
            "mRNA": "transcript_id",
            "transcript": "transcript_id",
            "CDS": ["gene_id", "transcript_id"],
            "exon": ["gene_id", "transcript_id"],
            "start_codon": ["gene_id", "transcript_id"],
            "stop_codon": ["gene_id", "transcript_id"],
        },
    )
    return transform_db(gmgtf, gff_transforms.gtf2gff3)


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
        with open(sdir / "introns.gff", "w") as fout:
            for cols in iter_gff3_rows(hints):
                if cols[2] == "intron":
                    fout.write("\t".join(cols) + "\n")
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

    # Convert GeneMark GTF to normalized GFF3
    gmdb = _read_genemark_gtf(sdir / "genemark.gtf")
    return normalize_to_gff3(
        gmdb, sdir / output,
        post_transform=_genemark_homogenize_source,
        fix_contig_names=True,
    )
