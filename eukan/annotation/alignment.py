"""Spliced protein alignment via spaln (intron-rich) or GenomeThreader (intron-poor)."""

from __future__ import annotations

import re
import statistics
from pathlib import Path

import gffutils
import pandas as pd

from eukan.annotation.validation import validate_fasta
from eukan.gff import GFF3_DIALECT, parser as gffparser
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_complete, step_dir
from eukan.infra.utils import symlink
from eukan.infra.logging import get_logger
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def align_proteins(
    config: PipelineConfig,
    gff3: Path,
    proteins: list[Path],
    intron_hints: Path | None = None,
) -> Path:
    """Spliced alignment of protein sequences to the genome."""
    output = "prot.gff3"
    existing = step_complete(config.work_dir, "prot_align", output)
    if existing:
        return existing

    sdir = step_dir(config.work_dir, "prot_align")
    log.info("Running protein alignment...")

    # Concatenate all protein files
    with open(sdir / "prots.faa", "w") as outfile:
        for prot_file in proteins:
            outfile.write(prot_file.read_text())

    validate_fasta(sdir / "prots.faa")
    symlink(config.genome, sdir / "genome")
    symlink(config.genome, sdir / "genome.gf")

    models = gffutils.create_db(str(gff3), ":memory:", merge_strategy="create_unique")

    # Decide between spaln (intron-rich) and gth (intron-poor)
    gene_count = 0
    intron_count = 0
    for gene in models.features_of_type("gene"):
        gene_count += 1
        intron_count += max(0, len(list(models.children(gene, featuretype="exon"))) - 1)

    if gene_count > 0 and intron_count / gene_count > 0.25:
        _run_spaln(config, sdir, models, intron_hints or gff3)
    else:
        _run_gth(sdir)

    # Create protein hints for AUGUSTUS
    prot_db = gffutils.create_db(
        str(sdir / output), ":memory:", transform=gffparser.prot2augustus_hints,
    )
    with open(sdir / "hints_protein.gff", "w") as fout:
        for feat in prot_db.all_features():
            fout.write(f"{feat}\n")

    return sdir / output


def _run_spaln(
    config: PipelineConfig,
    sdir: Path,
    models: gffutils.FeatureDB,
    intron_hints_path: Path,
) -> None:
    """Run spaln protein alignment (for intron-rich genomes)."""
    lengths = [abs(f.end - f.start) for f in models.features_of_type("gene")]
    mean_len = int(sum(lengths) / len(lengths))
    sd_len = int(statistics.stdev(lengths)) if len(lengths) >= 2 else 0
    maxlen = mean_len + 2 * sd_len

    run_cmd(["makdbs", "-KD", "genome"], cwd=sdir)
    run_cmd(["makblk.pl", "-Wgenome.bkp", "genome.gf"], cwd=sdir)

    # Compute intron length distribution
    introns = pd.read_csv(intron_hints_path, sep="\t", header=None)
    intron_lens = (introns[4] - introns[3]).value_counts().reset_index().sort_values("index")
    intron_lens.to_csv(sdir / "introns.ild", sep=" ", header=None, index=False)

    run_cmd(["fitild", "-d", "IldModel.txt", "introns.ild"], cwd=sdir, out_file="fitild.out")

    # Parse intron length distribution parameters
    int_dist_parts: list[str] = []
    with open(sdir / "fitild.out") as fp:
        for line in fp:
            if re.search(r"^introns\s", line):
                parts = line.replace("\t", " ").replace("\n", " ").split()[3:-3]
                int_dist_parts = [f"{float(x):f}" for x in parts]
    if int_dist_parts:
        del int_dist_parts[1]
    int_dist = " ".join(int_dist_parts)

    run_cmd(
        ["spaln", "-Wgenome.bkp", f"-XG{maxlen}", "-KP", "genome.gf"],
        cwd=sdir,
    )

    # Remove stale spaln output files to avoid interactive overwrite prompts
    for ext in (".grd", ".erd", ".qrd"):
        (sdir / f"prots{ext}").unlink(missing_ok=True)

    run_cmd(
        [
            "spaln", "-O12", "-Q7", "-yX", "-LS",
            f"-C{config.genetic_code}", f"-t{config.num_cpu}",
            f"-yI{int_dist}", "-dgenome", "prots.faa",
        ],
        cwd=sdir,
    )

    run_cmd(
        ["sortgrcd", "-F0", "-C60", "-O0", "prots.grd"],
        cwd=sdir, out_file="temp.gff3",
    )

    spaln_out = gffutils.create_db(
        str(sdir / "temp.gff3"), ":memory:",
        transform=gffparser.fix_spaln_cds_featuretype,
        merge_strategy="create_unique",
    )
    spaln_out = gffutils.create_db(
        spaln_out, ":memory:", dialect=GFF3_DIALECT,
        transform=gffparser.fix_spaln_ids,
    )
    spaln_out.update(
        gffparser.add_missing_feats_to_gff3(spaln_out),
        merge_strategy="create_unique",
    )
    featuredb2gff3_file(spaln_out, sdir / "prot.gff3")


def _run_gth(sdir: Path) -> None:
    """Run GenomeThreader protein alignment (for intron-poor genomes)."""
    run_cmd(
        [
            "gth", "-genomic", "genome.gf", "-protein", "prots.faa",
            "-intermediate", "-xmlout", "-gzip", "-o", "gth.xml.gz", "-force",
        ],
        cwd=sdir,
    )
    run_cmd(
        [
            "gthconsensus", "-mincoverage", "0.70", "-minalignmentscore", "0.30",
            "-intermediate", "-gff3out", "-force", "-o", "gth.gff3", "gth.xml.gz",
        ],
        cwd=sdir,
    )
    gth_out = gffutils.create_db(
        str(sdir / "gth.gff3"), ":memory:",
        transform=gffparser.filter_gth_gff3,
        merge_strategy="create_unique",
    )
    gth_out.update(
        gffparser.add_missing_feats_to_gff3(gth_out),
        merge_strategy="create_unique",
    )
    gth_out = gffutils.create_db(
        gth_out, ":memory:", dialect=GFF3_DIALECT,
        transform=gffparser.fix_gth_ids,
    )
    featuredb2gff3_file(gth_out, sdir / "prot.gff3")
