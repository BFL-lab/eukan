"""Spliced protein alignment via spaln (intron-rich) or GenomeThreader (intron-poor)."""

from __future__ import annotations

import re
import statistics
from collections import Counter
from pathlib import Path

import gffutils

from eukan.annotation.ssp import build_ssp
from eukan.annotation.validation import validate_fasta
from eukan.gff import create_gff_db
from eukan.gff import parser as gffparser
from eukan.gff.normalize import normalize_to_gff3
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.infra.utils import concat_files, symlink
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
    step_name = "prot_align_ssp" if config.spaln_ssp else "prot_align"
    sdir = step_dir(config.work_dir, step_name)
    log.info("Running protein alignment...")

    concat_files(proteins, sdir / "prots.faa")

    validate_fasta(sdir / "prots.faa")
    symlink(config.genome, sdir / "genome")
    symlink(config.genome, sdir / "genome.gf")

    models = create_gff_db(gff3)

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
    prot_db = create_gff_db(sdir / output, transform=gffparser.prot2augustus_hints)
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

    use_ssp = config.spaln_ssp and config.transcripts_fasta
    if use_ssp:
        ssp_name = build_ssp(config, sdir)
    else:
        int_dist = _run_fitild(intron_hints_path, sdir)

    run_cmd(
        ["spaln", "-Wgenome.bkp", f"-XG{maxlen}", "-KP", "genome.gf"],
        cwd=sdir,
    )

    # Remove stale spaln output files to avoid interactive overwrite prompts
    for ext in (".grd", ".erd", ".qrd"):
        (sdir / f"prots{ext}").unlink(missing_ok=True)

    if use_ssp:
        run_cmd(
            [
                "spaln", "-O12", "-Q7", "-LS",
                f"-C{config.genetic_code}", f"-t{config.num_cpu}",
                f"-T{ssp_name}", "-dgenome", "prots.faa",
            ],
            cwd=sdir,
        )
    else:
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

    normalize_to_gff3(
        sdir / "temp.gff3", sdir / "prot.gff3",
        parse_transform=gffparser.Spaln.fix_cds_featuretype,
        post_transform=gffparser.Spaln.fix_ids,
    )


def _run_fitild(intron_hints_path: Path, sdir: Path) -> str:
    """Compute intron length distribution via fitild. Returns the -yI parameter string."""
    counts: Counter[int] = Counter()
    with open(intron_hints_path) as fin:
        for line in fin:
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            try:
                length = int(cols[4]) - int(cols[3])
            except ValueError:
                continue
            counts[length] += 1
    with open(sdir / "introns.ild", "w") as fout:
        for length in sorted(counts):
            fout.write(f"{length} {counts[length]}\n")

    run_cmd(["fitild", "-d", "IldModel.txt", "introns.ild"], cwd=sdir, out_file="fitild.out")

    int_dist_parts: list[str] = []
    with open(sdir / "fitild.out") as fp:
        for line in fp:
            if re.search(r"^introns\s", line):
                parts = line.replace("\t", " ").replace("\n", " ").split()[3:-3]
                int_dist_parts = [f"{float(x):f}" for x in parts]
    if int_dist_parts:
        del int_dist_parts[1]
    return " ".join(int_dist_parts)


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
    normalize_to_gff3(
        sdir / "gth.gff3", sdir / "prot.gff3",
        parse_transform=gffparser.Gth.filter,
        post_transform=gffparser.Gth.fix_ids,
    )
