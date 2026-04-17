"""Spliced protein alignment via spaln (intron-rich) or GenomeThreader (intron-poor)."""

from __future__ import annotations

import os
import re
import shutil
import statistics
import time
from pathlib import Path

import gffutils
import pandas as pd

from eukan.annotation.validation import validate_fasta
from eukan.exceptions import ExternalToolError
from eukan.gff import create_gff_db, transform_db
from eukan.gff import parser as gffparser
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.runner import run_cmd, run_shell
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
    step_name = "prot_align_ssp" if config.spaln_ssp else "prot_align"
    existing = step_complete(config.work_dir, step_name, output)
    if existing:
        return existing

    sdir = step_dir(config.work_dir, step_name)
    log.info("Running protein alignment...")

    # Concatenate all protein files
    with open(sdir / "prots.faa", "w") as outfile:
        for prot_file in proteins:
            outfile.write(prot_file.read_text())

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
        ssp_name = _build_ssp(config, sdir)
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

    spaln_out = create_gff_db(sdir / "temp.gff3", transform=gffparser.Spaln.fix_cds_featuretype)
    spaln_out = transform_db(spaln_out, gffparser.Spaln.fix_ids)
    spaln_out.update(
        gffparser.add_missing_feats_to_gff3(spaln_out),
        merge_strategy="create_unique",
    )
    featuredb2gff3_file(spaln_out, sdir / "prot.gff3")


def _run_fitild(intron_hints_path: Path, sdir: Path) -> str:
    """Compute intron length distribution via fitild. Returns the -yI parameter string."""
    introns = pd.read_csv(intron_hints_path, sep="\t", header=None)
    intron_lens = (introns[4] - introns[3]).value_counts().reset_index().sort_values("index")
    intron_lens.to_csv(sdir / "introns.ild", sep=" ", header=None, index=False)

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


def _build_ssp(config: PipelineConfig, sdir: Path) -> str:
    """Build species-specific spaln parameters from transcript data.

    Creates a species parameter directory under $ALN_TAB and populates it
    using make_eij.pl + make_ssp.pl.  Returns the species name for spaln -T.

    The spaln makessp scripts expect genome databases named ``<species>_g``
    and transcript files named ``<species>_c.cf``, so we follow that
    convention using ``config.shortname`` as the species identifier.
    """
    ssp_name = config.shortname.lower()
    aln_tab = os.environ.get("ALN_TAB", "")
    if not aln_tab:
        raise RuntimeError(
            "ALN_TAB environment variable is not set. "
            "Ensure spaln is installed and eukan's environment is configured."
        )

    ssp_dir = Path(aln_tab) / ssp_name
    if ssp_dir.exists():
        shutil.rmtree(ssp_dir)
    ssp_dir.mkdir(parents=True)

    # spaln convention: genome = <ssp>_g, transcripts = <ssp>_c.cf
    genome_name = f"{ssp_name}_g"
    transcript_name = f"{ssp_name}_c"

    # Copy transcripts into the species parameter directory
    shutil.copy2(config.transcripts_fasta, ssp_dir / f"{transcript_name}.cf")

    # Symlink genome FASTA into ssp_dir and format for spaln
    symlink(config.genome, ssp_dir / f"{genome_name}.fna")
    run_cmd(
        ["spaln", "-W", "-KD", f"-t{config.num_cpu}", f"{genome_name}.fna"],
        cwd=ssp_dir,
    )

    log.info(
        "Building species-specific spaln parameters (%s) from %s",
        ssp_name, config.transcripts_fasta.name,
    )

    # Locate spaln perl scripts — prefer bundled versions in ALN_TAB/../alndbs
    # (which include Util.pm and match the compiled utility versions) over the
    # system-installed copies which may lag behind.
    # We invoke via "perl <script>" rather than executing the script directly
    # to ensure the conda perl is used, since the scripts' shebangs may
    # reference /usr/bin/perl (system perl) which can differ in version.
    alndbs = Path(aln_tab).parent / "alndbs"
    _make_eij = str(alndbs / "make_eij.pl") if (alndbs / "make_eij.pl").exists() else "make_eij.pl"
    _make_ssp = str(alndbs / "make_ssp.pl") if (alndbs / "make_ssp.pl").exists() else "make_ssp.pl"

    # Step 1: Map transcripts to genome, producing exon-intron junctions
    run_cmd(
        [
            "perl", _make_eij, f"-d{genome_name}",
            f"-t{config.num_cpu}", f"{transcript_name}.cf",
        ],
        cwd=ssp_dir,
    )

    eij_file = ssp_dir / f"{transcript_name}.eij"
    if not eij_file.exists():
        raise RuntimeError(
            f"make_eij.pl did not produce {eij_file}. "
            "Check that transcripts align to the genome."
        )

    # Count introns (header line starts with #)
    with open(eij_file) as fh:
        n_introns = sum(1 for line in fh if not line.startswith("#"))
    if n_introns < 5000:
        log.warning(
            "Only %d unique introns in %s (recommend >= 5000 for good parameters). "
            "Species-specific parameters may be unreliable.",
            n_introns, eij_file.name,
        )

    # Ensure file modification times are distinct — make_ssp.pl uses
    # Perl's -M operator to compare input/output ages and may skip
    # generation when files are created within the same second.
    time.sleep(1)

    # Step 2: Generate species-specific parameter files.
    # When -S is set, make_ssp.pl runs levels 9, 13, and 16.  Level 16
    # requires the fitild binary; if fitild is absent the script exits 1
    # but the critical output (Splice3/5, IntronPotTab) from levels 9 and
    # 13 has already been written.  We therefore tolerate a non-zero exit
    # code as long as the required files were produced.
    # Run via shell to ensure Perl's system() calls and shell redirects
    # inside make_ssp.pl work correctly.
    try:
        run_shell(
            f"perl {_make_ssp} -d{genome_name} -S {transcript_name}.eij",
            cwd=ssp_dir,
        )
    except ExternalToolError:
        pass  # check for output files below

    # Verify key output files were created (may be .dat, .dgz, or extensionless)
    for base in ("Splice3", "Splice5", "IntronPotTab"):
        if not any((ssp_dir / f"{base}{ext}").exists() for ext in ("", ".dat", ".dgz")):
            raise RuntimeError(
                f"make_ssp.pl did not produce {base} in {ssp_dir}. "
                "Species-specific parameter generation may have failed."
            )

    log.info("Species-specific spaln parameters written to %s", ssp_dir)
    return ssp_name


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
    gth_out = create_gff_db(sdir / "gth.gff3", transform=gffparser.Gth.filter)
    gth_out.update(
        gffparser.add_missing_feats_to_gff3(gth_out),
        merge_strategy="create_unique",
    )
    gth_out = transform_db(gth_out, gffparser.Gth.fix_ids)
    featuredb2gff3_file(gth_out, sdir / "prot.gff3")
