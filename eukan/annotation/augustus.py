"""AUGUSTUS gene prediction with training and parallel execution."""

from __future__ import annotations

import json
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
from fileinput import FileInput
from pathlib import Path

import gffutils
from Bio import SeqIO

from eukan.exceptions import ToolEnvError
from eukan.gff import create_gff_db
from eukan.gff import intersecter as gffintersecter
from eukan.gff import parser as gffparser
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_piped
from eukan.infra.steps import step_dir
from eukan.infra.utils import symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


# Minimum evidence thresholds for accepting non-canonical splice sites.
# A type must pass BOTH the absolute count AND the fraction of total junctions.
_MIN_JUNCTIONS = 10
_MIN_UNIQUE_READS = 50
_MIN_FRACTION = 0.01  # 1% of total junctions

# Splice site types that AUGUSTUS already allows by default — no need to add.
_AUGUSTUS_BUILTIN_SPLICE = {"gtag", "gcag"}


def _splice_type_to_augustus(name: str) -> str | None:
    """Convert a splice site type name (e.g. ``GC-AG``) to a 4-char AUGUSTUS value.

    Returns ``None`` for entries that don't look like valid dinucleotide pairs
    (e.g. ``"unknown"``).
    """
    parts = name.split("-")
    if len(parts) == 2 and len(parts[0]) == 2 and len(parts[1]) == 2:
        return (parts[0] + parts[1]).lower()
    return None


def _get_splice_sites_flag(config: PipelineConfig) -> list[str]:
    """Build ``--allow_hinted_splicesites`` flag from STAR splice site evidence.

    Reads ``splice_site_summary.json`` (written by the assembly pipeline)
    to determine which non-canonical splice site types have sufficient
    evidence.  A type must pass both an absolute count threshold and a
    proportional threshold (fraction of total junctions) to be accepted.
    All observed types that pass are included — not just the well-known
    semi-canonical pairs (GC-AG, AT-AC).

    Falls back to a blanket allowance when ``allow_noncanonical_splice``
    is set but no summary exists.
    """
    summary_path = config.work_dir / "splice_site_summary.json"
    allowed: set[str] = set()

    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)

        total_junctions = sum(s["count"] for s in summary.values())

        # Use lower thresholds when --splice-permissive is set
        if config.allow_noncanonical_splice:
            min_junctions = 1
            min_reads = 1
            min_fraction = 0.0
        else:
            min_junctions = _MIN_JUNCTIONS
            min_reads = _MIN_UNIQUE_READS
            min_fraction = _MIN_FRACTION

        for splice_type, stats in summary.items():
            count = stats["count"]
            fraction = count / total_junctions if total_junctions > 0 else 0.0
            if count < min_junctions or stats["unique_reads"] < min_reads or fraction < min_fraction:
                continue
            aug_name = _splice_type_to_augustus(splice_type)
            if aug_name and aug_name not in _AUGUSTUS_BUILTIN_SPLICE:
                allowed.add(aug_name)

    elif config.allow_noncanonical_splice:
        # No STAR evidence — blanket allowance for common non-canonical types
        allowed = {"atac"}

    if allowed:
        value = ",".join(sorted(allowed))
        log.info("Allowing non-canonical splice sites in AUGUSTUS: %s", value)
        return [f"--allow_hinted_splicesites={value}"]

    return []


def build_training_set(
    config: PipelineConfig, evidence: tuple[Path, ...], sdir: Path
) -> None:
    """Build evidence-driven training set for gene predictors."""
    training_db = gffintersecter.extract_supported_models(*evidence, output_dir=sdir)

    num_training = max(
        1, len([f.id for f in training_db.features_of_type("gene")]) // 4
    )

    gene_db = gffutils.create_db(str(evidence[0]), ":memory:")
    gene_lengths = [f.end - f.start for f in gene_db.features_of_type("gene")]
    flank = int(sum(gene_lengths) / len(gene_lengths) / 2)

    run_cmd(
        ["gff2gbSmallDNA.pl", "training_set.gff3", str(config.genome), str(flank), "genbank.gb"],
        cwd=sdir,
    )
    run_cmd(["randomSplit.pl", "genbank.gb", str(num_training)], cwd=sdir)


def run_augustus(config: PipelineConfig, *evidence: Path) -> Path:
    """Train and run AUGUSTUS gene prediction."""
    output = "augustus.gff3"
    sdir = step_dir(config.work_dir, "augustus")
    log.info("Running AUGUSTUS training and prediction...")

    # Gather hints
    hint_files = [config.work_dir / "prot_align" / "hints_protein.gff"]
    if config.rnaseq_hints:
        ext_cfg = "augustus.config"
        hint_files.append(config.rnaseq_hints)
    else:
        ext_cfg = "extrinsic.MPE.cfg"

    # Create AUGUSTUS species
    config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if not config_path:
        raise ToolEnvError("augustus", env_var="AUGUSTUS_CONFIG_PATH")

    # Install custom extrinsic config if not already present
    ext_dest = Path(config_path) / "extrinsic" / ext_cfg
    if not ext_dest.exists():
        pkg_cfg = Path(__file__).parent.parent.parent / "configs" / ext_cfg
        if pkg_cfg.exists():
            shutil.copy(pkg_cfg, ext_dest)
            log.info("Installed and using extrinsic config: %s", ext_cfg)
        else:
            log.warning("Custom extrinsic config %s not found, AUGUSTUS may fail", ext_cfg)

    aug_species_dir = Path(config_path) / "species" / config.name
    aug_cfg_file = aug_species_dir / f"{config.name}_parameters.cfg"

    # Skip training if species already exists and has been optimized
    if aug_cfg_file.exists() and (sdir / "genbank.gb.train").exists():
        log.info("[augustus] Species %s already trained, skipping to prediction.", config.name)
    else:
        if aug_species_dir.exists():
            shutil.rmtree(aug_species_dir)

        run_cmd(["new_species.pl", f"--species={config.name}"], cwd=sdir)

        # Build training set and train
        build_training_set(config, evidence, sdir)

        # Report training set size
        with open(sdir / "genbank.gb.train") as fh:
            train_count = sum(1 for line in fh if line.startswith("LOCUS"))
        with open(sdir / "genbank.gb.test") as fh:
            test_count = sum(1 for line in fh if line.startswith("LOCUS"))
        total = train_count + test_count
        if total <= 500:
            log.warning(
                "Low AUGUSTUS training set size: %d models (%d train, %d test) "
                "— prediction quality may be reduced",
                total, train_count, test_count,
            )
        else:
            log.info("AUGUSTUS training set: %d train, %d test models", train_count, test_count)

        run_cmd(["etraining", f"--species={config.name}", "genbank.gb.train"], cwd=sdir, out_file="etraining.out")

        # Set stop codon frequencies
        if config.genetic_code == "6":
            stop_probs = ["0", "0", "1"]
        else:
            stop_probs = []
            with open(sdir / "etraining.out") as fp:
                for line in fp:
                    if re.search(r"(tag:|taa:|tga:)", line):
                        stop_probs.append(re.sub(r"[()]", "", line).split()[-1])

        with FileInput(files=[str(aug_cfg_file)], inplace=True) as f:
            for line in f:
                line = line.rstrip()
                if re.search(r"(amberprob|ochreprob|opalprob)", line) and stop_probs:
                    parts = line.split()
                    parts[1] = stop_probs.pop(0)
                    print(" ".join(parts))
                else:
                    print(line)

        # Test and optimize
        run_cmd(["augustus", f"--species={config.name}", "genbank.gb.test"], cwd=sdir)
        run_cmd(
            [
                "optimize_augustus.pl", f"--species={config.name}",
                "--onlytrain=genbank.gb.train", f"--cpus={config.num_cpu}",
                "genbank.gb.test",
            ],
            cwd=sdir,
        )
        run_cmd(["etraining", f"--species={config.name}", "genbank.gb.train"], cwd=sdir)
        run_cmd(["augustus", f"--species={config.name}", "genbank.gb.test"], cwd=sdir)

    # Merge hints
    with open(sdir / "hints_all.gff", "w") as outfile:
        for fname in hint_files:
            with open(fname) as infile:
                outfile.write(infile.read())

    # Split genome for parallel prediction
    assembly_size = sum(len(rec) for rec in SeqIO.parse(str(config.genome), "fasta"))
    min_size = max(1, assembly_size // (config.num_cpu * 2))
    symlink(config.genome, sdir / "genome.fa")
    run_cmd(["splitMfasta.pl", "genome.fa", f"--minsize={min_size}"], cwd=sdir)

    splits = list(sdir.glob("genome.split.*.fa"))

    # Determine which non-canonical splice sites to allow based on evidence
    splice_sites_flag = _get_splice_sites_flag(config)

    base_cmd = [
        "augustus", f"--species={config.name}",
        f"--extrinsicCfgFile={Path(config_path) / 'extrinsic' / ext_cfg}",
        "--hintsfile=hints_all.gff", "--softmasking=1", "--UTR=off",
        *splice_sites_flag,
    ]

    # Run splits -- each writes stdout to its own .gff3 file
    def _run_split(split: Path) -> None:
        run_cmd(
            [*base_cmd, str(split)],
            cwd=sdir, out_file=f"{split.name}.gff3",
        )

    with ThreadPoolExecutor(max_workers=config.num_cpu) as pool:
        list(pool.map(_run_split, splits))

    # Join predictions
    cat_files = sorted(sdir.glob("genome.split.*.gff3"))

    run_piped(
        ["cat"] + [str(f) for f in cat_files],
        ["join_aug_pred.pl"],
        cwd=sdir, out_file="joined.gff",
    )
    run_piped(
        ["cat", "joined.gff"],
        ["gtf2gff.pl", "--printExon", "-gff3", f"--out={sdir / 'augustus.gff'}"],
        cwd=sdir,
    )

    aug_out = create_gff_db(sdir / "augustus.gff", transform=gffparser.Augustus.clean)
    featuredb2gff3_file(aug_out, sdir / output)

    # Cleanup splits
    for split in splits:
        split.unlink(missing_ok=True)
        (sdir / f"{split.name}.gff3").unlink(missing_ok=True)

    return sdir / output
