"""AUGUSTUS gene prediction with training and parallel execution."""

from __future__ import annotations

import json
import os
import re
import shutil
from fileinput import FileInput
from pathlib import Path

import gffutils
from Bio import SeqIO

from eukan.annotation.training import build_training_set
from eukan.exceptions import ConfigurationError, ToolEnvError
from eukan.gff import create_gff_db
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.artifacts import Artifact, find_or_warn
from eukan.infra.concurrency import parallel_map
from eukan.infra.logging import get_logger
from eukan.infra.manifest import load_manifest
from eukan.infra.runner import run_cmd, run_piped
from eukan.infra.steps import step_dir
from eukan.infra.utils import concat_files, package_resource, symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


# Minimum evidence thresholds for accepting non-canonical splice sites.
# By default a type must pass BOTH the absolute floors AND the fraction of total
# junctions; --splice-permissive keeps the floors and drops only the fraction.
_MIN_JUNCTIONS = 10
_MIN_UNIQUE_READS = 50
_MIN_FRACTION = 0.01  # 1% of total junctions

# Splice site types that AUGUSTUS already allows by default — no need to add.
_AUGUSTUS_BUILTIN_SPLICE = {"gtag", "gcag"}

# Written into a species dir right after new_species.pl creates it. config.name
# defaults to the genome stem, which can be the name of a stock AUGUSTUS species
# (arabidopsis, human, fly, ...), so only a dir carrying this marker is replaced.
_SPECIES_MARKER = ".created-by-eukan"


def _splice_type_to_augustus(name: str) -> str | None:
    """Convert a splice site type name (e.g. ``GC-AG``) to a 4-char AUGUSTUS value.

    Returns ``None`` for entries that aren't a pair of A/C/G/T dinucleotides:
    ``"unknown"``, malformed names, and types with other bases such as ``NN-AC``
    (junctions in assembly gaps). AUGUSTUS accepts only a/c/g/t patterns and
    rejects the whole ``--allow_hinted_splicesites`` option if any one of them
    contains anything else.
    """
    if re.fullmatch(r"[ACGT]{2}-[ACGT]{2}", name.upper()) is None:
        return None
    return name.replace("-", "").lower()


def _get_splice_sites_flag(config: PipelineConfig) -> list[str]:
    """Build ``--allow_hinted_splicesites`` flag from RNA-seq splice junction evidence.

    Reads ``splice_site_summary.json`` (written by the assembly pipeline)
    to determine which non-canonical splice site types have sufficient
    evidence.  A type must pass absolute floors (junctions and unique reads)
    and, by default, a proportional threshold (fraction of total junctions);
    ``allow_noncanonical_splice`` (``--splice-permissive``) drops only the
    proportional threshold, so rare but well-supported types pass.
    All observed types that pass are included — not just the well-known
    semi-canonical pairs (GC-AG, AT-AC) — except types with bases other than
    A/C/G/T (e.g. ``NN-AC`` from junctions in assembly gaps), which AUGUSTUS
    cannot accept; those are dropped with a warning.

    Falls back to a blanket allowance when ``allow_noncanonical_splice``
    is set but no summary exists.
    """
    summary_path = find_or_warn(
        config.work_dir, Artifact.SPLICE_SUMMARY, load_manifest(config.manifest_dir),
    )
    allowed: set[str] = set()

    if summary_path is not None:
        with open(summary_path) as f:
            summary = json.load(f)

        total_junctions = sum(s["count"] for s in summary.values())

        # --splice-permissive drops only the share-of-all-junctions gate; the
        # absolute floors still keep singleton (often spurious) junctions out.
        min_fraction = 0.0 if config.allow_noncanonical_splice else _MIN_FRACTION

        non_acgt: dict[str, int] = {}
        for splice_type, stats in summary.items():
            count = stats["count"]
            fraction = count / total_junctions if total_junctions > 0 else 0.0
            if (
                count < _MIN_JUNCTIONS
                or stats["unique_reads"] < _MIN_UNIQUE_READS
                or fraction < min_fraction
            ):
                continue
            aug_name = _splice_type_to_augustus(splice_type)
            if aug_name is None:
                # "unknown" = junction outside the genome, not a base problem
                if splice_type != "unknown":
                    non_acgt[splice_type] = count
            elif aug_name not in _AUGUSTUS_BUILTIN_SPLICE:
                allowed.add(aug_name)

        if non_acgt:
            log.warning(
                "Ignoring %d splice-site types with non-ACGT bases "
                "(junctions in assembly gaps?): %s",
                len(non_acgt),
                ", ".join(f"{t} ({n} junctions)" for t, n in non_acgt.items()),
            )

    elif config.allow_noncanonical_splice:
        # No junction evidence — blanket allowance for common non-canonical types
        allowed = {"atac"}

    if allowed:
        value = ",".join(sorted(allowed))
        log.info("Allowing non-canonical splice sites in AUGUSTUS: %s", value)
        return [f"--allow_hinted_splicesites={value}"]

    return []


def run_augustus(config: PipelineConfig, *evidence: Path) -> Path:
    """Train and run AUGUSTUS gene prediction."""
    output = "augustus.gff3"
    sdir = step_dir(config.work_dir, "augustus")
    log.info("Running AUGUSTUS training and prediction...")

    # Gather hints. The protein hints sit beside the protein evidence (evidence[1],
    # prot.gff3 in the alignment step dir: prot_align, or prot_align_ssp with --spsp).
    hint_files = [evidence[1].parent / "hints_protein.gff"]
    if config.rnaseq_hints:
        ext_cfg = "augustus.config"
        hint_files.append(config.rnaseq_hints)
    else:
        ext_cfg = "extrinsic.MPE.cfg"

    # Auto-discover RepeatMasker hints from `eukan mask-repeats`. The RM
    # extrinsic source has weights in both our augustus.config and AUGUSTUS's
    # stock extrinsic.MPE.cfg, so the file just needs to be appended.
    rm_hints = find_or_warn(
        config.work_dir, Artifact.REPEATMASK_HINTS, load_manifest(config.manifest_dir),
    )
    if rm_hints is not None:
        log.info("Including RepeatMasker hints: %s", rm_hints.name)
        hint_files.append(rm_hints)

    # Create AUGUSTUS species
    config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if not config_path:
        raise ToolEnvError("augustus", env_var="AUGUSTUS_CONFIG_PATH")

    # Install custom extrinsic config if not already present
    ext_dest = Path(config_path) / "extrinsic" / ext_cfg
    if not ext_dest.exists():
        pkg_cfg = package_resource(f"configs/{ext_cfg}")
        if pkg_cfg is not None:
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
            if not (aug_species_dir / _SPECIES_MARKER).exists():
                raise ConfigurationError(
                    f"AUGUSTUS species '{config.name}' already exists in "
                    f"{aug_species_dir.parent} and was not created by eukan; "
                    "refusing to replace it.",
                    hint="The species name defaults to the genome file name. Train "
                    "under another name by setting EUKAN_NAME=<name> (there is no "
                    "command-line option for it) or by renaming the genome FASTA. "
                    "Delete the directory only if an earlier eukan version created it.",
                )
            shutil.rmtree(aug_species_dir)

        run_cmd(["new_species.pl", f"--species={config.name}"], cwd=sdir)
        (aug_species_dir / _SPECIES_MARKER).touch()

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

    concat_files(hint_files, sdir / "hints_all.gff")

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

    parallel_map(_run_split, splits, max_workers=config.num_cpu)

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

    aug_out = create_gff_db(sdir / "augustus.gff", transform=_augustus_clean)
    featuredb2gff3_file(aug_out, sdir / output)

    # Cleanup splits
    for split in splits:
        split.unlink(missing_ok=True)
        (sdir / f"{split.name}.gff3").unlink(missing_ok=True)

    return sdir / output


# ---------------------------------------------------------------------------
# AUGUSTUS output transforms
# ---------------------------------------------------------------------------


def _augustus_clean(f: gffutils.Feature) -> gffutils.Feature | None:
    """Keep only gene/mRNA/CDS/exon features from AUGUSTUS output."""
    if f.featuretype in ("gene", "mRNA", "CDS", "exon"):
        f.source = "augustus"
        return f
    return None
