#!/usr/bin/env python3
"""Compare spaln fitild vs species-specific parameter alignment results.

Runs protein alignment in both modes on the existing test pipeline data
and prints a side-by-side comparison of the output GFF3 files.

Usage:
    conda run -n eukan python tests/compare_spaln_modes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from eukan.infra.environ import configure_process_env
from eukan.infra.logging import setup_logging

setup_logging(0)
configure_process_env()

import gffutils

from eukan.annotation.alignment import align_proteins
from eukan.settings import PipelineConfig

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "tests" / "data"
PIPELINE_RUN = BASE / "tests" / "pipeline-run"
ASSEMBLY = PIPELINE_RUN / "assembly"
ANNOTATION = PIPELINE_RUN / "annotation"
PROTEINS = DATA / "proteins.faa"


def count_features(gff3: Path) -> dict[str, int]:
    """Count features by type in a GFF3 file."""
    db = gffutils.create_db(str(gff3), ":memory:", merge_strategy="create_unique")
    counts: dict[str, int] = {}
    for feat in db.all_features():
        counts[feat.featuretype] = counts.get(feat.featuretype, 0) + 1
    return counts


def run_mode(spaln_ssp: bool) -> Path:
    """Run protein alignment in the specified mode and return the output path."""
    mode = "ssp" if spaln_ssp else "fitild"
    print(f"\n{'='*60}")
    print(f"Running protein alignment in {mode} mode...")
    print(f"{'='*60}")

    config = PipelineConfig(
        genome=ASSEMBLY / "genome.sanitized.fasta",
        proteins=[PROTEINS],
        work_dir=ANNOTATION,
        num_cpu=4,
        spaln_ssp=spaln_ssp,
        transcripts_fasta=ASSEMBLY / "nr_transcripts.fasta" if spaln_ssp else None,
    )
    genemark_gff = ANNOTATION / "genemark" / "genemark.gtf"
    intron_hints = ANNOTATION / "genemark" / "introns.gff"

    result = align_proteins(
        config,
        genemark_gff,
        config.proteins,
        intron_hints if intron_hints.exists() else None,
    )
    print(f"  Output: {result}")
    return result


def compare(fitild_gff3: Path, ssp_gff3: Path) -> None:
    """Print a comparison of the two GFF3 outputs."""
    fitild_counts = count_features(fitild_gff3)
    ssp_counts = count_features(ssp_gff3)

    all_types = sorted(set(fitild_counts) | set(ssp_counts))

    print(f"\n{'='*60}")
    print("Feature comparison: fitild vs species-specific params")
    print(f"{'='*60}")
    print(f"  {'Feature type':<15} {'fitild':>10} {'ssp':>10} {'diff':>10}")
    print(f"  {'-'*45}")
    for ft in all_types:
        f = fitild_counts.get(ft, 0)
        s = ssp_counts.get(ft, 0)
        diff = s - f
        sign = "+" if diff > 0 else ""
        print(f"  {ft:<15} {f:>10} {s:>10} {sign}{diff:>9}")

    print(f"\n  fitild output: {fitild_gff3}")
    print(f"  ssp output:    {ssp_gff3}")


if __name__ == "__main__":
    # The fitild result should already exist from a previous pipeline run
    fitild_gff3 = ANNOTATION / "prot_align" / "prot.gff3"
    if not fitild_gff3.exists():
        print("No existing fitild result found. Running fitild mode first...")
        fitild_gff3 = run_mode(spaln_ssp=False)
    else:
        print(f"Using existing fitild result: {fitild_gff3}")

    # Run the ssp mode
    ssp_gff3 = run_mode(spaln_ssp=True)

    # Compare
    compare(fitild_gff3, ssp_gff3)
