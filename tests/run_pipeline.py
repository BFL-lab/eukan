#!/usr/bin/env python3
"""Development CLI for pipeline integration testing.

Usage:
    python tests/run_pipeline.py setup-test-data [-o tests/data]
    python tests/run_pipeline.py test-pipeline [--kingdom fungus] [-n 8]
    python tests/run_pipeline.py clean-test-data [--all]

Full pipeline order:
    1. Database fetch (UniProt + Pfam)
    2. Transcriptome assembly (STAR + Trinity + PASA)
    3. Genome annotation (GeneMark + spaln + AUGUSTUS + SNAP + EVM)
    4. Functional annotation (phmmer + hmmscan on predicted proteins)

The test-pipeline command invokes `eukan` CLI subcommands via subprocess,
exercising the same code paths as real user workflows.
"""

from __future__ import annotations

import multiprocessing
import os
import shutil
import subprocess
import sys
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import click

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def _run_eukan(args: list[str], cwd: Path, label: str) -> subprocess.CompletedProcess:
    """Run an eukan CLI command via subprocess, mirroring real user usage."""
    cmd = ["eukan"] + args
    click.echo(f"  $ {' '.join(cmd)}")
    sys.stdout.flush()
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {result.returncode})")
    return result


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
def cli(verbose: bool) -> None:
    """Pipeline integration test utilities."""
    # Environment setup only needed for setup-test-data and clean-test-data
    # (test-pipeline shells out to eukan which does its own setup)
    from eukan.infra.environ import configure_process_env
    from eukan.infra.logging import setup_logging
    setup_logging(1 if verbose else 0)
    configure_process_env()


@cli.command("setup-test-data", short_help="Download S. pombe test data from NCBI.")
@click.option(
    "--output-dir", "-o", type=click.Path(path_type=Path), default="tests/data",
    show_default=True, help="Directory to download test data into.",
)
def setup_test_data_cmd(output_dir: Path) -> None:
    """Download S. pombe test data from NCBI for pipeline testing.

    Downloads:
    - Genome: S. pombe chromosome III (NC_003424.3)
    - Proteins: 10 close neighbor proteomes
    - RNA-seq: 5 SRA paired-end runs

    Requires NCBI datasets CLI and SRA Toolkit on PATH.
    """
    from tests.testdata import setup_test_data

    setup_test_data(output_dir.resolve())


@cli.command("test-pipeline", short_help="Run the full pipeline on S. pombe test data.")
@click.option(
    "--data-dir", "-d", type=click.Path(path_type=Path), default="tests/data",
    show_default=True, help="Directory containing test data.",
)
@click.option(
    "--work-dir", "-w", type=click.Path(path_type=Path), default="tests/pipeline-run",
    show_default=True, help="Working directory for the pipeline run.",
)
@click.option(
    "--kingdom", "-k", type=click.Choice(["fungus", "protist", "animal", "plant"]),
    default="fungus", show_default=True,
    help="Kingdom for test organism (S. pombe = fungus).",
)
@click.option(
    "--numcpu", "-n", type=int, default=multiprocessing.cpu_count(),
    show_default=True, help="Number of CPU threads.",
)
@click.option(
    "--protein-only", is_flag=True,
    help="Skip assembly, annotate using only protein + genome evidence.",
)
def test_pipeline_cmd(
    data_dir: Path, work_dir: Path, kingdom: str, numcpu: int, protein_only: bool,
) -> None:
    """Run the full pipeline on S. pombe test data.

    \b
    Runs all steps of the annotation workflow via the eukan CLI:
      1. eukan db-fetch        Download UniProt + Pfam databases
      2. eukan assemble        STAR mapping + Trinity + PASA
      3. eukan annotate        GeneMark + spaln + AUGUSTUS + SNAP + EVM
      4. eukan func-annot      phmmer + hmmscan on predicted proteins

    Requires test data from: python tests/run_pipeline.py setup-test-data
    """
    from tests.testdata import validate_test_data

    data_dir = data_dir.resolve()
    work_dir = work_dir.resolve()

    # --- Validate test data ---
    click.echo("Validating test data...")
    results = validate_test_data(data_dir)
    failures = [r for r in results if not r[2]]
    if failures:
        for name, msg, _ in failures:
            click.echo(f"  \u2717 {name}: {msg}")
        click.echo("\nRun `python tests/run_pipeline.py setup-test-data` first.")
        raise SystemExit(1)
    for name, msg, _ in results:
        click.echo(f"  \u2713 {name}: {msg}")

    # --- Locate test files ---
    genome = data_dir / "genome.fasta"
    proteins = data_dir / "proteins.faa"
    left_reads = sorted(data_dir.glob("SRR*_1.fastq.gz"))
    right_reads = sorted(data_dir.glob("SRR*_2.fastq.gz"))

    work_dir.mkdir(parents=True, exist_ok=True)
    db_dir = work_dir / "databases"

    click.echo(f"\nWork directory: {work_dir}")
    click.echo(f"Genome: {genome}")
    click.echo(f"Proteins: {proteins}")
    click.echo(f"Kingdom: {kingdom}")
    click.echo(f"CPUs: {numcpu}")

    # ================================================================
    # Step 1: Database fetch
    # ================================================================
    click.echo(f"\n{'=' * 60}")
    click.echo("STEP 1: Database fetch (UniProt + Pfam)")
    click.echo(f"{'=' * 60}")

    uniprot_db = db_dir / "uniprot_sprot.faa"
    pfam_db = db_dir / "Pfam-A.hmm"

    if uniprot_db.exists() and pfam_db.exists():
        click.echo("  Databases already present, skipping download.")
    else:
        try:
            _run_eukan(
                ["db-fetch", "-o", str(db_dir)],
                cwd=work_dir, label="Database fetch",
            )
            click.echo("  Database fetch complete.")
        except Exception as e:
            click.echo(f"\n  Database fetch failed: {e}")
            click.echo("  Functional annotation will be skipped.")

    # ================================================================
    # Step 2: Transcriptome assembly
    # ================================================================
    if not protein_only:
        if not left_reads or not right_reads:
            click.echo("\nNo paired-end reads (_1/_2) found. Skipping assembly.")
            protein_only = True

    if not protein_only:
        click.echo(f"\n{'=' * 60}")
        click.echo("STEP 2: Transcriptome assembly")
        click.echo(f"  {len(left_reads)} paired-end read pairs")
        click.echo(f"{'=' * 60}")

        # Concatenate all forward and reverse reads for assembly
        concat_left = work_dir / "all_left.fastq.gz"
        concat_right = work_dir / "all_right.fastq.gz"

        if not concat_left.exists():
            click.echo("  Concatenating forward reads...")
            _concat_files(left_reads, concat_left)
        if not concat_right.exists():
            click.echo("  Concatenating reverse reads...")
            _concat_files(right_reads, concat_right)

        try:
            _run_eukan(
                [
                    "assemble",
                    "-g", str(genome),
                    "-l", str(concat_left),
                    "-r", str(concat_right),
                    "-S", "RF",
                    "-n", str(numcpu),
                ],
                cwd=work_dir, label="Assembly",
            )
            click.echo("  Assembly complete.")
        except Exception as e:
            click.echo(f"\n  Assembly failed: {e}")
            click.echo("  Continuing to annotation without transcript evidence...")
            protein_only = True

    # ================================================================
    # Step 3: Genome annotation
    # ================================================================
    sys.stdout.flush()
    click.echo(f"\n{'=' * 60}")
    click.echo("STEP 3: Genome annotation")
    click.echo(f"{'=' * 60}")

    annotate_args = [
        "annotate",
        "-g", str(genome),
        "-p", str(proteins),
        "-k", kingdom,
        "-n", str(numcpu),
    ]

    if not protein_only:
        # Check for assembly outputs (written to work_dir by eukan assemble)
        nr_fasta = work_dir / "nr_transcripts.fasta"
        nr_gff3 = work_dir / "nr_transcripts.gff3"
        hints = work_dir / "hints_rnaseq.gff"

        if nr_fasta.exists() and nr_gff3.exists() and hints.exists():
            click.echo("  Using transcript evidence from assembly")
            annotate_args += [
                "-tf", str(nr_fasta),
                "-tg", str(nr_gff3),
                "-r", str(hints),
                "--strand-specific",
            ]
        else:
            click.echo("  Assembly outputs not found, running without transcript evidence")
            missing = [f for f in [nr_fasta, nr_gff3, hints] if not f.exists()]
            for f in missing:
                click.echo(f"    missing: {f.name}")
    else:
        click.echo("  Running without transcript evidence (protein-only)")

    try:
        _run_eukan(annotate_args, cwd=work_dir, label="Annotation")
        click.echo("  Annotation complete.")
    except Exception as e:
        click.echo(f"\n  Annotation failed: {e}")
        click.echo(f"\n  View run details: eukan status -d {work_dir}")
        raise SystemExit(1)

    # ================================================================
    # Step 4: Functional annotation
    # ================================================================
    sys.stdout.flush()
    click.echo(f"\n{'=' * 60}")
    click.echo("STEP 4: Functional annotation")
    click.echo(f"{'=' * 60}")

    if not uniprot_db.exists() or not pfam_db.exists():
        click.echo("  Skipping: databases not available.")
    else:
        # Extract predicted proteins if not already done
        predicted_proteins = work_dir / "predicted_proteins.faa"
        final_gff3 = work_dir / "final.gff3"

        if not predicted_proteins.exists() and final_gff3.exists():
            click.echo("  Extracting predicted proteins from annotation...")
            # gff3toseq writes to stdout; capture and save to file
            cmd = ["eukan", "gff3toseq", "-g", str(genome), "-i", str(final_gff3), "-o", "protein"]
            click.echo(f"  $ {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
            if result.returncode == 0 and result.stdout:
                predicted_proteins.write_text(result.stdout)
                n_seqs = result.stdout.count(">")
                click.echo(f"  Extracted {n_seqs} protein sequences.")
            else:
                click.echo("  Failed to extract proteins, skipping functional annotation.")
                predicted_proteins = None

        if predicted_proteins and predicted_proteins.exists():
            func_args = [
                "func-annot",
                "-p", str(predicted_proteins),
                "--uniprot", str(uniprot_db),
                "--pfam", str(pfam_db),
                "--gff3", str(final_gff3),
                "-n", str(numcpu),
            ]
            try:
                _run_eukan(func_args, cwd=work_dir, label="Functional annotation")
                click.echo("  Functional annotation complete.")
            except Exception as e:
                click.echo(f"\n  Functional annotation failed: {e}")
                raise SystemExit(1)

    # ================================================================
    # Summary
    # ================================================================
    click.echo(f"\n{'=' * 60}")
    click.echo("Pipeline complete.")
    click.echo(f"{'=' * 60}")
    final_gff3 = work_dir / "final.gff3"
    if final_gff3.exists():
        click.echo(f"  Annotation GFF3:  {final_gff3}")
    func_gff3 = work_dir / "final.mod.gff3"
    if func_gff3.exists():
        click.echo(f"  Functional GFF3:  {func_gff3}")
    func_faa = work_dir / "predicted_proteins.mod.faa"
    if func_faa.exists():
        click.echo(f"  Annotated proteins: {func_faa}")
    click.echo(f"\n  View run details: eukan status")


@cli.command("compare-annotations", short_help="Compare pipeline output against reference.")
@click.option(
    "--reference", "-r", type=click.Path(exists=True, path_type=Path),
    default="tests/data/reference.gff3", show_default=True,
    help="Reference GFF3 file.",
)
@click.option(
    "--predicted", "-p", type=click.Path(exists=True, path_type=Path),
    default="tests/pipeline-run/final.mod.gff3", show_default=True,
    help="Predicted GFF3 file to evaluate.",
)
def compare_annotations_cmd(reference: Path, predicted: Path) -> None:
    """Compare predicted gene models against reference annotations.

    \b
    Gene level: exact / inexact / missing / merged / fragmented
      (merged/fragmented use 50% overlap thresholds)
    mRNA / CDS / Intron levels: match / missing / FP
      (maximum pairwise overlap matching within matched parents)
    Metrics: sensitivity, precision, F1 (count- and overlap-based)

    Defaults compare tests/data/reference.gff3 vs the pipeline output.
    """
    from tests.annot_quality import compare_annotations, format_results

    click.echo(f"Loading reference: {reference}")
    click.echo(f"Loading predicted: {predicted}")
    click.echo()

    result = compare_annotations(reference.resolve(), predicted.resolve())
    click.echo(format_results(result))


@cli.command("clean-test-data")
@click.option(
    "--data-dir", "-d", type=click.Path(path_type=Path), default="tests/data",
    show_default=True, help="Directory containing downloaded test data.",
)
@click.option(
    "--work-dir", "-w", type=click.Path(path_type=Path), default="tests/pipeline-run",
    show_default=True, help="Pipeline run working directory.",
)
@click.option(
    "--all", "clean_all", is_flag=True,
    help="Also remove downloaded test data (genome, proteins, reads).",
)
def clean_test_data_cmd(data_dir: Path, work_dir: Path, clean_all: bool) -> None:
    """Remove pipeline test outputs and optionally downloaded data.

    \b
    By default, removes only the pipeline run directory (tests/pipeline-run),
    including assembly, annotation, and database outputs.
    With --all, also removes downloaded genome, proteins, and FASTQ files.
    Accession list files (tests/data/*.txt) are never deleted.
    """
    data_dir = data_dir.resolve()
    work_dir = work_dir.resolve()

    # Always clean pipeline run directory
    if work_dir.exists():
        # Some tools (Trinity) create dirs with restricted permissions
        for root, dirs, files in os.walk(work_dir):
            for d in dirs:
                os.chmod(os.path.join(root, d), 0o755)
        shutil.rmtree(work_dir)
        click.echo(f"Removed {work_dir}")
    else:
        click.echo(f"Nothing to clean: {work_dir} does not exist")

    # Optionally clean downloaded data (but keep accession .txt files)
    if clean_all:
        patterns = ["genome.fasta", "proteins.faa", "SRR*.fastq.gz", "SRR*.fastq", "SRR*.sra"]
        removed = 0
        for pattern in patterns:
            for f in data_dir.glob(pattern):
                f.unlink()
                removed += 1
        click.echo(f"Removed {removed} downloaded files from {data_dir}")


def _concat_files(inputs: list[Path], output: Path) -> None:
    """Concatenate multiple files into one (binary, for gzipped FASTQs)."""
    with open(output, "wb") as out:
        for f in inputs:
            with open(f, "rb") as inp:
                shutil.copyfileobj(inp, out)


if __name__ == "__main__":
    cli()
