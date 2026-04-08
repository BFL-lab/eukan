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
"""

from __future__ import annotations

import multiprocessing
import shutil
import sys
from pathlib import Path

# Ensure the project root is on sys.path so `from tests.testdata import ...` works
# regardless of how the script is invoked (poetry run, conda run, direct, etc.)
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import click

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
def cli(verbose: bool) -> None:
    """Pipeline integration test utilities."""
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
    Runs all steps of the annotation workflow:
      1. Database fetch (UniProt + Pfam for functional annotation)
      2. Transcriptome assembly (STAR mapping + Trinity + PASA)
      3. Genome annotation (GeneMark + spaln + AUGUSTUS + SNAP + EVM)
      4. Functional annotation (phmmer + hmmscan on predicted proteins)

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
    assembly_dir = work_dir / "assembly"
    annotation_dir = work_dir / "annotation"
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

    from eukan.functional.dbfetch import fetch_databases

    db_dir.mkdir(parents=True, exist_ok=True)
    uniprot_db = db_dir / "uniprot_sprot.faa"
    pfam_db = db_dir / "Pfam-A.hmm"

    if uniprot_db.exists() and pfam_db.exists():
        click.echo("  Databases already present, skipping download.")
    else:
        try:
            fetch_databases(db_dir)
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
            click.echo("Delete SRR*.fastq.gz and re-run: python tests/run_pipeline.py setup-test-data")
            protein_only = True

    if not protein_only:
        click.echo(f"\n{'=' * 60}")
        click.echo("STEP 2: Transcriptome assembly")
        click.echo(f"  {len(left_reads)} paired-end read pairs")
        click.echo(f"{'=' * 60}")

        from eukan.assembly import run_assembly
        from eukan.settings import AssemblyConfig

        # Sanitize genome headers before assembly
        from eukan.annotation.validation import sanitize_genome_fasta

        assembly_dir.mkdir(parents=True, exist_ok=True)
        genome = sanitize_genome_fasta(genome, assembly_dir)

        # Concatenate all forward and reverse reads for assembly
        concat_left = assembly_dir / "all_left.fastq.gz"
        concat_right = assembly_dir / "all_right.fastq.gz"

        if not concat_left.exists():
            click.echo("  Concatenating forward reads...")
            _concat_files(left_reads, concat_left)
        if not concat_right.exists():
            click.echo("  Concatenating reverse reads...")
            _concat_files(right_reads, concat_right)

        assembly_config = AssemblyConfig(
            genome=genome,
            work_dir=assembly_dir,
            left_reads=concat_left,
            right_reads=concat_right,
            num_cpu=numcpu,
            strand_specific="RF",
        )

        try:
            click.echo("  Running: read mapping (STAR)...")
            sys.stdout.flush()
            run_assembly(assembly_config, ["map"])
            click.echo("  Running: Trinity assembly (genome-guided + de novo)...")
            sys.stdout.flush()
            run_assembly(assembly_config, ["trinity"])
            click.echo("  Running: PASA alignment...")
            sys.stdout.flush()
            run_assembly(assembly_config, ["pasa"])
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

    from eukan.annotation import run_annotation_pipeline
    from eukan.settings import PipelineConfig

    annotation_dir.mkdir(parents=True, exist_ok=True)

    # Build config -- with or without transcript evidence
    config_kwargs: dict = {
        "genome": genome,
        "proteins": [proteins],
        "work_dir": annotation_dir,
        "kingdom": kingdom,
        "num_cpu": numcpu,
    }

    if not protein_only:
        # Check for assembly outputs
        nr_fasta = assembly_dir / "nr_transcripts.fasta"
        nr_gff3 = assembly_dir / "nr_transcripts.gff3"
        hints = assembly_dir / "hints_rnaseq.gff"

        if nr_fasta.exists() and nr_gff3.exists() and hints.exists():
            click.echo("  Using transcript evidence from assembly")
            config_kwargs["transcripts_fasta"] = nr_fasta
            config_kwargs["transcripts_gff"] = nr_gff3
            config_kwargs["rnaseq_hints"] = hints
        else:
            click.echo("  Assembly outputs not found, running without transcript evidence")
            missing = [f for f in [nr_fasta, nr_gff3, hints] if not f.exists()]
            for f in missing:
                click.echo(f"    missing: {f.name}")
    else:
        click.echo("  Running without transcript evidence (protein-only)")

    config = PipelineConfig(**config_kwargs)

    try:
        final_gff3 = run_annotation_pipeline(config)
        click.echo(f"\n  Annotation complete. Output: {final_gff3}")
    except Exception as e:
        click.echo(f"\n  Annotation failed: {e}")
        click.echo(f"\n  View run details: eukan status -d {annotation_dir}")
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
        from Bio import SeqIO

        from eukan.functional import run_functional_annotation
        from eukan.gff.io import extract_sequences

        # Extract predicted proteins from annotation
        predicted_proteins = annotation_dir / "predicted_proteins.faa"
        if not predicted_proteins.exists():
            click.echo("  Extracting predicted proteins from annotation...")
            records = list(extract_sequences(final_gff3, config.genome))
            SeqIO.write(records, str(predicted_proteins), "fasta")
            click.echo(f"  Extracted {len(records)} protein sequences.")
        else:
            num_seqs = sum(1 for _ in SeqIO.parse(str(predicted_proteins), "fasta"))
            click.echo(f"  Using existing predicted proteins ({num_seqs} sequences).")

        try:
            click.echo("  Running: phmmer (UniProt) + hmmscan (Pfam)...")
            sys.stdout.flush()
            run_functional_annotation(
                proteins=predicted_proteins,
                uniprot_db=uniprot_db,
                pfam_db=pfam_db,
                gff3_path=final_gff3,
                num_cpu=numcpu,
            )
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
    click.echo(f"  Annotation GFF3:  {final_gff3}")
    func_gff3 = final_gff3.with_suffix(".mod.gff3")
    if func_gff3.exists():
        click.echo(f"  Functional GFF3:  {func_gff3}")
    func_faa = annotation_dir / "predicted_proteins.mod.faa"
    if func_faa.exists():
        click.echo(f"  Annotated proteins: {func_faa}")
    click.echo(f"\n  View run details: eukan status -d {annotation_dir}")


@cli.command("compare-annotations", short_help="Compare pipeline output against reference.")
@click.option(
    "--reference", "-r", type=click.Path(exists=True, path_type=Path),
    default="tests/data/reference.gff3", show_default=True,
    help="Reference GFF3 file.",
)
@click.option(
    "--predicted", "-p", type=click.Path(exists=True, path_type=Path),
    default="tests/pipeline-run/annotation/final.mod.gff3", show_default=True,
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
    import os

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
