"""Annotation pipeline orchestration: step ordering, concurrency, and manifest tracking."""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from eukan.annotation.alignment import align_proteins
from eukan.annotation.augustus import run_augustus
from eukan.annotation.consensus import build_consensus_models
from eukan.annotation.genemark import run_genemark
from eukan.annotation.snap import run_codingquarry, run_snap
from eukan.annotation.validation import sanitize_genome_fasta, validate_fasta
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import count_gff3_features, get_logger
from eukan.infra.manifest import (
    RunManifest, init_manifest, load_manifest, save_manifest,
    pipeline_step, is_step_complete, is_step_interrupted, clean_interrupted_step,
)
from eukan.infra.steps import step_complete, step_dir
from eukan.annotation.orf import create_transcriptome_orf_db
from eukan.settings import PipelineConfig

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# ORF finding step
# ---------------------------------------------------------------------------


def find_orfs(config: PipelineConfig, trans_gff3: Path) -> Path:
    """Find ORFs in transcript assemblies."""
    from eukan.annotation.validation import validate_gff

    output = "transcript_orfs.gff3"
    existing = step_complete(config.work_dir, "orf_finder", output)
    if existing:
        return existing

    sdir = step_dir(config.work_dir, "orf_finder")
    log.info("Finding ORFs in transcript assemblies...")
    validate_gff(trans_gff3)

    orfs = create_transcriptome_orf_db(str(trans_gff3), str(config.genome))
    featuredb2gff3_file(orfs, sdir / output)
    return sdir / output


# ---------------------------------------------------------------------------
# Step helpers
# ---------------------------------------------------------------------------


def _run_step(
    config: PipelineConfig,
    manifest: RunManifest,
    step_name: str,
    fn: Callable,
    *args,
    **kwargs,
) -> Path:
    """Run a single pipeline step with manifest tracking via context manager.

    Skips if already completed in a previous run.
    """
    cached = is_step_complete(manifest, step_name)
    if cached:
        return cached

    # Clean up interrupted previous runs before starting
    if is_step_interrupted(config.work_dir, step_name):
        log.warning("[%s] Cleaning up interrupted previous run...", step_name)
        clean_interrupted_step(config.work_dir, step_name)

    with pipeline_step(config.work_dir, manifest, step_name) as step:
        result_path = fn(config, *args, **kwargs)
        step.output_file = str(result_path)
        return result_path


def _log_prediction_count(label: str, gff3_path: Path) -> None:
    """Log the number of gene predictions in a GFF3 file."""
    n = count_gff3_features(gff3_path)
    log.info("%s: %d gene predictions", label, n)


def _run_concurrent_steps(
    config: PipelineConfig,
    manifest: RunManifest,
    tasks: list[tuple[str, Callable, tuple, dict]],
) -> dict[str, Path]:
    """Run multiple independent steps concurrently.

    Args:
        tasks: List of (step_name, fn, args, kwargs) tuples.

    Returns:
        Dict of step_name -> output Path.
    """
    results: dict[str, Path] = {}

    with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
        future_to_name = {
            pool.submit(_run_step, config, manifest, name, fn, *args, **kwargs): name
            for name, fn, args, kwargs in tasks
        }
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            results[name] = future.result()

    return results


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


def run_annotation_pipeline(config: PipelineConfig) -> Path:
    """Run the full annotation pipeline.

    Features:
    - Writes eukan-run.json manifest for tracking/reproducibility
    - Skips completed steps on resume
    - Cleans up interrupted steps automatically
    - Runs independent steps concurrently where possible

    Returns:
        Path to the final GFF3 output.
    """
    validate_fasta(config.genome)

    # Sanitize genome headers (strip descriptions that break GFF tools)
    sanitized_genome = sanitize_genome_fasta(config.genome, config.work_dir)
    if sanitized_genome != config.genome:
        # Update config with sanitized genome path
        config = PipelineConfig(
            **{**config.model_dump(), "genome": sanitized_genome}
        )

    if config.has_transcripts:
        log.info("Transcript evidence: %s, %s, %s",
                 config.transcripts_fasta.name, config.transcripts_gff.name, config.rnaseq_hints.name)
    else:
        log.info("No transcript evidence found. Running protein-only annotation.")

    # Load or create manifest
    manifest = load_manifest(config.work_dir)
    if manifest and manifest.status == "completed":
        log.info("Previous run completed. Delete eukan-run.json to re-run.")
        final = config.work_dir / "final.gff3"
        if final.exists():
            return final
    manifest = init_manifest(config)
    save_manifest(config.work_dir, manifest)

    try:
        result = _execute_steps(config, manifest)
        manifest.status = "completed"
        manifest.finished_at = manifest.steps[
            max(manifest.steps, key=lambda k: manifest.steps[k].finished_at or "")
        ].finished_at
        save_manifest(config.work_dir, manifest)
        return result
    except Exception:
        manifest.status = "failed"
        save_manifest(config.work_dir, manifest)
        raise


def _execute_steps(config: PipelineConfig, manifest: RunManifest) -> Path:
    """Execute all pipeline steps with manifest tracking and concurrency."""
    ev: dict[str, Path] = {}

    if config.has_transcripts:
        # ORF finding and GeneMark are independent -- run concurrently
        concurrent = _run_concurrent_steps(config, manifest, [
            ("orf_finder", find_orfs, (config.transcripts_gff,), {}),
            ("genemark", run_genemark, (config.rnaseq_hints,), {}),
        ])
        ev["transcriptORFs"] = concurrent["orf_finder"]
        ev["genemark"] = concurrent["genemark"]
        _log_prediction_count("GeneMark", ev["genemark"])

        intron_hints = config.work_dir / "genemark" / "introns.gff"
        ev["spaln"] = _run_step(
            config, manifest, "prot_align", align_proteins,
            ev["genemark"], config.proteins,
            intron_hints if intron_hints.exists() else None,
        )
        ev["augustus"] = _run_step(
            config, manifest, "augustus", run_augustus,
            ev["genemark"], ev["spaln"], ev["transcriptORFs"],
        )
        _log_prediction_count("AUGUSTUS", ev["augustus"])

        if config.is_fungus:
            # SNAP and CodingQuarry are independent -- run concurrently
            concurrent = _run_concurrent_steps(config, manifest, [
                ("snap", run_snap, (ev["augustus"], ev["spaln"], ev["transcriptORFs"]), {}),
                ("codingquarry", run_codingquarry, (config.transcripts_gff,), {}),
            ])
            ev["snap"] = concurrent["snap"]
            ev["codingquarry"] = concurrent["codingquarry"]
            _log_prediction_count("SNAP", ev["snap"])
            _log_prediction_count("CodingQuarry", ev["codingquarry"])
            return _run_step(
                config, manifest, "evm_consensus_models", build_consensus_models,
                ev["spaln"], ev["augustus"], ev["snap"],
                ev["codingquarry"], config.transcripts_gff,
            )
        else:
            return _run_step(
                config, manifest, "evm_consensus_models", build_consensus_models,
                ev["spaln"], ev["augustus"], config.transcripts_gff,
            )
    else:
        ev["genemark"] = _run_step(config, manifest, "genemark", run_genemark)
        _log_prediction_count("GeneMark", ev["genemark"])
        ev["spaln"] = _run_step(
            config, manifest, "prot_align", align_proteins,
            ev["genemark"], config.proteins,
        )
        ev["augustus"] = _run_step(
            config, manifest, "augustus", run_augustus,
            ev["genemark"], ev["spaln"],
        )
        _log_prediction_count("AUGUSTUS", ev["augustus"])

        if config.is_fungus or config.is_protist:
            concurrent = _run_concurrent_steps(config, manifest, [
                ("snap", run_snap, (ev["augustus"], ev["spaln"]), {}),
                ("codingquarry", run_codingquarry, (ev["augustus"],), {}),
            ])
            ev["snap"] = concurrent["snap"]
            ev["codingquarry"] = concurrent["codingquarry"]
            _log_prediction_count("SNAP", ev["snap"])
            _log_prediction_count("CodingQuarry", ev["codingquarry"])
            return _run_step(
                config, manifest, "evm_consensus_models", build_consensus_models,
                ev["spaln"], ev["augustus"], ev["snap"], ev["codingquarry"],
            )
        else:
            ev["snap"] = _run_step(
                config, manifest, "snap", run_snap,
                ev["augustus"], ev["spaln"],
            )
            _log_prediction_count("SNAP", ev["snap"])
            return _run_step(
                config, manifest, "evm_consensus_models", build_consensus_models,
                ev["spaln"], ev["augustus"], ev["snap"], ev["genemark"],
            )
