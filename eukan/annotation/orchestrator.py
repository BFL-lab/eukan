"""Annotation pipeline orchestration: step ordering, concurrency, and manifest tracking."""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from eukan.annotation.alignment import align_proteins
from eukan.annotation.augustus import run_augustus
from eukan.annotation.consensus import build_consensus_models
from eukan.annotation.genemark import run_genemark
from eukan.annotation.orf import create_transcriptome_orf_db
from eukan.annotation.snap import run_codingquarry, run_snap
from eukan.annotation.validation import sanitize_genome_fasta, validate_fasta
from eukan.gff.io import featuredb2gff3_file
from eukan.infra.logging import count_gff3_features, get_logger
from eukan.infra.manifest import (
    ANNOTATION,
    RunManifest,
    get_or_create_manifest,
    run_orchestrated_step,
    save_manifest,
    step_key,
    validate_step_outputs,
)
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# ORF finding step
# ---------------------------------------------------------------------------


def find_orfs(config: PipelineConfig, trans_gff3: Path) -> Path:
    """Find ORFs in transcript assemblies."""
    from eukan.infra.logging import validate_gff

    output = "transcript_orfs.gff3"
    sdir = step_dir(config.work_dir, "orf_finder")
    log.info("Finding ORFs in transcript assemblies...")
    validate_gff(trans_gff3)

    orfs = create_transcriptome_orf_db(str(trans_gff3), str(config.genome), genetic_code=int(config.genetic_code))
    featuredb2gff3_file(orfs, sdir / output)
    return sdir / output


# ---------------------------------------------------------------------------
# Step helpers
# ---------------------------------------------------------------------------


def _run_step(
    config: PipelineConfig,
    manifest: RunManifest,
    name: str,
    fn: Callable,
    *args,
    **kwargs,
) -> Path:
    """Run an annotation step via the shared orchestrator helper.

    All annotation steps produce a GFF3 output, so the result is always
    a Path (never None).
    """
    result = run_orchestrated_step(
        config.manifest_dir, manifest, step_key(ANNOTATION, name),
        fn, config, *args,
        step_dir=config.work_dir / name,
        **kwargs,
    )
    assert result is not None, f"annotation step {name!r} returned no output"
    return result


def _log_prediction_count(label: str, gff3_path: Path) -> None:
    """Log the number of gene predictions in a GFF3 file."""
    log.info("%s: %d gene predictions", label, count_gff3_features(gff3_path))


def _run_concurrent_steps(
    config: PipelineConfig,
    manifest: RunManifest,
    tasks: list[tuple[str, Callable, tuple, dict]],
) -> dict[str, Path]:
    """Run multiple independent steps concurrently."""
    results: dict[str, Path] = {}
    with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
        future_to_name = {
            pool.submit(_run_step, config, manifest, name, fn, *args, **kwargs): name
            for name, fn, args, kwargs in tasks
        }
        for future in as_completed(future_to_name):
            results[future_to_name[future]] = future.result()
    return results


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


def run_annotation_pipeline(
    config: PipelineConfig,
    force_steps: list[str] | None = None,
) -> Path:
    """Run the full annotation pipeline.

    Features:
    - Writes eukan-run.json manifest for tracking/reproducibility
    - Skips completed steps on resume
    - Cleans up interrupted steps automatically
    - Runs independent steps concurrently where possible

    Args:
        config: Pipeline configuration.
        force_steps: Optional list of step names to force re-run, even if
            previously completed.  When provided, step records are removed
            from the manifest so they will be re-executed.

    Returns:
        Path to the final GFF3 output.
    """
    validate_fasta(config.genome)

    # Sanitize genome headers (strip descriptions that break GFF tools)
    sanitized_genome = sanitize_genome_fasta(config.genome, config.work_dir)
    if sanitized_genome != config.genome:
        config = config.model_copy(update={"genome": sanitized_genome})

    if config.has_transcripts:
        # has_transcripts is True iff all three are non-None
        assert (
            config.transcripts_fasta is not None
            and config.transcripts_gff is not None
            and config.rnaseq_hints is not None
        )
        log.info("Transcript evidence: %s, %s, %s",
                 config.transcripts_fasta.name, config.transcripts_gff.name, config.rnaseq_hints.name)
    else:
        log.info("No transcript evidence found. Running protein-only annotation.")

    # Load or create shared manifest (used by all pipelines in this work_dir)
    manifest = get_or_create_manifest(config.manifest_dir, config)

    if force_steps:
        # Explicit re-run: remove only the requested steps
        for step in force_steps:
            manifest.steps.pop(step, None)
        save_manifest(config.manifest_dir, manifest)
    else:
        # Validate manifest: check completed steps have valid output
        expected = _expected_steps(config)
        errors = _validate_step_outputs(manifest, expected)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

        # Check if there's any work to do
        pending = [s for s in expected if s not in manifest.steps]
        if not pending:
            final = config.work_dir / "final.gff3"
            if final.exists():
                log.info("All steps complete. Use --run-* flags to re-run specific steps.")
                return final
        save_manifest(config.manifest_dir, manifest)

    try:
        result = _execute_steps(config, manifest)
        manifest.status = "completed"
        manifest.finished_at = manifest.steps[
            max(manifest.steps, key=lambda k: manifest.steps[k].finished_at or "")
        ].finished_at
        save_manifest(config.manifest_dir, manifest)
        return result
    except Exception:
        manifest.status = "failed"
        save_manifest(config.manifest_dir, manifest)
        raise


# Map from manifest key to the --run-* CLI flag that forces re-run
_STEP_TO_FLAG = {
    step_key(ANNOTATION, name): flag
    for name, flag in {
        "genemark": "--run-genemark",
        "prot_align": "--run-prot-align",
        "prot_align_ssp": "--run-prot-align",
        "augustus": "--run-augustus",
        "snap": "--run-snap",
        "codingquarry": "--run-snap",
        "orf_finder": "--run-genemark",
        "evm_consensus_models": "--run-consensus",
    }.items()
}


def _validate_step_outputs(manifest: RunManifest, expected: list[str]) -> list[str]:
    """Validate that completed steps have non-empty output files."""
    return validate_step_outputs(manifest, expected, _STEP_TO_FLAG)


def _expected_steps(config: PipelineConfig) -> list[str]:
    """Return the list of manifest step keys expected for this config."""
    prot_align_step = "prot_align_ssp" if config.spaln_ssp else "prot_align"
    steps = ["genemark", prot_align_step, "augustus"]
    if config.has_transcripts:
        steps.insert(0, "orf_finder")
    if config.is_fungus or config.is_protist:
        steps.extend(["snap", "codingquarry"])
    else:
        steps.append("snap")
    steps.append("evm_consensus_models")
    return [step_key(ANNOTATION, s) for s in steps]


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
        prot_align_step = "prot_align_ssp" if config.spaln_ssp else "prot_align"
        ev["spaln"] = _run_step(
            config, manifest, prot_align_step, align_proteins,
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
        prot_align_step = "prot_align_ssp" if config.spaln_ssp else "prot_align"
        ev["spaln"] = _run_step(
            config, manifest, prot_align_step, align_proteins,
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
