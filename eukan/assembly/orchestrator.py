"""Assembly pipeline orchestration."""

from __future__ import annotations

from pathlib import Path

from eukan.assembly.pasa import run_pasa
from eukan.assembly.star import map_reads
from eukan.assembly.trinity import run_trinity
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    RunManifest, get_or_create_manifest, is_step_complete,
    is_step_interrupted, clean_interrupted_step, pipeline_step, save_manifest,
    validate_step_outputs,
)
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# Step names use a pipeline prefix for the shared manifest
_STEP_MAP = {
    "map": "assembly/star",
    "trinity": "assembly/trinity",
    "pasa": "assembly/pasa",
}

# Map from manifest key to CLI flag
_STEP_TO_FLAG = {
    "assembly/star": "-A / --run-star",
    "assembly/trinity": "-T / --run-trinity",
    "assembly/pasa": "-P / --run-pasa",
}

# Key output file for each step
_STEP_OUTPUTS = {
    "map": "STAR_Aligned.sortedByCoord.out.bam",
    "trinity": "trinity-gg.fasta",
    "pasa": "nr_transcripts.fasta",
}


def run_assembly(config: AssemblyConfig, steps: list[str], force: bool = False) -> None:
    """Run the specified assembly steps with manifest tracking."""
    manifest = get_or_create_manifest(config.manifest_dir)

    if not force:
        # Validate completed steps have valid output
        expected = [_STEP_MAP[s] for s in steps]
        errors = _validate_step_outputs(manifest, expected, config.work_dir)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

    save_manifest(config.manifest_dir, manifest)

    for step_key in ("map", "trinity", "pasa"):
        if step_key not in steps:
            continue
        step_name = _STEP_MAP[step_key]
        _run_assembly_step(config, manifest, step_name, step_key, force)


def _validate_step_outputs(
    manifest: RunManifest, expected: list[str], work_dir: Path,
) -> list[str]:
    """Validate that completed steps have non-empty output files."""
    return validate_step_outputs(manifest, expected, _STEP_TO_FLAG)


def _run_assembly_step(
    config: AssemblyConfig,
    manifest: RunManifest,
    step_name: str,
    step_key: str,
    force: bool,
) -> None:
    """Run a single assembly step with manifest tracking."""
    if not force:
        cached = is_step_complete(manifest, step_name)
        if cached:
            return

    if is_step_interrupted(config.manifest_dir, step_name, step_dir=config.work_dir / step_name):
        log.warning("[%s] Cleaning up interrupted previous run...", step_name)
        clean_interrupted_step(config.manifest_dir, step_name, step_dir=config.work_dir / step_name)

    with pipeline_step(config.manifest_dir, manifest, step_name, step_dir=config.work_dir / step_name) as step:
        step_fns = {"map": map_reads, "trinity": run_trinity, "pasa": run_pasa}
        step_fns[step_key](config, force=force)
        output_path = config.work_dir / _STEP_OUTPUTS[step_key]
        if output_path.exists():
            step.output_file = str(output_path)
