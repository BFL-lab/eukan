"""Assembly pipeline orchestration."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from eukan.assembly.pasa import run_pasa
from eukan.assembly.star import map_reads
from eukan.assembly.trinity import run_trinity
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    ASSEMBLY,
    get_or_create_manifest,
    run_orchestrated_step,
    save_manifest,
    step_key,
    validate_step_outputs,
)
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


@dataclass(frozen=True)
class _AssemblyStep:
    fn: Callable[..., None]
    output: str          # filename under work_dir
    flag: str            # re-run CLI flag shown to the user


# Keyed by the canonical step name; CLI accepts legacy aliases via _ALIASES.
_STEPS: dict[str, _AssemblyStep] = {
    "star":    _AssemblyStep(map_reads,   "STAR_Aligned.sortedByCoord.out.bam", "-A / --run-star"),
    "trinity": _AssemblyStep(run_trinity, "trinity-gg.fasta",                    "-T / --run-trinity"),
    "pasa":    _AssemblyStep(run_pasa,    "nr_transcripts.fasta",                "-P / --run-pasa"),
}

_ALIASES = {"map": "star"}

_STEP_TO_FLAG = {step_key(ASSEMBLY, n): s.flag for n, s in _STEPS.items()}


def run_assembly(config: AssemblyConfig, steps: list[str], force: bool = False) -> None:
    """Run the specified assembly steps with manifest tracking."""
    manifest = get_or_create_manifest(config.manifest_dir)
    requested = [_ALIASES.get(s, s) for s in steps]

    if not force:
        expected = [step_key(ASSEMBLY, s) for s in requested]
        errors = validate_step_outputs(manifest, expected, _STEP_TO_FLAG)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

    save_manifest(config.manifest_dir, manifest)

    for name, spec in _STEPS.items():
        if name not in requested:
            continue
        run_orchestrated_step(
            config.manifest_dir, manifest, step_key(ASSEMBLY, name),
            spec.fn, config,
            step_dir=config.work_dir / name,
            force=force,
            output_file=config.work_dir / spec.output,
        )
