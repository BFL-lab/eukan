"""Assembly pipeline orchestration."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from eukan.assembly.pasa import run_pasa
from eukan.assembly.star import map_reads
from eukan.assembly.trinity import run_trinity
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    ASSEMBLY,
    get_or_create_manifest,
    run_orchestrated_step,
    save_manifest,
    step_key,
    validate_or_raise,
)
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


@dataclass(frozen=True)
class _AssemblyStep:
    fn: Callable[..., None]
    output: str          # filename under work_dir
    flag: str            # re-run CLI flag shown to the user


_STEPS: dict[str, _AssemblyStep] = {
    "star":    _AssemblyStep(map_reads,   "STAR_Aligned.sortedByCoord.out.bam",  "-A / --run-star"),
    "trinity": _AssemblyStep(run_trinity, "trinity-gg.fasta",                    "-T / --run-trinity"),
    "pasa":    _AssemblyStep(run_pasa,    Artifact.NR_TRANSCRIPTS_FASTA.value,   "-P / --run-pasa"),
}

_STEP_TO_FLAG = {step_key(ASSEMBLY, n): s.flag for n, s in _STEPS.items()}
_ALL_STEP_KEYS = [step_key(ASSEMBLY, n) for n in _STEPS]


def force_steps_from_run_flags(
    *,
    run_star: bool = False,
    run_trinity: bool = False,
    run_pasa: bool = False,
    force: bool = False,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force.

    Returns full ``assembly/<step>`` keys, harmonized with the annotation
    pipeline's ``force_steps_from_run_flags``. Semantics:

    * ``[]`` — no flags set: pipeline runs all pending steps; cached steps
      are skipped.
    * ``[some keys]`` — ``--run-X`` flags set: pipeline runs only those
      steps and forces re-execution.
    * ``[all keys]`` — ``--force`` set with no ``--run-X``: pipeline
      re-runs every step from scratch.
    """
    flag_to_step = {
        "run_star": "star",
        "run_trinity": "trinity",
        "run_pasa": "pasa",
    }
    flag_states = {"run_star": run_star, "run_trinity": run_trinity, "run_pasa": run_pasa}
    selected = [step for flag, step in flag_to_step.items() if flag_states[flag]]
    if selected:
        return [step_key(ASSEMBLY, s) for s in selected]
    if force:
        return list(_ALL_STEP_KEYS)
    return []


def run_assembly(
    config: AssemblyConfig,
    *,
    force_steps: list[str] | None = None,
) -> None:
    """Run the assembly pipeline with manifest tracking.

    Args:
        force_steps: Manifest keys to re-run from scratch. ``None`` or
            empty means "run all pending steps; skip cached". A non-empty
            list narrows the active step set to just those keys *and*
            forces re-execution.
    """
    manifest = get_or_create_manifest(config.manifest_dir, config)
    forced = set(force_steps or ())

    if forced:
        active = [name for name in _STEPS if step_key(ASSEMBLY, name) in forced]
    else:
        active = list(_STEPS)
        # Validate cached outputs only when not re-running everything.
        expected = [step_key(ASSEMBLY, s) for s in active]
        validate_or_raise(manifest, expected, _STEP_TO_FLAG)

    save_manifest(config.manifest_dir, manifest)

    for name in active:
        spec = _STEPS[name]
        run_orchestrated_step(
            config.manifest_dir, manifest, step_key(ASSEMBLY, name),
            spec.fn, config,
            step_dir=config.work_dir / name,
            force=step_key(ASSEMBLY, name) in forced,
            output_file=config.work_dir / spec.output,
        )
