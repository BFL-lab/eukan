"""Repeat-masking pipeline orchestration."""

from __future__ import annotations

from pathlib import Path

from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    REPEATS,
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.pipeline import run_orchestrated_step
from eukan.infra.steps import is_step_complete, validate_or_raise
from eukan.repeats.masker import run_masker
from eukan.repeats.modeler import run_modeler
from eukan.settings import RepeatsConfig

log = get_logger(__name__)


_MODELER = "modeler"
_MASKER = "masker"
_ALL_NAMES = (_MODELER, _MASKER)

_STEP_TO_FLAG = {
    step_key(REPEATS, _MODELER): "--run-modeler",
    step_key(REPEATS, _MASKER): "--run-masker",
}
_ALL_STEP_KEYS = [step_key(REPEATS, n) for n in _ALL_NAMES]


def force_steps_from_run_flags(
    *,
    run_modeler: bool = False,
    run_masker: bool = False,
    force: bool = False,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force.

    Returns full ``repeats/<step>`` keys, harmonized with the annotation
    and assembly pipelines' ``force_steps_from_run_flags``.
    """
    selected: list[str] = []
    if run_modeler:
        selected.append(_MODELER)
    if run_masker:
        selected.append(_MASKER)
    if selected:
        return [step_key(REPEATS, s) for s in selected]
    if force:
        return list(_ALL_STEP_KEYS)
    return []


def _modeler_output(config: RepeatsConfig) -> Path:
    """Path to the families library produced by the modeler step."""
    return config.work_dir / _MODELER / f"{config.name}.replib-families.fa"


def run_repeats(
    config: RepeatsConfig,
    *,
    force_steps: list[str] | None = None,
) -> Path:
    """Run the repeat-masking pipeline.

    Returns the path to the softmasked genome.

    Args:
        force_steps: Manifest keys to re-run. ``None`` or empty means
            "run all pending steps; skip cached". A non-empty list
            narrows the active step set to just those keys *and* forces
            re-execution.
    """
    manifest = get_or_create_manifest(config.manifest_dir, config)
    forced = set(force_steps or ())

    if forced:
        active = [name for name in _ALL_NAMES if step_key(REPEATS, name) in forced]
    else:
        active = list(_ALL_NAMES)
        expected = [step_key(REPEATS, s) for s in active if s != _MODELER or not config.lib]
        validate_or_raise(manifest, expected, _STEP_TO_FLAG)

    save_manifest(config.manifest_dir, manifest)

    families: Path | None = None
    if _MODELER in active:
        if config.lib:
            log.info("Using user-provided repeat library: %s — skipping RepeatModeler.", config.lib)
            families = config.lib
        else:
            modeler_step = step_key(REPEATS, _MODELER)
            cached = is_step_complete(manifest, modeler_step) if modeler_step not in forced else None
            if cached:
                families = cached
            else:
                run_orchestrated_step(
                    config.manifest_dir, manifest, modeler_step,
                    _run_modeler_step, config,
                    step_dir=config.work_dir / _MODELER,
                    force=modeler_step in forced,
                    output_file=_modeler_output(config),
                )
                families = _modeler_output(config)

    masked: Path | None = None
    if _MASKER in active:
        if families is None:
            # Masker requested without modeler — pull the library from
            # the manifest record or from --lib.
            if config.lib:
                families = config.lib
            else:
                cached = is_step_complete(manifest, step_key(REPEATS, _MODELER))
                if cached is None:
                    from eukan.exceptions import ConfigurationError
                    raise ConfigurationError(
                        "RepeatMasker step requested but no families library found.",
                        hint="Run --run-modeler first or pass --lib.",
                    )
                families = cached

        masker_step = step_key(REPEATS, _MASKER)
        masked_output = config.work_dir / f"{config.name}.masked.fasta"
        run_orchestrated_step(
            config.manifest_dir, manifest, masker_step,
            _run_masker_step, config, families,
            step_dir=config.work_dir / _MASKER,
            force=masker_step in forced,
            output_file=masked_output,
        )
        masked = masked_output

    return masked or (config.work_dir / f"{config.name}.masked.fasta")


def _run_modeler_step(config: RepeatsConfig) -> Path:
    return run_modeler(config)


def _run_masker_step(config: RepeatsConfig, families: Path) -> Path:
    masked, _hints = run_masker(config, families)
    return masked
