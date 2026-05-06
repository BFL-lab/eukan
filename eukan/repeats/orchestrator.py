"""Repeat-masking pipeline orchestration."""

from __future__ import annotations

from pathlib import Path

from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    REPEATS,
    get_or_create_manifest,
    is_step_complete,
    run_orchestrated_step,
    save_manifest,
    step_key,
    validate_step_outputs,
)
from eukan.repeats.masker import run_masker
from eukan.repeats.modeler import run_modeler
from eukan.settings import RepeatsConfig

log = get_logger(__name__)


_MODELER = "modeler"
_MASKER = "masker"

_STEP_TO_FLAG = {
    step_key(REPEATS, _MODELER): "--run-modeler",
    step_key(REPEATS, _MASKER): "--run-masker",
}


def steps_and_force_from_run_flags(
    *,
    run_modeler: bool = False,
    run_masker: bool = False,
    force: bool = False,
) -> tuple[list[str], bool]:
    """Translate per-flag CLI booleans into ``(steps, force)`` for ``run_repeats``.

    Mirrors ``eukan.assembly.orchestrator.steps_and_force_from_run_flags``.
    """
    selected: list[str] = []
    if run_modeler:
        selected.append(_MODELER)
    if run_masker:
        selected.append(_MASKER)
    if selected:
        return selected, True
    return [_MODELER, _MASKER], force


def _modeler_output(config: RepeatsConfig) -> Path:
    """Path to the families library produced by the modeler step."""
    return config.work_dir / "repeats" / "modeler" / f"{config.name}.replib-families.fa"


def run_repeats(config: RepeatsConfig, steps: list[str], force: bool = False) -> Path:
    """Run repeat masking. Returns the path to the softmasked genome."""
    manifest = get_or_create_manifest(config.manifest_dir)

    if not force:
        expected = [step_key(REPEATS, s) for s in steps if s != _MODELER or not config.lib]
        errors = validate_step_outputs(manifest, expected, _STEP_TO_FLAG)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

    save_manifest(config.manifest_dir, manifest)

    families: Path | None = None
    if _MODELER in steps:
        if config.lib:
            log.info("Using user-provided repeat library: %s — skipping RepeatModeler.", config.lib)
            families = config.lib
        else:
            modeler_step = step_key(REPEATS, _MODELER)
            cached = is_step_complete(manifest, modeler_step) if not force else None
            if cached:
                families = cached
            else:
                run_orchestrated_step(
                    config.manifest_dir, manifest, modeler_step,
                    _run_modeler_step, config,
                    step_dir=config.work_dir / "repeats" / _MODELER,
                    force=force,
                    output_file=_modeler_output(config),
                )
                families = _modeler_output(config)

    masked: Path | None = None
    if _MASKER in steps:
        if families is None:
            # Masker requested without modeler — pull the library from the
            # manifest record or from --lib.
            if config.lib:
                families = config.lib
            else:
                cached = is_step_complete(manifest, step_key(REPEATS, _MODELER))
                if cached is None:
                    log.error(
                        "RepeatMasker step requested but no families library found. "
                        "Run --run-modeler first or pass --lib."
                    )
                    raise SystemExit(1)
                families = cached

        masker_step = step_key(REPEATS, _MASKER)
        masked_output = config.work_dir / f"{config.name}.masked.fasta"
        run_orchestrated_step(
            config.manifest_dir, manifest, masker_step,
            _run_masker_step, config, families,
            step_dir=config.work_dir / "repeats" / _MASKER,
            force=force,
            output_file=masked_output,
        )
        masked = masked_output

    return masked or (config.work_dir / f"{config.name}.masked.fasta")


def _run_modeler_step(config: RepeatsConfig) -> Path:
    return run_modeler(config)


def _run_masker_step(config: RepeatsConfig, families: Path) -> Path:
    masked, _hints = run_masker(config, families)
    return masked
