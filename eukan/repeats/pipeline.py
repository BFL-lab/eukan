"""Repeat-masking pipeline: RepeatModeler library → RepeatMasker softmasking."""

from __future__ import annotations

from pathlib import Path

from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    REPEATS,
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.pipeline import (
    StepSpec,
    force_steps_from_run_flags as _force_steps_from_run_flags,
    run_orchestrated_step,
)
from eukan.infra.steps import is_step_complete, validate_or_raise
from eukan.repeats.masker import run_masker
from eukan.repeats.modeler import run_modeler
from eukan.settings import RepeatsConfig

log = get_logger(__name__)


def _run_modeler_step(config: RepeatsConfig) -> Path:
    return run_modeler(config)


def _run_masker_step(config: RepeatsConfig, families: Path) -> Path:
    masked, _hints = run_masker(config, families)
    return masked


# The masker step is declared without a concrete fn since it's invoked
# with an extra positional argument (the families path) that the simple
# driver doesn't supply. The custom run loop below threads that arg.
_STEPS: list[StepSpec] = [
    StepSpec("modeler", _run_modeler_step, flag="--run-modeler"),
    StepSpec("masker",  _run_masker_step,  flag="--run-masker"),
]


def force_steps_from_run_flags(
    *,
    run_modeler: bool = False,
    run_masker: bool = False,
    force: bool = False,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force."""
    return _force_steps_from_run_flags(
        REPEATS, _STEPS,
        force=force,
        run_modeler=run_modeler, run_masker=run_masker,
    )


def _modeler_output(config: RepeatsConfig) -> Path:
    """Path to the families library produced by the modeler step."""
    return config.work_dir / "modeler" / f"{config.name}.replib-families.fa"


def run_repeats(
    config: RepeatsConfig,
    *,
    force_steps: list[str] | None = None,
) -> Path:
    """Run the repeat-masking pipeline.

    The simple driver doesn't fit because the modeler step is conditional
    on ``config.lib`` and the masker takes the modeler's output as an
    extra argument. This loop is a small custom driver around the same
    ``run_orchestrated_step`` primitive.
    """
    manifest = get_or_create_manifest(config.manifest_dir, config)
    forced = set(force_steps or ())
    flag_map = {step_key(REPEATS, s.name): s.flag for s in _STEPS if s.flag}

    if forced:
        active = [s.name for s in _STEPS if step_key(REPEATS, s.name) in forced]
    else:
        active = [s.name for s in _STEPS]
        expected = [step_key(REPEATS, n) for n in active if n != "modeler" or not config.lib]
        validate_or_raise(manifest, expected, flag_map)

    save_manifest(config.manifest_dir, manifest)

    families: Path | None = None
    if "modeler" in active:
        if config.lib:
            log.info("Using user-provided repeat library: %s — skipping RepeatModeler.", config.lib)
            families = config.lib
        else:
            modeler_key = step_key(REPEATS, "modeler")
            cached = is_step_complete(manifest, modeler_key) if modeler_key not in forced else None
            if cached:
                families = cached
            else:
                run_orchestrated_step(
                    config.manifest_dir, manifest, modeler_key,
                    _run_modeler_step, config,
                    step_dir=config.work_dir / "modeler",
                    force=modeler_key in forced,
                    output_file=_modeler_output(config),
                )
                families = _modeler_output(config)

    masked_output = config.work_dir / f"{config.name}.masked.fasta"
    if "masker" in active:
        if families is None:
            if config.lib:
                families = config.lib
            else:
                cached = is_step_complete(manifest, step_key(REPEATS, "modeler"))
                if cached is None:
                    from eukan.exceptions import ConfigurationError
                    raise ConfigurationError(
                        "RepeatMasker step requested but no families library found.",
                        hint="Run --run-modeler first or pass --lib.",
                    )
                families = cached

        masker_key = step_key(REPEATS, "masker")
        run_orchestrated_step(
            config.manifest_dir, manifest, masker_key,
            _run_masker_step, config, families,
            step_dir=config.work_dir / "masker",
            force=masker_key in forced,
            output_file=masked_output,
        )

    return masked_output
