"""eukan mask-repeats — RepeatModeler + RepeatMasker softmasking pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PreformattedEpilogCommand,
    force_option,
    genome_option,
    given,
    numcpu_option,
)
from eukan.infra.layout import step_work_dir


@click.command("mask-repeats", cls=PreformattedEpilogCommand)
@optgroup.group("Required input")
@genome_option("Genome sequence in FASTA format.")
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--engine", type=click.Choice(["rmblast", "ncbi"], case_sensitive=False),
    default="rmblast", show_default=True,
    help="Search engine for BuildDatabase / RepeatMasker.",
)
@optgroup.option(
    "--lib", type=click.Path(exists=True, resolve_path=True, path_type=Path), default=None,
    help="Pre-built repeat-family library FASTA. When set, RepeatModeler is skipped.",
)
@optgroup.group("Re-run steps")
@optgroup.option(
    "--run-modeler", is_flag=True,
    help="Force re-run BuildDatabase + RepeatModeler (also re-runs RepeatMasker).",
)
@optgroup.option("--run-masker", is_flag=True, help="Force re-run RepeatMasker.")
@force_option
def mask_repeats(**opts: Any) -> None:
    """Soft-mask repeats with RepeatModeler + RepeatMasker.

    \b
    Produces, in the working directory:
      <stem>.masked.fasta     soft-masked genome (lower-case in repeats)
      <stem>.repeats.gff      raw RepeatMasker GFF
      hints_repeatmask.gff    AUGUSTUS-format hints (auto-discovered by
                              `eukan annotate`)

    \b
    Pass the masked genome to `eukan annotate -g <stem>.masked.fasta`.
    """
    from eukan.repeats import run_repeats  # lazy: startup cost
    from eukan.repeats.pipeline import force_steps_from_run_flags  # lazy: startup cost
    from eukan.settings import RepeatsConfig  # lazy: startup cost

    run_flags = {k: opts.pop(k) for k in list(opts) if k.startswith("run_") or k == "force"}
    config = RepeatsConfig(
        **given(opts), work_dir=step_work_dir("mask-repeats"), manifest_dir=Path.cwd(),
    )

    force_steps = force_steps_from_run_flags(**run_flags)
    masked = run_repeats(config, force_steps=force_steps or None)
    click.echo(f"Done. Masked genome: {masked}")
    click.echo(
        f"Pass it to the next stage with: eukan annotate -g {masked.name} ..."
    )
