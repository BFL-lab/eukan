"""eukan mask-repeats — RepeatModeler + RepeatMasker softmasking pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PreformattedEpilogCommand,
    drop_none,
    force_option,
    genome_option,
    numcpu_option,
)


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
    "--lib", type=click.Path(exists=True, path_type=Path), default=None,
    help="Pre-built repeat-family library FASTA. When set, RepeatModeler is skipped.",
)
@optgroup.group("Re-run steps")
@optgroup.option("--run-modeler", is_flag=True, help="Force re-run BuildDatabase + RepeatModeler.")
@optgroup.option("--run-masker", is_flag=True, help="Force re-run RepeatMasker.")
@force_option
def mask_repeats(
    genome: Path,
    numcpu: int,
    engine: str,
    lib: Path | None,
    run_modeler: bool,
    run_masker: bool,
    force: bool,
) -> None:
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
    from eukan.repeats import run_repeats
    from eukan.repeats.pipeline import force_steps_from_run_flags
    from eukan.settings import RepeatsConfig

    config = RepeatsConfig(**drop_none(
        genome=genome.resolve(),
        work_dir=Path.cwd(),
        num_cpu=numcpu,
        engine=engine.lower(),
        lib=lib.resolve() if lib else None,
    ))

    force_steps = force_steps_from_run_flags(
        run_modeler=run_modeler, run_masker=run_masker, force=force,
    )
    masked = run_repeats(config, force_steps=force_steps or None)
    click.echo(f"Done. Masked genome: {masked}")
    click.echo(
        f"Pass it to the next stage with: eukan annotate -g {masked.name} ..."
    )
