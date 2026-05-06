"""eukan assemble — transcriptome assembly pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PASA_CODE_TABLE,
    PreformattedEpilogCommand,
    drop_none,
    force_option,
    genome_option,
    numcpu_option,
)


@click.command(cls=PreformattedEpilogCommand, epilog=PASA_CODE_TABLE)
@optgroup.group("Required input")
@genome_option("Genome FASTA file.")
@optgroup.option("--left", "-l", type=click.Path(exists=True, path_type=Path), help="Left paired-end reads.")
@optgroup.option("--right", "-r", type=click.Path(exists=True, path_type=Path), help="Right paired-end reads.")
@optgroup.option("--single", "-s", type=click.Path(exists=True, path_type=Path), help="Single-end reads.")
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--strand-specific", "-S", type=click.Choice(["RF", "FR", "R", "F"]), default=None,
    help="Strand-specific library type.",
)
@optgroup.option("--align-mode", "-t", type=click.Choice(["EndToEnd", "Local"]), default="Local", show_default=True)
@optgroup.option(
    "--splice-permissive", is_flag=True, default=False,
    help="Allow non-canonical splice sites (GC-AG, AT-AC). "
    "Sets PASA splice boundary stringency to 0 and retains non-canonical junctions.",
)
@optgroup.option(
    "--genetic-code", "-c",
    type=click.Choice(["1", "6", "10", "12"]),
    default="1", show_default=True,
    help="Genetic code for PASA. Supported: 1=standard, 6=Tetrahymena, 10=Euplotes, 12=Candida.",
)
@optgroup.option("--min-intron", "-m", type=int, default=20, show_default=True, help="Minimum intron length.")
@optgroup.option("--max-intron", "-M", type=int, default=5000, show_default=True, help="Maximum intron length.")
@optgroup.option("--phred", type=click.Choice(["33", "64"]), default="33", show_default=True, help="Phred quality score.")
@optgroup.option("--jaccard-clip", "-j", is_flag=True, help="Enable jaccard clipping.")
@optgroup.option(
    "--memory-gb", type=int, default=None,
    help="Trinity --max_memory cap in GiB. Defaults to 60 percent of "
         "currently-available memory (floored at 4 GiB).",
)
@optgroup.group("Re-run steps")
@optgroup.option("--run-star", "-A", is_flag=True, help="Force re-run STAR read mapping.")
@optgroup.option("--run-trinity", "-T", is_flag=True, help="Force re-run Trinity assembly.")
@optgroup.option("--run-pasa", "-P", is_flag=True, help="Force re-run PASA alignment.")
@force_option
def assemble(
    genome: Path,
    left: Path | None,
    right: Path | None,
    single: Path | None,
    min_intron: int,
    max_intron: int,
    phred: str,
    numcpu: int,
    strand_specific: str | None,
    align_mode: str,
    run_star: bool,
    run_trinity: bool,
    run_pasa: bool,
    jaccard_clip: bool,
    splice_permissive: bool,
    genetic_code: str,
    memory_gb: int | None,
    force: bool,
) -> None:
    """Assemble transcriptome from RNA-seq reads.

    \b
    Provide either paired-end reads (--left and --right together) or
    single-end reads (--single). If using paired-end reads, both --left
    and --right are required.
    """
    from eukan.assembly import run_assembly
    from eukan.assembly.orchestrator import force_steps_from_run_flags
    from eukan.settings import AssemblyConfig

    if not left and not right and not single:
        raise click.UsageError("Provide --left/--right (paired) or --single reads.")
    if (left or right) and not (left and right):
        raise click.UsageError("Paired-end mode requires both --left and --right.")

    if strand_specific:
        if single and strand_specific in ("RF", "FR"):
            raise click.UsageError(
                "Paired-end strand types (RF/FR) cannot be used with single-end reads."
            )
        if (left or right) and strand_specific in ("R", "F"):
            raise click.UsageError(
                "Single-end strand types (R/F) cannot be used with paired-end reads."
            )

    if memory_gb is not None and memory_gb < 1:
        raise click.UsageError("--memory-gb must be at least 1 GiB.")

    config = AssemblyConfig(**drop_none(
        genome=genome.resolve(),
        work_dir=Path.cwd(),
        min_intron_len=min_intron,
        max_intron_len=max_intron,
        phred_quality=int(phred),
        num_cpu=numcpu,
        align_mode=align_mode,
        jaccard_clip=jaccard_clip,
        splice_permissive=splice_permissive,
        genetic_code=genetic_code,
        left_reads=left.resolve() if left else None,
        right_reads=right.resolve() if right else None,
        single_reads=single.resolve() if single else None,
        strand_specific=strand_specific,
        memory_gb=memory_gb,
    ))

    force_steps = force_steps_from_run_flags(
        run_star=run_star, run_trinity=run_trinity, run_pasa=run_pasa, force=force,
    )
    run_assembly(config, force_steps=force_steps or None)
    click.echo("Done.")
