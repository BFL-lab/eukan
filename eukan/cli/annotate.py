"""eukan annotate — genome annotation pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    FULL_CODE_TABLE,
    PreformattedEpilogCommand,
    drop_none,
    genome_option,
    numcpu_option,
)


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@optgroup.group("Required input")
@genome_option(
    "Genome sequence in FASTA format. Must not contain lower-case nucleotides "
    "(the pipeline soft-masks repeats by converting to lower-case)."
)
@optgroup.option(
    "--proteins", "-p", required=True, multiple=True,
    type=click.Path(exists=True, path_type=Path),
    help="One or more protein FASTA files.",
)
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--kingdom", "-k",
    type=click.Choice(["fungus", "protist", "animal", "plant"], case_sensitive=False),
    help="Target organism kingdom (tunes predictor parameters).",
)
@numcpu_option
@optgroup.option(
    "--existing-augustus", type=str, default=None,
    help="Use pre-trained AUGUSTUS species parameters.",
)
@optgroup.option(
    "--weights", "-w", type=int, multiple=True, default=(2, 1, 3),
    show_default=True,
    help="Weights for evidence sources: protein, gene predictions, transcripts.",
)
@optgroup.option(
    "--code", "-C", type=int, default=11, show_default=True,
    help="NCBI genetic code table number.",
)
@optgroup.group("Override options")
@optgroup.option(
    "--transcripts-fasta", "-tf", type=click.Path(exists=True, path_type=Path),
    help="Override auto-discovered transcript FASTA.",
)
@optgroup.option(
    "--transcripts-gff", "-tg", type=click.Path(exists=True, path_type=Path),
    help="Override auto-discovered transcript GFF3.",
)
@optgroup.option(
    "--rnaseq-hints", "-r", type=click.Path(exists=True, path_type=Path),
    help="Override auto-discovered RNA-seq hints GFF.",
)
@optgroup.option("--strand-specific", is_flag=True, help="Transcripts are strand-oriented.")
@optgroup.option(
    "--utrs", type=click.Path(exists=True, path_type=Path),
    help="PASA SQLite database path for adding UTRs.",
)
@optgroup.option(
    "--splice-permissive", is_flag=True, default=False,
    help="Allow non-canonical splice sites (GC-AG, AT-AC). "
    "When assembly evidence exists, observed splice types are used automatically; "
    "otherwise enables blanket allowance in AUGUSTUS.",
)
@optgroup.group("Experimental")
@optgroup.option(
    "--spsp", is_flag=True, default=False,
    help="Build species-specific spaln parameters from transcripts (alternative to fitild).",
)
@optgroup.group("Re-run steps")
@optgroup.option("--run-genemark", is_flag=True, help="Force re-run GeneMark gene prediction.")
@optgroup.option("--run-prot-align", is_flag=True, help="Force re-run protein alignment (spaln/gth).")
@optgroup.option("--run-augustus", is_flag=True, help="Force re-run AUGUSTUS training and prediction.")
@optgroup.option("--run-snap", is_flag=True, help="Force re-run SNAP (and CodingQuarry) prediction.")
@optgroup.option("--run-consensus", is_flag=True, help="Force re-run EVM consensus model building.")
def annotate(
    genome: Path,
    proteins: tuple[Path, ...],
    transcripts_fasta: Path | None,
    transcripts_gff: Path | None,
    rnaseq_hints: Path | None,
    existing_augustus: str | None,
    strand_specific: bool,
    splice_permissive: bool,
    spsp: bool,
    numcpu: int,
    weights: tuple[int, ...],
    code: int,
    utrs: Path | None,
    kingdom: str | None,
    run_genemark: bool,
    run_prot_align: bool,
    run_augustus: bool,
    run_snap: bool,
    run_consensus: bool,
) -> None:
    """Run the genome annotation pipeline.

    \b
    When run in the same directory as `eukan assemble`, transcript evidence
    (FASTA, GFF3, RNA-seq hints) and strand-specificity are discovered
    automatically. A PASA database for UTR addition is also detected if
    present. Use the override options to supply your own files or to
    replace the auto-discovered values.
    """
    from eukan.annotation import run_annotation_pipeline
    from eukan.annotation.pipeline import force_steps_from_run_flags
    from eukan.infra.layout import step_work_dir
    from eukan.settings import PipelineConfig

    # Only pass fields explicitly set by the user; pydantic-settings
    # fills the rest from pyproject.toml / env vars / defaults.
    config = PipelineConfig(**drop_none(
        genome=genome.resolve(),
        proteins=[p.resolve() for p in proteins],
        work_dir=step_work_dir("annotate"),
        num_cpu=numcpu,
        genetic_code=str(code),
        weights=list(weights),
        strand_specific=strand_specific,
        allow_noncanonical_splice=splice_permissive,
        spaln_ssp=spsp,
        kingdom=kingdom or None,
        transcripts_fasta=transcripts_fasta.resolve() if transcripts_fasta else None,
        transcripts_gff=transcripts_gff.resolve() if transcripts_gff else None,
        rnaseq_hints=rnaseq_hints.resolve() if rnaseq_hints else None,
        utrs_db=utrs.resolve() if utrs else None,
    ))

    force_steps = force_steps_from_run_flags(
        spaln_ssp=spsp,
        run_genemark=run_genemark,
        run_prot_align=run_prot_align,
        run_augustus=run_augustus,
        run_snap=run_snap,
        run_consensus=run_consensus,
    )

    result = run_annotation_pipeline(config, force_steps=force_steps or None)
    click.echo(f"Done. Final annotation: {result}")
