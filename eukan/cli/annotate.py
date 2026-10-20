"""eukan annotate — genome annotation pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    FULL_CODE_TABLE,
    PreformattedEpilogCommand,
    code_option,
    genome_option,
    given,
    numcpu_option,
)
from eukan.infra.layout import step_work_dir


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@optgroup.group("Required input")
@genome_option(
    "Genome sequence in FASTA format. Must not contain lower-case nucleotides "
    "(the pipeline soft-masks repeats by converting to lower-case)."
)
@optgroup.option(
    "--proteins", "-p", required=True, multiple=True,
    type=click.Path(exists=True, resolve_path=True, path_type=Path),
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
    "--weights", "-w", type=int, multiple=True, default=(2, 1, 3),
    show_default=True,
    help="Weights for evidence sources: protein, gene predictions, transcripts.",
)
@optgroup.option(
    "--combinr-stringent-overlap", type=float, default=None,
    help="combinr --stringent-overlap for the consensus --alt-splice isoform grouping: "
    "two transcript isoforms attach to one gene "
    "only when their span overlap is >= this percent of the shorter (default 0 = any "
    "overlap). Raise it (e.g. 30) to keep tip-overlapping collinear neighbours separate.",
)
@optgroup.option(
    "--orf-patch", type=click.Choice(["supported", "all"]), default=None,
    help="Patch transcript ORFs that overlap no consensus gene into the final models, "
    "labelled eukan-orf-patch / support=transcript_orf: 'supported' only those whose CDS "
    "overlaps a protein-alignment CDS on the same strand, 'all' every one (default: off).",
)
@code_option(default=11)
@optgroup.group("Override options")
@optgroup.option(
    "--transcripts-fasta", "-tf", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help="Override auto-discovered transcript FASTA.",
)
@optgroup.option(
    "--transcripts-gff", "-tg", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help="Override auto-discovered transcript GFF3.",
)
@optgroup.option(
    "--rnaseq-hints", "-r", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help="Override auto-discovered RNA-seq hints GFF.",
)
@optgroup.option(
    "--combinr-path", type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    help="Path to the combinr binary (default: 'combinr' on PATH).",
)
@optgroup.option(
    "--splice-permissive", "allow_noncanonical_splice", is_flag=True, default=False,
    help="Allow rarer non-canonical splice sites in AUGUSTUS. With assembly evidence, "
    "every observed donor-acceptor pair seen at >= 10 junctions by >= 50 unique reads "
    "is allowed regardless of its share of all junctions (the default also requires "
    ">= 1%); without evidence, allows AT-AC.",
)
@optgroup.group("Experimental")
@optgroup.option(
    "--spsp", "spaln_ssp", is_flag=True, default=False,
    help="Build species-specific spaln parameters from transcripts (alternative to fitild).",
)
@optgroup.group("Re-run steps")
@optgroup.option("--run-genemark", is_flag=True, help="Force re-run GeneMark gene prediction.")
@optgroup.option("--run-prot-align", is_flag=True, help="Force re-run protein alignment (spaln/gth).")
@optgroup.option("--run-augustus", is_flag=True, help="Force re-run AUGUSTUS training and prediction.")
@optgroup.option("--run-snap", is_flag=True, help="Force re-run SNAP (and CodingQuarry) prediction.")
@optgroup.option("--run-consensus", is_flag=True, help="Force re-run consensus model building.")
def annotate(**opts: Any) -> None:
    """Run the genome annotation pipeline.

    \b
    When run in the same directory as `eukan assemble`, transcript evidence
    (FASTA, GFF3, RNA-seq hints) is discovered automatically. Use the override
    options to supply your own files or to replace the auto-discovered values.
    """
    from eukan.annotation import run_annotation_pipeline  # lazy: startup cost
    from eukan.annotation.pipeline import force_steps_from_run_flags  # lazy: startup cost
    from eukan.settings import PipelineConfig  # lazy: startup cost

    run_flags = {k: opts.pop(k) for k in list(opts) if k.startswith("run_")}
    # Only the options given on the command line; the config fills the rest
    # from EUKAN_* env vars or its defaults.
    config = PipelineConfig(
        **given(opts), work_dir=step_work_dir("annotate"), manifest_dir=Path.cwd(),
    )
    config.discover()

    force_steps = force_steps_from_run_flags(spaln_ssp=config.spaln_ssp, **run_flags)

    result = run_annotation_pipeline(config, force_steps=force_steps or None)
    click.echo(f"Done. Final annotation: {result}")
