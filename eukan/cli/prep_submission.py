"""eukan prep-submission — package annotated genome for NCBI via table2asn."""

from __future__ import annotations

import shlex
from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import PreformattedEpilogCommand


@click.command("prep-submission", cls=PreformattedEpilogCommand)
@optgroup.group("Required input")
@optgroup.option(
    "--template", "-t", required=True,
    type=click.Path(exists=True, path_type=Path),
    help="NCBI submission template (.sbt). Generate one at "
         "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
)
@optgroup.group("Source qualifiers")
@optgroup.option(
    "--organism", type=str, default=None,
    help="Organism scientific name (e.g. 'Homo sapiens'). "
         "Required unless --source-info is given.",
)
@optgroup.option(
    "--isolate", type=str, default=None,
    help="Isolate / strain identifier (optional).",
)
@optgroup.option(
    "--source-info", "-j", type=str, default=None,
    help="Raw -j string for table2asn (e.g. '[organism=Foo] [isolate=Bar] "
         "[country=Canada]'). Overrides --organism / --isolate when set.",
)
@optgroup.option(
    "--locus-tag-prefix", type=str, default=None,
    help="NCBI-registered locus tag prefix (required for new-genome submissions).",
)
@optgroup.group("Override options")
@optgroup.option(
    "--genome", "-g", type=click.Path(exists=True, path_type=Path), default=None,
    help="Genome FASTA. Auto-discovered from eukan-run.json when omitted.",
)
@optgroup.option(
    "--gff3", "-i", type=click.Path(exists=True, path_type=Path), default=None,
    help="Annotated GFF3. Defaults to final.mod.gff3 (or final.gff3) "
         "in the working directory.",
)
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--cleanup", "-c", type=str, default="befw", show_default=True,
    help="table2asn -c cleanup flags.",
)
@optgroup.option(
    "--mode", "-M", type=str, default="n", show_default=True,
    help="table2asn -M flatfile mode.",
)
@optgroup.option(
    "--assembly-type", "-a", type=str, default="r10k", show_default=True,
    help="table2asn -a assembly type / gap configuration.",
)
@optgroup.option(
    "--extra-args", type=str, default="",
    help="Extra table2asn arguments, shell-quoted "
         "(e.g. --extra-args '-split-dr -huge').",
)
@optgroup.option(
    "--cleanup-gff3/--no-cleanup-gff3", "cleanup_gff3", default=True, show_default=True,
    help="Pre-process the GFF3 (strip UniProt cruft, drop CDS-less mRNAs, "
         "cap inferences) before handing it to table2asn.",
)
@optgroup.group("Output options")
@optgroup.option(
    "--output", "-o", type=click.Path(path_type=Path), default=None,
    help="Output .sqn path. Defaults to <outdir>/<genome-stem>.sqn.",
)
@optgroup.option(
    "--outdir", "-d", type=click.Path(path_type=Path), default=None,
    help="Output directory for .sqn and validator reports. Defaults to ./submission.",
)
@optgroup.option(
    "--print-command", is_flag=True,
    help="Print the resolved table2asn command and exit (no outputs written).",
)
@optgroup.option(
    "--dry-run", is_flag=True,
    help="Print the command and create the output directory, but don't run table2asn.",
)
def prep_submission(
    template: Path,
    organism: str | None,
    isolate: str | None,
    source_info: str | None,
    locus_tag_prefix: str | None,
    genome: Path | None,
    gff3: Path | None,
    cleanup: str,
    mode: str,
    assembly_type: str,
    extra_args: str,
    cleanup_gff3: bool,
    output: Path | None,
    outdir: Path | None,
    print_command: bool,
    dry_run: bool,
) -> None:
    """Validate + package an annotated genome for NCBI submission via table2asn.

    \b
    Runs the standard NCBI submission recipe (-split-logs -W -J -Z -euk
    -T -V b plus -c/-M/-a) over the genome FASTA and annotated GFF3,
    producing a .sqn file ready for upload along with .val (validator),
    .dr (discrepancy), and .stats reports for iterative GFF3 refinement.

    \b
    Auto-discovers inputs from the current working directory:
      - genome:   from eukan-run.json (manifest)
      - gff3:     final.mod.gff3 (preferred) or final.gff3 (fallback)

    \b
    The .sbt submission template must be created via NCBI's web form:
      https://submit.ncbi.nlm.nih.gov/genbank/template/submission/

    \b
    Use --print-command to inspect the exact table2asn invocation, or
    --extra-args to append flags not exposed here.
    """
    from eukan.cli._framework import drop_none
    from eukan.settings import SubmissionConfig
    from eukan.submission import run_prep_submission

    config = SubmissionConfig(**drop_none(
        work_dir=Path.cwd(),
        template=template.resolve(),
        organism=organism,
        isolate=isolate,
        source_info=source_info,
        locus_tag_prefix=locus_tag_prefix,
        genome=genome.resolve() if genome else None,
        gff3=gff3.resolve() if gff3 else None,
        cleanup=cleanup,
        mode=mode,
        assembly_type=assembly_type,
        extra_args=shlex.split(extra_args) if extra_args else None,
        cleanup_gff3=cleanup_gff3,
        output=output.resolve() if output else None,
        outdir=outdir.resolve() if outdir else None,
    ))

    sqn = run_prep_submission(config, print_only=print_command, dry_run=dry_run)

    if print_command:
        from eukan.submission import build_command, shell_repr
        click.echo(shell_repr(build_command(config)))
        return

    if dry_run:
        click.echo(f"Dry-run: would write {sqn}")
        return

    click.echo(f"Done. Wrote {sqn}")
