"""Click CLI for the eukan genome annotation pipeline."""

from __future__ import annotations

import multiprocessing
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import click
from click_option_group import optgroup

try:
    _EUKAN_VERSION = version("eukan")
except PackageNotFoundError:
    _EUKAN_VERSION = "unknown"


class _PreformattedEpilogCommand(click.Command):
    """Click command with preformatted epilog shown only on --help, not -h.

    Also renders option-group help without the wrapping "Options:" header
    and without the extra indentation that click-option-group adds to
    grouped options.

    Epilog values are resolved through :func:`_resolve_epilog` so heavy
    text (e.g. the genetic-code tables, which import BioPython) can be
    deferred until ``--help`` is actually requested.
    """

    def format_epilog(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        if self.epilog and _show_full_help:
            text = _resolve_epilog(self.epilog)
            formatter.write("\n")
            for line in text.split("\n"):
                formatter.write(f"  {line}\n")

    def format_options(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        from click_option_group import OptionGroup

        # Group help records into sections separated by option-group headers.
        sections: list[list[tuple[str, str]]] = [[]]
        for param in self.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is None:
                continue
            if type(param).__name__ == "_GroupTitleFakeOption":
                # Start a new section for each group header.
                if sections[-1]:
                    sections.append([])
                sections[-1].append(rv)
            elif isinstance(getattr(param, "group", None), OptionGroup):
                sections[-1].append((rv[0].lstrip(), rv[1]))
            else:
                # Ungrouped options (e.g. --help) get their own section.
                if sections[-1]:
                    sections.append([])
                sections[-1].append(rv)

        for section in sections:
            if section:
                formatter.write("\n")
                formatter.write_dl(section)


# Track whether --help (verbose) or -h (brief) was used
_show_full_help = False


class _EukanGroup(click.Group):
    """Click group with help mode tracking and structured error formatting.

    Catches :class:`~eukan.exceptions.EukanError` at the CLI boundary and
    prints user-friendly messages instead of raw Python tracebacks.
    """

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        global _show_full_help
        _show_full_help = "--help" in args
        return super().parse_args(ctx, args)

    def invoke(self, ctx: click.Context):
        try:
            return super().invoke(ctx)
        except click.exceptions.Exit:
            raise
        except click.exceptions.Abort:
            raise
        except SystemExit:
            raise
        except Exception as exc:
            from pydantic import ValidationError as PydanticValidationError

            if isinstance(exc, PydanticValidationError):
                click.secho("Error: invalid configuration", fg="red", err=True)
                for err in exc.errors():
                    loc = " → ".join(str(part) for part in err["loc"])
                    click.echo(f"  {loc}: {err['msg']}", err=True)
                raise SystemExit(1) from exc

            from eukan.exceptions import (
                ConfigurationError,
                DependencyError,
                EukanError,
                ExternalToolError,
                ValidationError,
            )
            if not isinstance(exc, EukanError):
                raise

            if isinstance(exc, ExternalToolError):
                click.secho(f"Error: {exc.tool} failed (exit {exc.returncode})", fg="red", err=True)
                if exc.step:
                    click.echo(f"  Step: {exc.step}", err=True)
                if exc.stderr_snippet:
                    click.echo(f"  stderr: {exc.stderr_snippet[:300]}", err=True)
                click.echo("  Run with -v for full command and stderr output.", err=True)
            elif isinstance(exc, (ValidationError, DependencyError, ConfigurationError)):
                click.secho(f"Error: {exc}", fg="red", err=True)
            else:
                click.secho(f"Error: {exc}", fg="red", err=True)

            if exc.hint:
                click.echo(f"  Hint: {exc.hint}", err=True)

            raise SystemExit(1) from exc


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


# Shared option decorators reused across subcommands.

def _numcpu_option(func):
    return optgroup.option(
        "--numcpu", "-n", type=int, default=multiprocessing.cpu_count(),
        show_default=True, help="Number of CPU threads.",
    )(func)


def _force_option(func):
    return optgroup.option(
        "--force", "-f", is_flag=True,
        help="Force re-run all steps (ignore cached outputs).",
    )(func)


def _genome_option(help_text: str = "Genome sequence in FASTA format."):
    """Build a --genome/-g required-path option with custom help text."""
    def decorator(func):
        return optgroup.option(
            "--genome", "-g", required=True,
            type=click.Path(exists=True, path_type=Path),
            help=help_text,
        )(func)
    return decorator


# ---------------------------------------------------------------------------
# Genetic code tables (deferred -- importing BioPython is slow)
# ---------------------------------------------------------------------------

# Sentinel epilog markers; the real text is built lazily on first --help.
_FULL_CODE_TABLE = "__EUKAN_FULL_CODE_TABLE__"
_PASA_CODE_TABLE = "__EUKAN_PASA_CODE_TABLE__"


def _resolve_epilog(value: str) -> str:
    """Resolve an epilog sentinel to its rendered text. Pass-through otherwise."""
    if value == _FULL_CODE_TABLE:
        return _full_code_table_text()
    if value == _PASA_CODE_TABLE:
        return _pasa_code_table_text()
    return value


def _full_code_table_text() -> str:
    from Bio.Data import CodonTable

    from eukan.gencode import _PASA_NAMES

    lines = ["Genetic codes (NCBI translation tables):", ""]
    for cid, table in sorted(CodonTable.unambiguous_dna_by_id.items()):
        marker = " *" if cid in _PASA_NAMES else ""
        lines.append(f"  {cid:>2}  {table.names[0]}{marker}")
    lines.append("")
    lines.append("  * = also supported by PASA (eukan assemble)")
    return "\n".join(lines)


def _pasa_code_table_text() -> str:
    from eukan.gencode import _PASA_NAMES, GeneticCode

    lines = ["Genetic codes supported by PASA:", ""]
    for cid in sorted(_PASA_NAMES):
        gc = GeneticCode(cid)
        ncbi_name = gc.codon_table.names[0]
        lines.append(f"  {cid:>2}  {ncbi_name} ({gc.pasa_name})")
    return "\n".join(lines)


@click.group(cls=_EukanGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(_EUKAN_VERSION)
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress info logging (warnings only).")
def cli(verbose: bool, quiet: bool) -> None:
    """Eukan: Eukaryotic nuclear genome annotation pipeline.

    \b
    Typical workflow:
      1. eukan check        Verify installation and tools
         eukan db-fetch     Download UniProt + Pfam databases
      2. eukan assemble     Build transcriptome from RNA-seq (optional but recommended)
      3. eukan annotate     Annotate genome using proteins + assembly (if available)
      4. eukan func-annot   Add functional info to predicted proteins
    \b
    Helpers:
         eukan compare      Compare annotations against a reference/previous annotation
         eukan gff3toseq    Extract sequences from annotation
         eukan status       View progress of a pipeline run
    """
    from eukan.infra.environ import configure_process_env
    from eukan.infra.logging import setup_logging

    verbosity = 1 if verbose else (-1 if quiet else 0)
    setup_logging(verbosity)
    configure_process_env()


# ---------------------------------------------------------------------------
# eukan annotate
# ---------------------------------------------------------------------------


@cli.command(cls=_PreformattedEpilogCommand, epilog=_FULL_CODE_TABLE)
@optgroup.group("Required input")
@_genome_option(
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
@_numcpu_option
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
    from eukan.settings import PipelineConfig

    # Only pass fields explicitly set by the user; pydantic-settings
    # fills the rest from pyproject.toml / env vars / defaults.
    kwargs: dict = {
        "genome": genome.resolve(),
        "proteins": [p.resolve() for p in proteins],
        "work_dir": Path.cwd(),
        "num_cpu": numcpu,
        "genetic_code": str(code),
        "weights": list(weights),
        "strand_specific": strand_specific,
        "allow_noncanonical_splice": splice_permissive,
        "spaln_ssp": spsp,
    }
    if kingdom:
        kwargs["kingdom"] = kingdom
    if transcripts_fasta:
        kwargs["transcripts_fasta"] = transcripts_fasta.resolve()
    if transcripts_gff:
        kwargs["transcripts_gff"] = transcripts_gff.resolve()
    if rnaseq_hints:
        kwargs["rnaseq_hints"] = rnaseq_hints.resolve()
    if utrs:
        kwargs["utrs_db"] = utrs.resolve()

    config = PipelineConfig(**kwargs)

    from eukan.annotation.orchestrator import force_steps_from_run_flags

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


# ---------------------------------------------------------------------------
# eukan assemble
# ---------------------------------------------------------------------------


@cli.command(cls=_PreformattedEpilogCommand, epilog=_PASA_CODE_TABLE)
@optgroup.group("Required input")
@_genome_option("Genome FASTA file.")
@optgroup.option("--left", "-l", type=click.Path(exists=True, path_type=Path), help="Left paired-end reads.")
@optgroup.option("--right", "-r", type=click.Path(exists=True, path_type=Path), help="Right paired-end reads.")
@optgroup.option("--single", "-s", type=click.Path(exists=True, path_type=Path), help="Single-end reads.")
@optgroup.group("Pipeline parameters")
@_numcpu_option
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
@optgroup.group("Re-run steps")
@optgroup.option("--run-star", "-A", is_flag=True, help="Force re-run STAR read mapping.")
@optgroup.option("--run-trinity", "-T", is_flag=True, help="Force re-run Trinity assembly.")
@optgroup.option("--run-pasa", "-P", is_flag=True, help="Force re-run PASA alignment.")
@_force_option
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
    force: bool,
) -> None:
    """Assemble transcriptome from RNA-seq reads.

    \b
    Provide either paired-end reads (--left and --right together) or
    single-end reads (--single). If using paired-end reads, both --left
    and --right are required.
    """
    from eukan.assembly import run_assembly
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

    kwargs: dict = {
        "genome": genome.resolve(),
        "work_dir": Path.cwd(),
        "min_intron_len": min_intron,
        "max_intron_len": max_intron,
        "phred_quality": int(phred),
        "num_cpu": numcpu,
        "align_mode": align_mode,
        "jaccard_clip": jaccard_clip,
        "splice_permissive": splice_permissive,
        "genetic_code": genetic_code,
    }
    if left:
        kwargs["left_reads"] = left.resolve()
    if right:
        kwargs["right_reads"] = right.resolve()
    if single:
        kwargs["single_reads"] = single.resolve()
    if strand_specific:
        kwargs["strand_specific"] = strand_specific

    config = AssemblyConfig(**kwargs)

    steps = []
    if run_star:
        steps.append("map")
    if run_trinity:
        steps.append("trinity")
    if run_pasa:
        steps.append("pasa")

    # If no specific steps or --force, run all
    if not steps:
        steps = ["map", "trinity", "pasa"]

    run_assembly(config, steps, force=force)
    click.echo("Done.")


# ---------------------------------------------------------------------------
# eukan func-annot
# ---------------------------------------------------------------------------


@cli.command("func-annot", cls=_PreformattedEpilogCommand)
@optgroup.group("Pipeline parameters")
@_numcpu_option
@optgroup.option("--evalue", "-e", type=str, default="1e-1", show_default=True, help="E-value cutoff.")
@optgroup.group("Override options")
@optgroup.option(
    "--proteins", "-p", type=click.Path(exists=True, path_type=Path),
    help="Amino acid sequences in FASTA format.",
)
@optgroup.option(
    "--uniprot", type=click.Path(exists=True, path_type=Path),
    default=None, help="UniProt-SwissProt database FASTA.",
)
@optgroup.option(
    "--pfam", type=click.Path(exists=True, path_type=Path),
    default=None, help="Pfam HMM database.",
)
@optgroup.option(
    "--gff3", type=click.Path(exists=True, path_type=Path),
    default=None, help="GFF3 file to annotate with functional info.",
)
@_force_option
def func_annot(
    proteins: Path,
    uniprot: Path | None,
    pfam: Path | None,
    gff3: Path | None,
    numcpu: int,
    evalue: str,
    force: bool,
) -> None:
    """Add functional annotations (UniProt + Pfam) to proteins.

    \b
    When run after `eukan annotate` and `eukan db-fetch`, the predicted
    protein sequences, UniProt, and Pfam databases are discovered
    automatically. Use the override options to point to different files
    or to run functional annotation independently of the main pipeline.
    """
    from eukan.functional import run_functional_annotation

    if proteins is None:
        raise click.UsageError(
            "No protein file found. Provide --proteins or run `eukan annotate` first."
        )

    run_functional_annotation(
        proteins=proteins.resolve(),
        uniprot_db=uniprot.resolve() if uniprot else None,
        pfam_db=pfam.resolve() if pfam else None,
        gff3_path=gff3.resolve() if gff3 else None,
        num_cpu=numcpu,
        evalue=evalue,
        force=force,
    )
    click.echo("Done.")


# ---------------------------------------------------------------------------
# eukan gff3toseq
# ---------------------------------------------------------------------------


@cli.command(cls=_PreformattedEpilogCommand, epilog=_FULL_CODE_TABLE)
@click.option(
    "--genome", "-g", required=True, type=click.Path(exists=True, path_type=Path),
    help="Genome assembly in FASTA format.",
)
@click.option(
    "--gff3", "-i", required=True, type=click.Path(exists=True, path_type=Path),
    help="GFF3 file with gene models.",
)
@click.option(
    "--output", "-o", type=click.Choice(["protein", "cdna"]), default="protein",
    show_default=True, help="Output sequence type.",
)
@click.option(
    "--code", "-c", type=int, default=1, show_default=True,
    help="NCBI genetic code table number.",
)
def gff3toseq(genome: Path, gff3: Path, output: str, code: int) -> None:
    """Extract protein or cDNA sequences from GFF3 + genome."""

    from eukan.gff.io import extract_sequences

    for record in extract_sequences(gff3, genome, extract_to=output, genetic_code=code):
        sys.stdout.write(record.format("fasta"))


# ---------------------------------------------------------------------------
# eukan db-fetch
# ---------------------------------------------------------------------------


@cli.command("db-fetch")
@click.option(
    "--output-dir", "-o", type=click.Path(path_type=Path), default="databases",
    show_default=True, help="Directory to download databases into.",
)
@click.option(
    "--force", "-f", is_flag=True, help="Re-download even if databases are up to date.",
)
@click.option(
    "--database", "-d", multiple=True,
    type=click.Choice(["uniprot", "pfam"], case_sensitive=False),
    help="Specific database(s) to fetch. If omitted, fetch all.",
)
def db_fetch(output_dir: Path, force: bool, database: tuple[str, ...]) -> None:
    """Download reference databases (UniProt, Pfam)."""
    from eukan.functional.dbfetch import fetch_databases

    output_dir = output_dir.resolve()
    fetch_databases(
        output_dir,
        force=force,
        databases=list(database) if database else None,
    )
    click.echo("Done.")


# ---------------------------------------------------------------------------
# eukan check
# ---------------------------------------------------------------------------


@cli.command()
@click.option(
    "--for", "subcommands", multiple=True,
    type=click.Choice(["annotate", "assemble", "func-annot", "db-fetch"], case_sensitive=False),
    help="Only check tools needed by these subcommands. If omitted, check all.",
)
@click.option(
    "--db-dir", type=click.Path(path_type=Path), default="databases",
    show_default=True, help="Database directory to check.",
)
def check(subcommands: tuple[str, ...], db_dir: Path) -> None:
    """Verify Python deps, external tools, and databases."""
    from eukan.check import format_results, run_checks

    passed, failed, db_results, python_results = run_checks(
        list(subcommands) if subcommands else None,
        db_dir=db_dir.resolve(),
    )
    click.echo(format_results(passed, failed, db_results, python_results))

    py_failures = any(not r.ok for r in python_results)
    db_failures = any(not ok for _, _, ok in db_results)
    if failed or db_failures or py_failures:
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# eukan status
# ---------------------------------------------------------------------------


@cli.command()
@click.option(
    "--work-dir", "-d", type=click.Path(exists=True, path_type=Path), default=".",
    help="Working directory containing eukan-run.json.",
)
def status(work_dir: Path) -> None:
    """Show the status of a pipeline run."""
    from eukan.infra.manifest import format_status, load_manifest

    manifest = load_manifest(work_dir.resolve())
    if not manifest:
        click.echo("No pipeline run found in this directory (no eukan-run.json).")
        raise SystemExit(1)

    click.echo(format_status(manifest))


# ---------------------------------------------------------------------------
# eukan compare
# ---------------------------------------------------------------------------


@cli.command(cls=_PreformattedEpilogCommand)
@optgroup.group("Required input")
@optgroup.option(
    "--reference", "-r", required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference GFF3 file.",
)
@optgroup.option(
    "--predicted", "-p", required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Predicted GFF3 file to evaluate.",
)
@optgroup.group("Output options")
@optgroup.option(
    "--output-file", "-o", type=click.Path(path_type=Path), default=None,
    help="Write per-feature details to a TSV file.",
)
def compare(reference: Path, predicted: Path, output_file: Path | None) -> None:
    """Compare predicted gene models against a reference GFF3. The reference can
    be manually curated models, or predictions from another tool or pipeline.

    \b
    Computes quality metrics at gene, mRNA, CDS, and intron levels.
    Use --output-file to write a per-feature TSV for further analysis.
    """
    from eukan.stats import compare_annotations, format_results, write_details_tsv

    result = compare_annotations(reference.resolve(), predicted.resolve())
    click.echo(format_results(result))
    if output_file is not None:
        tsv_path = write_details_tsv(result, output_file.resolve())
        click.echo(f"\nDetails written to {tsv_path}")
