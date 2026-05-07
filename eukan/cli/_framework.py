"""Shared CLI infrastructure: option-group rendering, help-mode tracking,
structured error formatting, and reusable option decorators.

The ``cli`` group itself lives in :mod:`eukan.cli`. Per-subcommand modules
import only the helpers they need from here.
"""

from __future__ import annotations

import multiprocessing
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import click
from click_option_group import optgroup

try:
    EUKAN_VERSION = version("eukan")
except PackageNotFoundError:
    EUKAN_VERSION = "unknown"


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


# ---------------------------------------------------------------------------
# Help mode + epilog plumbing
# ---------------------------------------------------------------------------


# Sentinel epilog markers; the real text is built lazily on first --help so
# importing BioPython for the genetic-code tables doesn't tax -h users.
FULL_CODE_TABLE = "__EUKAN_FULL_CODE_TABLE__"
PASA_CODE_TABLE = "__EUKAN_PASA_CODE_TABLE__"


def _resolve_epilog(value: str) -> str:
    """Resolve an epilog sentinel to its rendered text. Pass-through otherwise."""
    if value == FULL_CODE_TABLE:
        return _full_code_table_text()
    if value == PASA_CODE_TABLE:
        return _pasa_code_table_text()
    return value


def _full_code_table_text() -> str:
    from Bio.Data import CodonTable

    from eukan.infra.genetic_code import _PASA_NAMES

    lines = ["Genetic codes (NCBI translation tables):", ""]
    for cid, table in sorted(CodonTable.unambiguous_dna_by_id.items()):
        marker = " *" if cid in _PASA_NAMES else ""
        lines.append(f"  {cid:>2}  {table.names[0]}{marker}")
    lines.append("")
    lines.append("  * = also supported by PASA (eukan assemble)")
    return "\n".join(lines)


def _pasa_code_table_text() -> str:
    from eukan.infra.genetic_code import _PASA_NAMES, GeneticCode

    lines = ["Genetic codes supported by PASA:", ""]
    for cid in sorted(_PASA_NAMES):
        gc = GeneticCode(cid)
        ncbi_name = gc.codon_table.names[0]
        lines.append(f"  {cid:>2}  {ncbi_name} ({gc.pasa_name})")
    return "\n".join(lines)


@dataclass
class CLIState:
    """Per-invocation CLI state, attached to ``ctx.obj``.

    Replaces a previous module global so concurrent invocations (tests,
    programmatic use) don't clobber each other.
    """

    show_full_help: bool = False


class PreformattedEpilogCommand(click.Command):
    """Click command with preformatted epilog shown only on --help, not -h.

    Also renders option-group help without the wrapping "Options:" header
    and without the extra indentation that click-option-group adds to
    grouped options.

    Epilog values are resolved through :func:`_resolve_epilog` so heavy
    text (e.g. the genetic-code tables, which import BioPython) can be
    deferred until ``--help`` is actually requested.
    """

    def format_epilog(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        # ctx.obj.show_full_help is set by EukanGroup.parse_args based on
        # whether --help (verbose) or -h (brief) was used. Walk up the
        # context chain to find it (subcommand contexts inherit obj).
        ctx_walk: click.Context | None = ctx
        while ctx_walk is not None:
            obj = getattr(ctx_walk, "obj", None)
            if isinstance(obj, CLIState) and obj.show_full_help:
                text = _resolve_epilog(self.epilog or "")
                formatter.write("\n")
                for line in text.split("\n"):
                    formatter.write(f"  {line}\n")
                return
            ctx_walk = ctx_walk.parent

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


class EukanGroup(click.Group):
    """Click group with help mode tracking and structured error formatting.

    Catches :class:`~eukan.exceptions.EukanError` at the CLI boundary and
    prints user-friendly messages instead of raw Python tracebacks.
    """

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        # ensure_object so ctx.obj is set even before any callback runs.
        ctx.ensure_object(CLIState)
        ctx.obj.show_full_help = "--help" in args
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

            from eukan.exceptions import EukanError
            if not isinstance(exc, EukanError):
                raise

            title, details = exc.format_for_cli()
            click.secho(title, fg="red", err=True)
            for line in details:
                click.echo(f"  {line}", err=True)
            if exc.hint:
                click.echo(f"  Hint: {exc.hint}", err=True)

            raise SystemExit(1) from exc


# ---------------------------------------------------------------------------
# Helpers shared across subcommands
# ---------------------------------------------------------------------------


def drop_none(**kwargs) -> dict:
    """Build a kwargs dict, dropping keys whose value is ``None``.

    Keeps the CLI → config plumbing terse without losing the ability for
    pydantic-settings to discover unset fields from pyproject.toml or
    environment variables: an explicit ``None`` from a missing CLI flag
    must not overwrite a value those sources would have provided.
    """
    return {k: v for k, v in kwargs.items() if v is not None}


def numcpu_option(func):
    return optgroup.option(
        "--numcpu", "-n", type=int, default=multiprocessing.cpu_count(),
        show_default=True, help="Number of CPU threads.",
    )(func)


def force_option(func):
    return optgroup.option(
        "--force", "-f", is_flag=True,
        help="Force re-run all steps (ignore cached outputs).",
    )(func)


def genome_option(help_text: str = "Genome sequence in FASTA format."):
    """Build a --genome/-g required-path option with custom help text."""
    def decorator(func):
        return optgroup.option(
            "--genome", "-g", required=True,
            type=click.Path(exists=True, path_type=Path),
            help=help_text,
        )(func)
    return decorator
