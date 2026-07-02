"""Shared CLI infrastructure: option-group rendering, help-mode tracking,
structured error formatting, and reusable option decorators.

The ``cli`` group itself lives in :mod:`eukan.cli`. Per-subcommand modules
import only the helpers they need from here.
"""

from __future__ import annotations

import errno
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


def _resolve_epilog(value: str) -> str:
    """Resolve an epilog sentinel to its rendered text. Pass-through otherwise."""
    if value == FULL_CODE_TABLE:
        return _full_code_table_text()
    return value


def _full_code_table_text() -> str:
    from Bio.Data import CodonTable

    lines = ["Genetic codes (NCBI translation tables):", ""]
    for cid, table in sorted(CodonTable.unambiguous_dna_by_id.items()):
        lines.append(f"  {cid:>2}  {table.names[0]}")
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

            disk_full = _format_disk_full(exc)
            if disk_full is not None:
                title, details = disk_full
                click.secho(title, fg="red", err=True)
                for line in details:
                    click.echo(f"  {line}", err=True)
                raise SystemExit(1) from exc

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


_DISK_FULL_PATTERNS = (
    "no space left on device",
    "disk quota exceeded",
)


def _format_disk_full(exc: BaseException) -> tuple[str, list[str]] | None:
    """Detect ENOSPC conditions and return a CLI-formatted (title, details).

    Triggers on a direct ``OSError(errno=ENOSPC)`` from Python file ops,
    and on an ``ExternalToolError`` whose captured stderr indicates the
    underlying tool ran out of disk space. Returns ``None`` when the
    exception is unrelated to disk-full.
    """
    if isinstance(exc, OSError) and exc.errno in (errno.ENOSPC, errno.EDQUOT):
        title = "Error: no space left on device"
        details: list[str] = []
        if exc.filename:
            details.append(f"Path: {exc.filename}")
        details.append("Free up disk space on the device hosting the work directory and re-run.")
        return title, details

    from eukan.exceptions import ExternalToolError
    if isinstance(exc, ExternalToolError):
        haystack = (exc.stderr_snippet or "").lower()
        if any(p in haystack for p in _DISK_FULL_PATTERNS):
            details = [f"Tool: {exc.tool} (exit {exc.returncode})"]
            if exc.step:
                details.append(f"Step: {exc.step}")
            details.append("Free up disk space on the device hosting the work directory and re-run.")
            return "Error: no space left on device", details

    return None


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


def resolve_optional_path(path: Path | None) -> Path | None:
    """Resolve an optional CLI path to an absolute path, preserving ``None``.

    Centralizes the ``x.resolve() if x else None`` idiom repeated when
    plumbing optional path arguments into configs.
    """
    return path.resolve() if path else None


def numcpu_option(func):
    return optgroup.option(
        "--numcpu", "-n", type=int, default=multiprocessing.cpu_count(),
        show_default=True, help="Number of CPU threads.",
    )(func)


def force_option(
    func=None,
    *,
    help_text: str = "Force re-run all steps (ignore cached outputs).",
):
    """Build a --force/-f flag option, optionally with custom help text.

    Usable as a bare decorator (``@force_option``) for the default help,
    or as a factory (``@force_option(help_text="…")``) for command-specific
    wording (e.g. db-fetch reads "Re-download…", not "Force re-run…").
    """
    def decorator(f):
        return optgroup.option(
            "--force", "-f", is_flag=True, help=help_text,
        )(f)

    if func is None:
        return decorator
    return decorator(func)


def genome_option(help_text: str = "Genome sequence in FASTA format."):
    """Build a --genome/-g required-path option with custom help text."""
    def decorator(func):
        return optgroup.option(
            "--genome", "-g", required=True,
            type=click.Path(exists=True, path_type=Path),
            help=help_text,
        )(func)
    return decorator


def code_option(default: int = 1):
    """Build a --code/-c NCBI genetic-code option with a command-specific default."""
    def decorator(func):
        return optgroup.option(
            "--code", "-c", type=int, default=default, show_default=True,
            help="NCBI genetic code table number.",
        )(func)
    return decorator
