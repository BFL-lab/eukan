"""Click CLI for the eukan genome annotation pipeline.

This package wires the per-subcommand modules to a top-level ``cli``
group. The pyproject entry point is ``eukan = "eukan.cli:cli"``, which
resolves to :data:`cli` below.
"""

from __future__ import annotations

import click

from eukan.cli._framework import (
    CONTEXT_SETTINGS,
    EUKAN_VERSION,
    EukanGroup,
)
from eukan.cli.annotate import annotate
from eukan.cli.assemble import assemble
from eukan.cli.check import check
from eukan.cli.compare import compare
from eukan.cli.db_fetch import db_fetch
from eukan.cli.func_annot import func_annot
from eukan.cli.gff3toseq import gff3toseq
from eukan.cli.mask_repeats import mask_repeats
from eukan.cli.prep_submission import prep_submission
from eukan.cli.status import status


@click.group(cls=EukanGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(EUKAN_VERSION)
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress info logging (warnings only).")
def cli(verbose: bool, quiet: bool) -> None:
    """Eukan: Eukaryotic nuclear genome annotation pipeline.

    \b
    Typical workflow:
      1. eukan check           Verify installation and tools
         eukan db-fetch        Download UniProt + Pfam databases
      2. eukan mask-repeats    Soft-mask repeats (optional but recommended)
      3. eukan assemble        Build transcriptome from RNA-seq (optional but recommended)
      4. eukan annotate        Annotate genome using proteins + assembly (if available)
      5. eukan func-annot      Add functional info to predicted proteins
      6. eukan prep-submission Validate + package for NCBI submission via table2asn
    \b
    Helpers:
         eukan compare      Compare annotations against a reference/previous annotation
         eukan gff3toseq    Extract sequences from annotation
         eukan status       View progress of a pipeline run
    """
    import signal

    from eukan.infra.environ import configure_process_env
    from eukan.infra.logging import setup_logging
    from eukan.infra.runner import terminate_all_children

    verbosity = 1 if verbose else (-1 if quiet else 0)
    setup_logging(verbosity)
    configure_process_env()

    def _on_sigint(signum, frame):
        # Tear down child subprocesses before propagating the interrupt.
        terminate_all_children()
        raise KeyboardInterrupt

    signal.signal(signal.SIGINT, _on_sigint)


for _cmd in (
    annotate, assemble, mask_repeats, func_annot, prep_submission,
    gff3toseq, db_fetch, check, status, compare,
):
    cli.add_command(_cmd)


__all__ = ["cli"]
