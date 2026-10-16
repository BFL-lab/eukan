"""NCBI submission preparation via table2asn.

Wraps the table2asn validator with the standard recipe used to produce a
``.sqn`` file from a genome FASTA + annotation GFF3 + submitter template
(``.sbt``). Iterating on table2asn's ``.val`` / ``.dr`` reports is the
intended workflow for refining the GFF3 prior to NCBI upload.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

from eukan.exceptions import ConfigurationError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import SubmissionConfig
from eukan.submission.cleanup import clean_gff3_for_submission

log = get_logger(__name__)


# table2asn flags we always set — the eukan submission recipe.
# Users who need to override these can pass `--extra-args -- <flag> <value>`.
_FIXED_FLAGS: tuple[str, ...] = (
    "-split-logs",   # split warnings/errors into separate log files
    "-W",            # write warnings report
    "-J",            # contigs only (no chromosomes)
    "-Z",            # write discrepancy report
    "-euk",          # eukaryotic
    "-T",            # taxonomy lookup
    "-V", "b",       # validation: basic
)


def _build_source_info(config: SubmissionConfig) -> str:
    """Compose the ``-j`` source-info string from organism/isolate, or pass through.

    Raises ConfigurationError if neither ``source_info`` nor ``organism`` is set.
    """
    if config.source_info:
        return config.source_info

    if not config.organism:
        raise ConfigurationError(
            "Missing required source qualifier.",
            hint="Pass --organism (and optionally --isolate), "
                 "or supply --source-info '[organism=...] [isolate=...]'.",
        )

    parts = [f"[organism={config.organism}]"]
    if config.isolate:
        parts.append(f"[isolate={config.isolate}]")
    return " ".join(parts)


def build_command(config: SubmissionConfig) -> list[str]:
    """Render the table2asn invocation for *config* as a list of argv tokens.

    Pure function — no filesystem or subprocess side-effects. Validates
    that required source qualifiers are present.
    """
    genome, gff3, sqn = config.required_paths()
    cmd: list[str] = ["table2asn", *_FIXED_FLAGS]
    cmd += ["-M", config.mode]
    cmd += ["-c", config.cleanup]
    cmd += ["-a", config.assembly_type]
    cmd += ["-t", str(config.template)]
    cmd += ["-j", _build_source_info(config)]
    cmd += ["-i", str(genome)]
    cmd += ["-f", str(gff3)]
    cmd += ["-o", str(sqn)]
    cmd += ["-outdir", str(config.output_dir)]
    if config.locus_tag_prefix:
        cmd += ["-locus-tag-prefix", config.locus_tag_prefix]
    cmd += list(config.extra_args)
    return cmd


def shell_repr(cmd: list[str]) -> str:
    """Return a copy-pasteable shell representation of *cmd*."""
    return shlex.join(cmd)


_SEVERITY_RE = re.compile(r"^(FATAL|ERROR|WARNING|INFO):", re.MULTILINE)


def _parse_val_report(val_path: Path) -> dict[str, int]:
    """Count severity occurrences in a table2asn ``.val`` report.

    Returns a dict with FATAL/ERROR/WARNING/INFO keys (zero when absent).
    Returns all-zero counts if the file cannot be read.
    """
    counts = {"FATAL": 0, "ERROR": 0, "WARNING": 0, "INFO": 0}
    try:
        text = val_path.read_text(errors="replace")
    except OSError:
        return counts
    for match in _SEVERITY_RE.finditer(text):
        counts[match.group(1)] += 1
    return counts


def run_prep_submission(
    config: SubmissionConfig,
    *,
    print_only: bool = False,
    dry_run: bool = False,
) -> Path:
    """Run table2asn for *config* and return the resulting ``.sqn`` path.

    Args:
        print_only: Print the command (no outdir, no subprocess), exit 0.
        dry_run: Print the command and create outdir, but skip the run.

    Returns:
        Path to the produced ``.sqn`` file (or the planned path under
        ``--print-command``/``--dry-run``).
    """
    _, gff3, sqn = config.required_paths()

    if print_only:
        cmd = build_command(config)
        log.info("table2asn command: %s", shell_repr(cmd))
        return sqn

    config.output_dir.mkdir(parents=True, exist_ok=True)

    if config.cleanup_gff3:
        cleaned = config.output_dir / f"{gff3.stem}.cleaned.gff3"
        report = clean_gff3_for_submission(gff3, cleaned)
        log.info("GFF3 cleanup: %s", report.summary())
        log.info("Cleaned GFF3: %s", cleaned)
        config = config.model_copy(update={"gff3": cleaned})

    cmd = build_command(config)
    log.info("table2asn command: %s", shell_repr(cmd))

    if dry_run:
        log.info("Dry-run requested; skipping table2asn invocation.")
        return sqn

    # table2asn doesn't honour absolute paths in -outdir consistently; cd
    # to the parent so its log files land next to the .sqn output.
    run_cmd(cmd, cwd=config.work_dir)

    log.info("Wrote %s", sqn)

    val = sqn.with_suffix(".val")  # table2asn writes it next to the .sqn
    counts = _parse_val_report(val)
    if any(counts.values()):
        log.info(
            "Validator: %d FATAL, %d ERROR, %d WARNING, %d INFO (full report: %s)",
            counts["FATAL"], counts["ERROR"], counts["WARNING"], counts["INFO"], val,
        )
    if counts["FATAL"]:
        raise ConfigurationError(
            f"table2asn reported {counts['FATAL']} FATAL validation error(s).",
            hint=f"See {val} for details.",
        )

    return sqn
