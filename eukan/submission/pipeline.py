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
    assert config.genome is not None  # post-discovery invariant
    assert config.gff3 is not None
    assert config.output is not None

    cmd: list[str] = ["table2asn", *_FIXED_FLAGS]
    cmd += ["-M", config.mode]
    cmd += ["-c", config.cleanup]
    cmd += ["-a", config.assembly_type]
    cmd += ["-t", str(config.template)]
    cmd += ["-j", _build_source_info(config)]
    cmd += ["-i", str(config.genome)]
    cmd += ["-f", str(config.gff3)]
    cmd += ["-o", str(config.output)]
    cmd += ["-outdir", str(config.outdir)]
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


def _val_path(config: SubmissionConfig) -> Path:
    """Path to the .val report table2asn writes alongside the .sqn output."""
    assert config.output is not None
    return config.output.with_suffix(".val")


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
    assert config.output is not None  # post-discovery invariant

    if print_only:
        cmd = build_command(config)
        log.info("table2asn command: %s", shell_repr(cmd))
        return config.output

    config.outdir.mkdir(parents=True, exist_ok=True)

    if config.cleanup_gff3:
        assert config.gff3 is not None
        cleaned = config.outdir / f"{config.gff3.stem}.cleaned.gff3"
        report = clean_gff3_for_submission(config.gff3, cleaned)
        log.info("GFF3 cleanup: %s", report.summary())
        log.info("Cleaned GFF3: %s", cleaned)
        object.__setattr__(config, "gff3", cleaned)

    cmd = build_command(config)
    log.info("table2asn command: %s", shell_repr(cmd))

    if dry_run:
        log.info("Dry-run requested; skipping table2asn invocation.")
        return config.output

    # table2asn doesn't honour absolute paths in -outdir consistently; cd
    # to the parent so its log files land next to the .sqn output.
    run_cmd(cmd, cwd=config.work_dir)

    log.info("Wrote %s", config.output)

    counts = _parse_val_report(_val_path(config))
    if any(counts.values()):
        log.info(
            "Validator: %d FATAL, %d ERROR, %d WARNING, %d INFO (full report: %s)",
            counts["FATAL"], counts["ERROR"], counts["WARNING"], counts["INFO"],
            _val_path(config),
        )
    if counts["FATAL"]:
        raise ConfigurationError(
            f"table2asn reported {counts['FATAL']} FATAL validation error(s).",
            hint=f"See {_val_path(config)} for details.",
        )

    return config.output
