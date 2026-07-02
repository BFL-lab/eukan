"""rnaSPAdes de novo transcriptome assembly."""

from __future__ import annotations

import shutil

from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def _reads_args(config: AssemblyConfig) -> list[str]:
    """rnaSPAdes read arguments: ``-1``/``-2`` for paired, ``-s`` for single."""
    if config.left_reads and config.right_reads:
        return ["-1", str(config.left_reads), "-2", str(config.right_reads)]
    if config.single_reads:
        return ["-s", str(config.single_reads)]
    raise ValueError("No read files provided")


def run_rnaspades(config: AssemblyConfig) -> None:
    """Run rnaSPAdes de novo assembly, normalizing output to ``rnaspades.fasta``.

    Skips when ``rnaspades.fasta`` already exists. rnaSPAdes writes a tree of
    intermediate files under its ``-o`` directory; we keep only the final
    ``transcripts.fasta`` (renamed to ``rnaspades.fasta``) and discard the rest.
    """
    wd = config.work_dir
    final = wd / "rnaspades.fasta"
    if final.exists():
        return

    out_dir = wd / "rnaspades_out"
    log.info("Running rnaSPAdes de novo assembly...")
    run_cmd(
        [
            "rnaspades.py",
            *_reads_args(config),
            "-t", str(config.num_cpu),
            "-m", str(config.memory_gb),
            "--phred-offset", str(config.phred_quality),
            "-o", str(out_dir),
        ],
        cwd=wd,
    )
    produced = out_dir / "transcripts.fasta"
    if produced.exists():
        shutil.move(str(produced), str(final))
    shutil.rmtree(out_dir, ignore_errors=True)
