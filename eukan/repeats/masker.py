"""RepeatMasker softmasking + RepeatMasker → AUGUSTUS hint conversion."""

from __future__ import annotations

import shutil
from pathlib import Path

from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import RepeatsConfig

log = get_logger(__name__)


def gff_to_hints(repeat_gff: Path, hints_out: Path) -> None:
    """Convert RepeatMasker GFF into AUGUSTUS ``nonexonpart`` hints.

    Each non-comment, non-blank line becomes a ``nonexonpart`` feature
    sourced from RepeatMasker (``src=RM``), keyed for the ``RM`` extrinsic
    source already declared in ``data/configs/augustus.config``.
    """
    with open(repeat_gff) as fin, open(hints_out, "w") as fout:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                # RepeatMasker GFF rows always have 9 tab-delimited columns;
                # anything shorter is malformed and should be skipped.
                continue
            cols[2] = "nonexonpart"
            cols[8] = "src=RM"
            fout.write("\t".join(cols) + "\n")


def run_masker(config: RepeatsConfig, families: Path) -> tuple[Path, Path]:
    """Softmask the genome with RepeatMasker, using *families* as the library.

    Returns ``(masked_genome, hints_gff)`` paths in ``work_dir``.
    """
    sdir = config.work_dir / "repeats" / "masker"
    sdir.mkdir(parents=True, exist_ok=True)

    # RepeatMasker writes its outputs alongside the input genome by default.
    # Stage a symlink-free copy in sdir so we don't pollute work_dir with
    # *.cat / *.tbl artifacts and so reruns can be cleaned by deleting sdir.
    staged_genome = sdir / config.genome.name
    if not staged_genome.exists() or staged_genome.is_symlink():
        staged_genome.unlink(missing_ok=True)
        shutil.copy2(config.genome, staged_genome)

    log.info("Running RepeatMasker (this may take many hours)...")
    run_cmd(
        [
            "RepeatMasker",
            "-xsmall",
            "-gff",
            "-engine", config.engine,
            "-s",
            "-pa", str(config.num_cpu),
            "-lib", str(families),
            staged_genome.name,
        ],
        cwd=sdir,
    )

    masked_src = sdir / f"{config.genome.name}.masked"
    repeats_gff_src = sdir / f"{config.genome.name}.out.gff"

    masked_dst = config.work_dir / f"{config.name}.masked.fasta"
    repeats_gff_dst = config.work_dir / f"{config.name}.repeats.gff"
    hints_dst = config.work_dir / "hints_repeatmask.gff"

    shutil.copy2(masked_src, masked_dst)
    shutil.copy2(repeats_gff_src, repeats_gff_dst)
    gff_to_hints(repeats_gff_src, hints_dst)

    return masked_dst, hints_dst
