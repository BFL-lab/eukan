"""De novo repeat-family inference via BuildDatabase + RepeatModeler.

The genome is sorted by sequence length (descending) and uppercased
before being indexed — RepeatModeler's sampler benefits from large
contigs appearing first, and the uppercase pass strips any prior
softmasking that would otherwise be invisible to the modeler.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import RepeatsConfig

log = get_logger(__name__)


def sort_and_uppercase(genome: Path, out_path: Path) -> None:
    """Write a length-sorted, uppercased copy of *genome* to *out_path*.

    Two-pass: first pass indexes the FASTA (no full load), second pass
    streams records out in length-descending order. Keeps memory bounded
    on multi-GB plant genomes.
    """
    index = SeqIO.index(str(genome), "fasta")
    try:
        order = sorted(index, key=lambda rid: len(index[rid]), reverse=True)
        with open(out_path, "w") as fh:
            for rid in order:
                rec = index[rid]
                rec.seq = rec.seq.upper()
                # description=id avoids duplicate header tokens from Biopython
                rec.description = ""
                SeqIO.write([rec], fh, "fasta")
    finally:
        index.close()


def run_modeler(config: RepeatsConfig) -> Path:
    """Build an rmblast database and infer repeat families.

    Returns the path to the families FASTA produced by RepeatModeler.
    """
    sdir = config.work_dir / "repeats" / "modeler"
    sdir.mkdir(parents=True, exist_ok=True)

    sorted_fa = sdir / f"{config.name}.sorted.fasta"
    db_name = f"{config.name}.replib"

    log.info("Sorting and uppercasing genome for RepeatModeler input...")
    sort_and_uppercase(config.genome, sorted_fa)

    log.info("Building rmblast database...")
    # RepeatModeler 2.x's BuildDatabase only supports rmblast and dropped the
    # -engine flag entirely; passing it raises "Unknown option: engine".
    run_cmd(
        ["BuildDatabase", "-name", db_name, sorted_fa.name],
        cwd=sdir,
    )

    log.info("Running RepeatModeler (this may take many hours)...")
    run_cmd(
        ["RepeatModeler", "-database", db_name, "-threads", str(config.num_cpu)],
        cwd=sdir,
    )

    families = sdir / f"{db_name}-families.fa"
    if not families.exists():
        # RepeatModeler 2.x default; fall back to legacy locations if needed.
        candidates = list(sdir.glob(f"{db_name}-families.fa")) + list(
            sdir.glob("RM_*/consensi.fa.classified")
        )
        if candidates:
            families = candidates[0]
    return families
