"""De novo repeat-family inference via BuildDatabase + RepeatModeler.

The genome is sorted by sequence length (descending) and uppercased
before being indexed — RepeatModeler's sampler benefits from large
contigs appearing first, and the uppercase pass strips any prior
softmasking that would otherwise be invisible to the modeler.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from Bio import SeqIO

from eukan.exceptions import ExternalToolError
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


def _salvage_partial_library(sdir: Path) -> Path | None:
    """Find the most-complete partial library after a RepeatModeler crash.

    RepeatModeler runs in numbered rounds; round-1 typically finishes
    even on small genomes that crash later in RECON. We probe in order
    of "most complete":

    1. ``RM_*/consensi.fa`` -- cumulative library appended after each
       successful round, classified up to the last completed round.
    2. ``RM_*/round-N/consensi.fa.classified`` -- per-round classified
       library (latest round first).
    3. ``RM_*/round-N/consensi.fa`` -- per-round unclassified consensus.
    4. ``RM_*/round-N/refined-cons.fa`` -- earliest viable form,
       refiner output before classification.

    Returns the first non-empty candidate, or ``None`` if nothing
    salvageable exists.
    """
    rm_dirs = sorted(sdir.glob("RM_*"))
    if not rm_dirs:
        return None
    rm = rm_dirs[-1]

    candidates: list[Path] = [rm / "consensi.fa"]
    for round_dir in sorted(rm.glob("round-*"), reverse=True):
        candidates += [
            round_dir / "consensi.fa.classified",
            round_dir / "consensi.fa",
            round_dir / "refined-cons.fa",
        ]
    for c in candidates:
        if c.exists() and c.stat().st_size > 0:
            return c
    return None


def run_modeler(config: RepeatsConfig) -> Path:
    """Build an rmblast database and infer repeat families.

    Returns the path to the families FASTA produced by RepeatModeler.
    On RepeatModeler failure, attempts to salvage a partial library from
    completed rounds so downstream RepeatMasker can still softmask.
    """
    sdir = config.work_dir / "modeler"
    sdir.mkdir(parents=True, exist_ok=True)

    sorted_fa = sdir / f"{config.name}.sorted.fasta"
    db_name = f"{config.name}.replib"
    families = sdir / f"{db_name}-families.fa"

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
    try:
        run_cmd(
            ["RepeatModeler", "-database", db_name, "-threads", str(config.num_cpu)],
            cwd=sdir,
        )
    except ExternalToolError as exc:
        partial = _salvage_partial_library(sdir)
        if partial is None:
            raise
        log.warning(
            "RepeatModeler failed (%s). Salvaged partial library from %s; "
            "softmasking will use whatever families completed before the crash.",
            exc, partial.relative_to(sdir),
        )
        shutil.copy2(partial, families)
        return families

    if not families.exists():
        # RepeatModeler 2.x default landing spot; older builds occasionally
        # leave the families library nested under RM_*/ instead.
        for candidate in sdir.glob("RM_*/consensi.fa.classified"):
            if candidate.exists():
                families = candidate
                break
    return families
