"""Assembly pipeline: STAR mapping → Trinity assembly → PASA alignment."""

from __future__ import annotations

from eukan.assembly.pasa import run_pasa
from eukan.assembly.star import map_reads
from eukan.assembly.trinity import run_trinity
from eukan.infra.artifacts import Artifact
from eukan.infra.manifest import ASSEMBLY
from eukan.infra.pipeline import (
    StepSpec,
    force_steps_from_run_flags as _force_steps_from_run_flags,
    run_simple_pipeline,
)
from eukan.settings import AssemblyConfig

_STEPS: list[StepSpec] = [
    StepSpec("star",    map_reads,   "STAR_Aligned.sortedByCoord.out.bam",  "-A / --run-star"),
    StepSpec("trinity", run_trinity, "trinity-gg.fasta",                    "-T / --run-trinity"),
    StepSpec("pasa",    run_pasa,    Artifact.NR_TRANSCRIPTS_FASTA.value,   "-P / --run-pasa"),
]


def force_steps_from_run_flags(
    *,
    run_star: bool = False,
    run_trinity: bool = False,
    run_pasa: bool = False,
    force: bool = False,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force."""
    return _force_steps_from_run_flags(
        ASSEMBLY, _STEPS,
        force=force,
        run_star=run_star, run_trinity=run_trinity, run_pasa=run_pasa,
    )


def run_assembly(
    config: AssemblyConfig,
    *,
    force_steps: list[str] | None = None,
) -> None:
    """Run the assembly pipeline with manifest tracking."""
    run_simple_pipeline(ASSEMBLY, _STEPS, config, force_steps=force_steps)
