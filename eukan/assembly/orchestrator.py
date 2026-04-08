"""Assembly pipeline orchestration."""

from __future__ import annotations

from eukan.assembly.pasa import run_pasa
from eukan.assembly.star import map_reads
from eukan.assembly.trinity import run_trinity
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def run_assembly(config: AssemblyConfig, steps: list[str], force: bool = False) -> None:
    """Run the specified assembly steps."""
    if "map" in steps:
        map_reads(config, force=force)
    if "trinity" in steps:
        run_trinity(config, force=force)
    if "pasa" in steps:
        run_pasa(config, force=force)
