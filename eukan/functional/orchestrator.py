"""Functional annotation pipeline orchestration."""

from __future__ import annotations

from pathlib import Path

from eukan.functional.search import (
    annotate_fasta,
    annotate_gff3,
    run_homology_search,
)
from eukan.infra.logging import get_logger

log = get_logger(__name__)


def run_functional_annotation(
    proteins: Path,
    uniprot_db: Path | None = None,
    pfam_db: Path | None = None,
    gff3_path: Path | None = None,
    num_cpu: int = 1,
    evalue: str = "1e-1",
) -> None:
    """Run the full functional annotation pipeline.

    If uniprot_db or pfam_db are None, defaults are loaded from
    FunctionalConfig (pyproject.toml / EUKAN_FUNC_ env vars).
    """
    from eukan.settings import FunctionalConfig

    defaults = FunctionalConfig(proteins=proteins)
    uniprot_db = uniprot_db or defaults.uniprot_db
    pfam_db = pfam_db or defaults.pfam_db

    phmmer_res, hmmscan_res = run_homology_search(
        proteins, uniprot_db, pfam_db, num_cpu, evalue,
    )

    annotate_fasta(proteins, phmmer_res, hmmscan_res)

    if gff3_path:
        annotate_gff3(gff3_path, phmmer_res, hmmscan_res)
