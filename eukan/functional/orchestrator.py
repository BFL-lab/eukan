"""Functional annotation pipeline orchestration."""

from __future__ import annotations

import json
from pathlib import Path

from eukan.functional.search import (
    annotate_fasta,
    annotate_gff3,
    run_homology_search,
)
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    FUNCTIONAL,
    get_or_create_manifest,
    run_orchestrated_step,
    save_manifest,
    step_key,
    validate_or_raise,
)
from eukan.settings import FunctionalConfig

log = get_logger(__name__)


_FUNCTIONAL_STEPS = [
    step_key(FUNCTIONAL, "search"),
    step_key(FUNCTIONAL, "annotate_fasta"),
    step_key(FUNCTIONAL, "annotate_gff3"),
]

_FUNCTIONAL_FLAGS = {s: "--force" for s in _FUNCTIONAL_STEPS}


def _search_and_cache(
    proteins: Path, uniprot_db: Path, pfam_db: Path, num_cpu: int, evalue: str,
    phmmer_json: Path, hmmscan_json: Path,
) -> Path:
    phmmer_res, hmmscan_res = run_homology_search(
        proteins, uniprot_db, pfam_db, num_cpu, evalue,
    )
    phmmer_json.write_text(json.dumps(phmmer_res))
    hmmscan_json.write_text(json.dumps(hmmscan_res))
    return phmmer_json


def run_functional_annotation(
    config: FunctionalConfig, *, force: bool = False,
) -> None:
    """Run the full functional annotation pipeline."""
    work_dir = config.work_dir
    manifest = get_or_create_manifest(work_dir, config)

    if not force:
        validate_or_raise(manifest, _FUNCTIONAL_STEPS, _FUNCTIONAL_FLAGS)

    save_manifest(work_dir, manifest)

    phmmer_json = config.proteins.parent / f"{config.proteins.stem}.phmmer.json"
    hmmscan_json = config.proteins.parent / f"{config.proteins.stem}.hmmscan.json"

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "search"),
        _search_and_cache,
        config.proteins, config.uniprot_db, config.pfam_db,
        config.num_cpu, config.evalue, phmmer_json, hmmscan_json,
        step_dir=work_dir / "search",
        force=force,
    )

    phmmer_res = json.loads(phmmer_json.read_text())
    hmmscan_res = json.loads(hmmscan_json.read_text())

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "annotate_fasta"),
        annotate_fasta, config.proteins, phmmer_res, hmmscan_res,
        step_dir=work_dir / "annotate_fasta",
        force=force,
    )

    if config.gff3_path:
        run_orchestrated_step(
            work_dir, manifest, step_key(FUNCTIONAL, "annotate_gff3"),
            annotate_gff3, config.gff3_path, phmmer_res, hmmscan_res,
            step_dir=work_dir / "annotate_gff3",
            force=force,
        )
