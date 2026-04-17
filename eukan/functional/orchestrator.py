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
    FUNCTIONAL, RunManifest, get_or_create_manifest, is_step_complete,
    run_orchestrated_step, save_manifest, step_key, validate_step_outputs,
)

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
    proteins: Path,
    uniprot_db: Path | None = None,
    pfam_db: Path | None = None,
    gff3_path: Path | None = None,
    num_cpu: int = 1,
    evalue: str = "1e-1",
    force: bool = False,
    manifest_dir: Path | None = None,
) -> None:
    """Run the full functional annotation pipeline.

    If uniprot_db or pfam_db are None, defaults are loaded from
    FunctionalConfig (pyproject.toml / EUKAN_FUNC_ env vars).

    Args:
        manifest_dir: Directory for eukan-run.json. Defaults to cwd.
    """
    from eukan.settings import FunctionalConfig

    defaults = FunctionalConfig(proteins=proteins)
    uniprot_db = uniprot_db or defaults.uniprot_db
    pfam_db = pfam_db or defaults.pfam_db

    work_dir = manifest_dir or Path.cwd()
    manifest = get_or_create_manifest(work_dir)

    if not force:
        errors = validate_step_outputs(manifest, _FUNCTIONAL_STEPS, _FUNCTIONAL_FLAGS)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

    save_manifest(work_dir, manifest)

    phmmer_json = proteins.parent / f"{proteins.stem}.phmmer.json"
    hmmscan_json = proteins.parent / f"{proteins.stem}.hmmscan.json"

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "search"),
        _search_and_cache,
        proteins, uniprot_db, pfam_db, num_cpu, evalue, phmmer_json, hmmscan_json,
        step_dir=work_dir / "search",
        force=force,
    )

    phmmer_res = json.loads(phmmer_json.read_text())
    hmmscan_res = json.loads(hmmscan_json.read_text())

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "annotate_fasta"),
        annotate_fasta, proteins, phmmer_res, hmmscan_res,
        step_dir=work_dir / "annotate_fasta",
        force=force,
    )

    if gff3_path:
        run_orchestrated_step(
            work_dir, manifest, step_key(FUNCTIONAL, "annotate_gff3"),
            annotate_gff3, gff3_path, phmmer_res, hmmscan_res,
            step_dir=work_dir / "annotate_gff3",
            force=force,
        )
