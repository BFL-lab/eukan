"""Functional annotation pipeline orchestration."""

from __future__ import annotations

from pathlib import Path

from eukan.functional.search import (
    annotate_fasta,
    annotate_gff3,
    run_homology_search,
)
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    RunManifest, StepStatus, get_or_create_manifest,
    is_step_complete, pipeline_step, save_manifest,
    validate_step_outputs,
)

log = get_logger(__name__)


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
        errors = _validate_step_outputs(manifest)
        if errors:
            for msg in errors:
                log.error(msg)
            raise SystemExit(1)

    save_manifest(work_dir, manifest)

    # Step 1: Homology search (phmmer + hmmscan)
    # Results are in-memory dicts; serialize to JSON for caching/manifest.
    import json
    phmmer_json = proteins.parent / f"{proteins.stem}.phmmer.json"
    hmmscan_json = proteins.parent / f"{proteins.stem}.hmmscan.json"

    step_name = "functional/search"
    if force or not is_step_complete(manifest, step_name):
        with pipeline_step(work_dir, manifest, step_name) as step:
            phmmer_res, hmmscan_res = run_homology_search(
                proteins, uniprot_db, pfam_db, num_cpu, evalue,
            )
            phmmer_json.write_text(json.dumps(phmmer_res))
            hmmscan_json.write_text(json.dumps(hmmscan_res))
            step.output_file = str(phmmer_json)
    else:
        phmmer_res = json.loads(phmmer_json.read_text())
        hmmscan_res = json.loads(hmmscan_json.read_text())

    # Step 2: Annotate FASTA
    step_name = "functional/annotate_fasta"
    if force or not is_step_complete(manifest, step_name):
        with pipeline_step(work_dir, manifest, step_name) as step:
            output = annotate_fasta(proteins, phmmer_res, hmmscan_res)
            step.output_file = str(output)

    # Step 3: Annotate GFF3 (optional)
    if gff3_path:
        step_name = "functional/annotate_gff3"
        if force or not is_step_complete(manifest, step_name):
            with pipeline_step(work_dir, manifest, step_name) as step:
                output = annotate_gff3(gff3_path, phmmer_res, hmmscan_res)
                step.output_file = str(output)


_FUNCTIONAL_STEPS = [
    "functional/search", "functional/annotate_fasta", "functional/annotate_gff3",
]

_FUNCTIONAL_FLAGS = {s: "--force" for s in _FUNCTIONAL_STEPS}


def _validate_step_outputs(manifest: RunManifest) -> list[str]:
    """Validate that completed functional steps have non-empty output."""
    return validate_step_outputs(manifest, _FUNCTIONAL_STEPS, _FUNCTIONAL_FLAGS)
