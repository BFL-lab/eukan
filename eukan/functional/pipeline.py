"""Functional annotation pipeline: phmmer search → FASTA + GFF3 annotation.

Doesn't fit ``run_simple_pipeline``: the search step writes two JSON
caches that the FASTA/GFF3 annotation steps read back, so the steps
aren't independent in the way the linear driver assumes. StepSpec is
still used for the step declarations to keep the shape consistent with
the other pipelines.
"""

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
    save_manifest,
    step_key,
)
from eukan.infra.pipeline import StepSpec, run_orchestrated_step
from eukan.infra.steps import validate_or_raise
from eukan.settings import FunctionalConfig

log = get_logger(__name__)


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


# Step specs used for the validate_or_raise stale-output check. fn fields
# are the inner search/annotate functions; the actual call sites below
# wrap them with the JSON-cache plumbing that doesn't fit StepSpec's
# (config) → output contract.
_STEPS: list[StepSpec] = [
    StepSpec("search",         _search_and_cache, flag="--force"),
    StepSpec("annotate_fasta", annotate_fasta,    flag="--force"),
    StepSpec("annotate_gff3",  annotate_gff3,     flag="--force"),
]


def run_functional_annotation(
    config: FunctionalConfig, *, force: bool = False,
) -> None:
    """Run the full functional annotation pipeline."""
    work_dir = config.work_dir
    manifest = get_or_create_manifest(work_dir, config)

    expected = [step_key(FUNCTIONAL, s.name) for s in _STEPS]
    flag_map = {key: "--force" for key in expected}
    if not force:
        validate_or_raise(manifest, expected, flag_map)

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
