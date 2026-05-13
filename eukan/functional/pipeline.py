"""Functional annotation pipeline: homology search → FASTA + GFF3 annotation.

Doesn't fit ``run_simple_pipeline``: the search step writes two JSON
caches that the FASTA/GFF3 annotation steps read back, so the steps
aren't independent in the way the linear driver assumes. StepSpec is
still used for the step declarations to keep the shape consistent with
the other pipelines.

The homology search has two modes controlled by ``config.homology_db``:

* ``"uniprot"`` — phmmer of the proteome vs UniProt-SwissProt (default).
* ``"kofam"``   — hmmscan of the proteome vs the pressed KOfam HMM
  database, with per-KO score thresholds read from ``ko_list``.

Pfam hmmscan runs in both modes. The cache filename used for homology
results changes per mode (``phmmer.json`` vs ``kofam.json``) so that
switching modes doesn't quietly reuse stale results.
"""

from __future__ import annotations

import gc
import json
from pathlib import Path

from eukan.functional.search import (
    HitResults,
    annotate_fasta,
    annotate_gff3,
    run_hmmscan_search,
    run_phmmer_search,
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


def _homology_cache_path(config: FunctionalConfig) -> Path:
    """Return the JSON cache path for the active homology DB."""
    stem = config.proteins.stem
    suffix = "kofam" if config.homology_db == "kofam" else "phmmer"
    return config.proteins.parent / f"{stem}.{suffix}.json"


def _run_uniprot_phmmer(
    proteins: Path, uniprot_db: Path, num_cpu: int, evalue: float,
) -> HitResults:
    from eukan.functional.search import _load_digital_sequences
    log.info("Loading proteome from %s", proteins)
    queries = _load_digital_sequences(proteins)
    log.info("Loading UniProt database from %s", uniprot_db)
    targets = _load_digital_sequences(uniprot_db)
    log.info(
        "Running phmmer (%d queries vs %d targets, %d CPUs)...",
        len(queries), len(targets), num_cpu,
    )
    res = run_phmmer_search(queries, targets, num_cpu, evalue)
    del targets
    gc.collect()
    return res


def _run_pfam_hmmscan(
    proteins: Path, pfam_db: Path, num_cpu: int, evalue: float,
) -> HitResults:
    from eukan.functional.search import _load_digital_sequences, _load_hmm_db
    log.info("Loading proteome from %s", proteins)
    queries = _load_digital_sequences(proteins)
    log.info("Loading Pfam HMMs from %s", pfam_db)
    hmms = _load_hmm_db(pfam_db)
    log.info(
        "Running hmmscan (%d queries vs %d profiles, %d CPUs)...",
        len(queries), len(hmms), num_cpu,
    )
    res = run_hmmscan_search(queries, hmms, num_cpu, evalue)
    del hmms
    gc.collect()
    return res


def _search_and_cache(
    config: FunctionalConfig,
    homology_json: Path,
    hmmscan_json: Path,
) -> Path:
    """Run the active homology search + Pfam hmmscan; write both caches.

    Stages run sequentially with the target database released between
    them — keeping all profiles + SwissProt resident at once was
    OOM-prone on container runtimes.
    """
    evalue_f = float(config.evalue)

    if config.homology_db == "kofam":
        from eukan.functional.kofam import run_kofam_search
        homology_res = run_kofam_search(
            config.proteins, config.kofam_db, config.ko_list_path,
            config.num_cpu, evalue_f,
        )
    else:
        homology_res = _run_uniprot_phmmer(
            config.proteins, config.uniprot_db, config.num_cpu, evalue_f,
        )

    pfam_res = _run_pfam_hmmscan(
        config.proteins, config.pfam_db, config.num_cpu, evalue_f,
    )

    homology_json.write_text(json.dumps(homology_res))
    hmmscan_json.write_text(json.dumps(pfam_res))
    return homology_json


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

    homology_json = _homology_cache_path(config)
    hmmscan_json = config.proteins.parent / f"{config.proteins.stem}.hmmscan.json"

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "search"),
        _search_and_cache,
        config, homology_json, hmmscan_json,
        step_dir=work_dir / "search",
        force=force,
    )

    homology_res = json.loads(homology_json.read_text())
    hmmscan_res = json.loads(hmmscan_json.read_text())

    run_orchestrated_step(
        work_dir, manifest, step_key(FUNCTIONAL, "annotate_fasta"),
        annotate_fasta, config.proteins, homology_res, hmmscan_res,
        config.homology_db,
        step_dir=work_dir / "annotate_fasta",
        force=force,
    )

    if config.gff3_path:
        run_orchestrated_step(
            work_dir, manifest, step_key(FUNCTIONAL, "annotate_gff3"),
            annotate_gff3, config.gff3_path, homology_res, hmmscan_res,
            work_dir, config.homology_db,
            step_dir=work_dir / "annotate_gff3",
            force=force,
        )
