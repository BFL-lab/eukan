"""Functional annotation pipeline: homology search → FASTA + GFF3 annotation.

Doesn't fit ``run_simple_pipeline``: the search step writes two JSON
caches that the FASTA/GFF3 annotation steps read back, so the steps
aren't independent in the way the linear driver assumes. Each step is
driven directly via ``run_orchestrated_step`` with a per-proteome+mode
scoped manifest key (see :func:`_step_scope`) and explicit input
fingerprints (see :func:`_func_step_fingerprints`).

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
import hashlib
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
from eukan.infra.pipeline import run_orchestrated_step
from eukan.infra.steps import is_step_complete
from eukan.settings import FunctionalConfig

log = get_logger(__name__)


def _homology_cache_path(config: FunctionalConfig) -> Path:
    """Return the JSON cache path for the active homology DB."""
    stem = config.proteins.stem
    suffix = "kofam" if config.homology_db == "kofam" else "phmmer"
    return config.proteins.parent / f"{stem}.{suffix}.json"


def _step_scope(config: FunctionalConfig) -> str:
    """Per-(proteome, mode) discriminator that namespaces manifest step keys.

    The search/annotate steps name their cache files and outputs after the
    proteome stem, and the homology-DB mode selects the cache suffix, so two
    different proteomes — or the same proteome run in two modes — already
    write independent files into a shared work_dir. The manifest's
    skip-if-complete check, however, keys only on the step name: without a
    discriminator a second proteome inherits the first's "search complete"
    record, the search is skipped, and the pipeline then crashes reading a
    cache file that was never written for it.

    Scoping every functional step key by ``(stem, mode, path-hash)`` keeps
    resume correct for each proteome+mode while still skipping a genuine
    re-run of the same one. The path hash disambiguates two proteomes that
    share a stem but live in different directories (they have distinct cache
    files, so they must not share a manifest record either). Scoping handles
    *which* proteome+mode; the input fingerprints (see
    :func:`_func_step_fingerprints`) handle *content* changes within one scope.
    """
    suffix = "kofam" if config.homology_db == "kofam" else "phmmer"
    digest = hashlib.md5(str(config.proteins.resolve()).encode()).hexdigest()[:8]
    return f"{config.proteins.stem}.{suffix}.{digest}"


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


def _func_step_fingerprints(
    config: FunctionalConfig,
    homology_json: Path,
    hmmscan_json: Path,
) -> dict[str, tuple[list[Path], list[str]]]:
    """Declared ``(input_files, input_scalars)`` per step, keyed by bare name.

    The single source of truth shared by the resume guard and the
    ``run_orchestrated_step`` calls below, so their fingerprints can't drift.
    Step keys are already scoped per proteome+mode (see :func:`_step_scope`),
    so this layer's job is to detect *content* changes within a scope:
    ``fingerprint_inputs`` keys each file by ``path=md5``, so the same proteome
    edited in place (same path, new bytes) flips the digest and re-runs the
    step, while a byte-identical input reuses the cached result. Databases
    enter as path-string scalars, not content-hashed: they're multi-GB and
    re-hashing them every run would be prohibitive. ``num_cpu`` is deliberately
    excluded — thread count doesn't change results. The annotate steps depend
    on both JSON caches, so a re-run of ``search`` (which rewrites them)
    cascades into re-annotation; this also gives the otherwise-untracked
    ``hmmscan.json`` an integrity-relevant fingerprint.
    """
    if config.homology_db == "kofam":
        db_scalars = [
            f"kofam_db={config.kofam_db}",
            f"ko_list={config.ko_list_path}",
        ]
    else:
        db_scalars = [f"uniprot_db={config.uniprot_db}"]
    search_scalars = [
        f"homology_db={config.homology_db}",
        f"evalue={config.evalue}",
        f"pfam_db={config.pfam_db}",
        *db_scalars,
    ]
    steps: dict[str, tuple[list[Path], list[str]]] = {
        "search": ([config.proteins], search_scalars),
        "annotate_fasta": (
            [config.proteins, homology_json, hmmscan_json],
            [f"homology_db={config.homology_db}"],
        ),
    }
    if config.gff3_path:
        steps["annotate_gff3"] = (
            [config.gff3_path, homology_json, hmmscan_json],
            [f"homology_db={config.homology_db}"],
        )
    return steps


def run_functional_annotation(
    config: FunctionalConfig, *, force: bool = False,
) -> None:
    """Run the full functional annotation pipeline."""
    work_dir = config.work_dir
    manifest = get_or_create_manifest(work_dir, config)
    save_manifest(work_dir, manifest)

    homology_json = _homology_cache_path(config)
    hmmscan_json = config.proteins.parent / f"{config.proteins.stem}.hmmscan.json"

    # Namespace step keys per proteome+mode so several proteomes can be
    # annotated in one work_dir with independent manifest records (see
    # _step_scope, which also fixes a second-proteome crash). The step *dirs*
    # stay unscoped — they only hold a transient .running sentinel. Within a
    # scope, the input fingerprints below re-run a step on a content/mode
    # change and skip a byte-identical re-run.
    scope = _step_scope(config)

    def _key(name: str) -> str:
        return step_key(FUNCTIONAL, f"{name}.{scope}")

    fingerprints = _func_step_fingerprints(config, homology_json, hmmscan_json)

    # Identical input already fully annotated -> a friendly no-op. A different
    # proteins/GFF3 file (name or content) flips at least one fingerprint, so
    # the guard falls through to a normal run without needing -f.
    if not force and all(
        is_step_complete(
            manifest, _key(name),
            input_files=ifiles, input_scalars=iscalars,
        ) is not None
        for name, (ifiles, iscalars) in fingerprints.items()
    ):
        fasta_out = config.proteins.parent / f"{config.proteins.stem}.mod{config.proteins.suffix}"
        log.info(
            "Already annotated %s -> %s; re-run with -f to recompute.",
            config.proteins.name, fasta_out.name,
        )
        return

    search_files, search_scalars = fingerprints["search"]
    run_orchestrated_step(
        work_dir, manifest, _key("search"),
        _search_and_cache,
        config, homology_json, hmmscan_json,
        step_dir=work_dir / "search",
        force=force,
        input_files=search_files,
        input_scalars=search_scalars,
    )

    homology_res = json.loads(homology_json.read_text())
    hmmscan_res = json.loads(hmmscan_json.read_text())

    fasta_files, fasta_scalars = fingerprints["annotate_fasta"]
    run_orchestrated_step(
        work_dir, manifest, _key("annotate_fasta"),
        annotate_fasta, config.proteins, homology_res, hmmscan_res,
        config.homology_db,
        step_dir=work_dir / "annotate_fasta",
        force=force,
        input_files=fasta_files,
        input_scalars=fasta_scalars,
    )

    if config.gff3_path:
        gff3_files, gff3_scalars = fingerprints["annotate_gff3"]
        run_orchestrated_step(
            work_dir, manifest, _key("annotate_gff3"),
            annotate_gff3, config.gff3_path, homology_res, hmmscan_res,
            work_dir, config.homology_db,
            step_dir=work_dir / "annotate_gff3",
            force=force,
            input_files=gff3_files,
            input_scalars=gff3_scalars,
        )
