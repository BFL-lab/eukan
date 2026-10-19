"""Functional annotation pipeline: homology search → FASTA + GFF3 annotation.

Doesn't fit ``run_simple_pipeline``: the FASTA/GFF3 annotation steps take
the search step's hits as arguments (read back from its JSON cache when the
search was skipped), so the steps aren't independent in the way the linear
driver assumes. Each step is driven directly via ``run_orchestrated_step``
with a per-proteome+mode scoped manifest key (see :func:`_step_scope`) and
explicit input fingerprints (see :func:`_func_step_fingerprints`).

The homology search has two modes controlled by ``config.homology_db``, one
``_MODES`` entry each:

* ``"uniprot"`` — phmmer of the proteome vs UniProt-SwissProt (default).
* ``"kofam"``   — hmmscan of the proteome vs the pressed KOfam HMM
  database, with per-KO score thresholds read from ``ko_list``.

Pfam hmmscan runs in both modes. The search caches the homology and Pfam
hits together in one file per proteome in the work dir, whose name changes per
mode (``phmmer.search.json`` vs ``kofam.search.json``) so that switching modes
doesn't quietly reuse stale results. The annotated FASTA goes to the work dir
too: nothing is written next to the input proteome.
"""

from __future__ import annotations

import gc
import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

from eukan.functional.kofam import run_kofam_search
from eukan.functional.search import (
    HitApplier,
    HitFormatter,
    HitResults,
    _load_digital_sequences,
    _open_hmm_db,
    annotate_fasta,
    annotate_gff3,
    apply_kofam_hits,
    apply_uniprot_hits,
    format_kofam_hit,
    format_uniprot_hit,
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
from eukan.infra.steps import is_step_complete, run_dir_relative
from eukan.settings import FunctionalConfig

log = get_logger(__name__)


@dataclass(frozen=True)
class _HomologyMode:
    """Everything that depends on ``config.homology_db``."""

    cache_suffix: str  # names the search cache and scopes the step keys
    search: Callable[[FunctionalConfig, float], HitResults]  # (config, E-value cutoff)
    db_scalars: Callable[[FunctionalConfig], list[str]]  # its databases, for the search fingerprint
    format_hit: HitFormatter  # one homology hit in a FASTA header
    apply_hits: HitApplier  # a query's homology hits onto its GFF3 feature


_MODES: dict[str, _HomologyMode] = {
    "uniprot": _HomologyMode(
        cache_suffix="phmmer",
        search=lambda c, evalue: _run_uniprot_phmmer(c.proteins, c.uniprot_db, c.num_cpu, evalue),
        db_scalars=lambda c: [f"uniprot_db={c.uniprot_db}"],
        format_hit=format_uniprot_hit,
        apply_hits=apply_uniprot_hits,
    ),
    "kofam": _HomologyMode(
        cache_suffix="kofam",
        search=lambda c, evalue: run_kofam_search(
            c.proteins, c.kofam_db, c.ko_list_path, c.num_cpu, evalue,
        ),
        db_scalars=lambda c: [f"kofam_db={c.kofam_db}", f"ko_list={c.ko_list_path}"],
        format_hit=format_kofam_hit,
        apply_hits=apply_kofam_hits,
    ),
}


def _search_cache_path(config: FunctionalConfig) -> Path:
    """Return the search step's JSON cache (homology and Pfam hits) for the active
    homology DB, in the work dir."""
    suffix = _MODES[config.homology_db].cache_suffix
    return config.work_dir / f"{config.proteins.stem}.{suffix}.search.json"


def _step_scope(config: FunctionalConfig) -> str:
    """Per-(proteome, mode) discriminator that namespaces manifest step keys.

    The search/annotate steps name their cache files and outputs after the
    proteome stem, and the homology-DB mode selects the cache suffix, so two
    different proteomes — or the same proteome run in two modes — write
    independent files into the shared work_dir. The manifest's
    skip-if-complete check, however, keys only on the step name: without a
    discriminator a second proteome inherits the first's "search complete"
    record, the search is skipped, and the pipeline then crashes reading a
    cache file that was never written for it.

    Scoping every functional step key by ``(stem, mode, path-hash)`` keeps
    resume correct for each proteome+mode while still skipping a genuine
    re-run of the same one. The path hash disambiguates two proteomes that
    share a stem but live in different directories. Those write the same file
    names into the work_dir, the later run replacing the earlier one's, so each
    needs its own record: going back to the first then finds its cache
    rewritten (its checksum no longer matches) and searches it again. It
    hashes the path as the manifest records it, relative to the run dir when
    the proteome lies inside it, so a moved run dir keeps its keys. Scoping handles
    *which* proteome+mode; the input fingerprints (see
    :func:`_func_step_fingerprints`) handle *content* changes within one scope.
    """
    suffix = _MODES[config.homology_db].cache_suffix
    path = run_dir_relative(config.proteins, config.manifest_dir)
    digest = hashlib.md5(str(path).encode()).hexdigest()[:8]
    return f"{config.proteins.stem}.{suffix}.{digest}"


def _run_uniprot_phmmer(
    proteins: Path, uniprot_db: Path, num_cpu: int, evalue: float,
) -> HitResults:
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
    log.info("Loading proteome from %s", proteins)
    queries = _load_digital_sequences(proteins)
    log.info("Loading Pfam HMMs from %s", pfam_db)
    with _open_hmm_db(pfam_db) as profiles:
        # A pressed database arrives as one block of profiles; hmmscan reads an
        # unpressed one itself, so its size is not known here.
        if isinstance(profiles, pyhmmer.plan7.OptimizedProfileBlock):
            target = f"{len(profiles)} profiles"
        else:
            target = "the unpressed HMMs"
        log.info("Running hmmscan (%d queries vs %s, %d CPUs)...", len(queries), target, num_cpu)
        return run_hmmscan_search(queries, profiles, num_cpu, evalue)


def _search_and_cache(config: FunctionalConfig, cache: Path) -> dict[str, HitResults]:
    """Run the active homology search + Pfam hmmscan; cache and return their hits.

    The hits go to *cache* as ``{"homology": ..., "pfam": ...}``. Stages run
    sequentially with the target database released between them — keeping
    all profiles + SwissProt resident at once was OOM-prone on container
    runtimes.
    """
    evalue_f = float(config.evalue)
    homology_res = _MODES[config.homology_db].search(config, evalue_f)
    pfam_res = _run_pfam_hmmscan(
        config.proteins, config.pfam_db, config.num_cpu, evalue_f,
    )

    results = {"homology": homology_res, "pfam": pfam_res}
    cache.write_text(json.dumps(results))
    return results


def _func_step_fingerprints(
    config: FunctionalConfig, cache: Path,
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
    on the search cache, so a re-run of ``search`` (which rewrites it)
    cascades into re-annotation. The cache's place in the run dir enters the
    search fingerprint (``cache_path=``), so a record from an earlier cache
    layout re-runs the search: the homology and Pfam hits in two files (no
    cache scalar), or the cache beside the proteome (fingerprinted by its
    name, as ``cache=``).
    """
    search_scalars = [
        f"homology_db={config.homology_db}",
        f"evalue={config.evalue}",
        f"pfam_db={config.pfam_db}",
        *_MODES[config.homology_db].db_scalars(config),
        f"cache_path={run_dir_relative(cache, config.manifest_dir)}",
    ]
    steps: dict[str, tuple[list[Path], list[str]]] = {
        "search": ([config.proteins], search_scalars),
        "annotate_fasta": (
            [config.proteins, cache],
            [f"homology_db={config.homology_db}"],
        ),
    }
    if config.gff3_path:
        steps["annotate_gff3"] = (
            [config.gff3_path, cache],
            [f"homology_db={config.homology_db}"],
        )
    return steps


def run_functional_annotation(
    config: FunctionalConfig, *, force: bool = False,
) -> None:
    """Run the full functional annotation pipeline."""
    config.ensure_dirs()
    work_dir = config.work_dir
    # The manifest lives in manifest_dir (the run-dir root), shared with the other
    # pipelines; step outputs and sentinels stay under work_dir.
    manifest_dir = config.manifest_dir
    manifest = get_or_create_manifest(manifest_dir, config)
    save_manifest(manifest_dir, manifest)

    mode = _MODES[config.homology_db]
    cache = _search_cache_path(config)

    # Namespace step keys per proteome+mode so several proteomes can be
    # annotated in one work_dir with independent manifest records (see
    # _step_scope, which also fixes a second-proteome crash). The step *dirs*
    # stay unscoped — they only hold a transient .running sentinel. Within a
    # scope, the input fingerprints below re-run a step on a content/mode
    # change and skip a byte-identical re-run.
    scope = _step_scope(config)

    def _key(name: str) -> str:
        return step_key(FUNCTIONAL, f"{name}.{scope}")

    fingerprints = _func_step_fingerprints(config, cache)

    # Identical input already fully annotated -> a friendly no-op. A different
    # proteins/GFF3 file (name or content) flips at least one fingerprint, so
    # the guard falls through to a normal run without needing -f.
    if not force and all(
        is_step_complete(
            manifest_dir, manifest, _key(name),
            input_files=ifiles, input_scalars=iscalars,
        ) is not None
        for name, (ifiles, iscalars) in fingerprints.items()
    ):
        fasta_out = work_dir / f"{config.proteins.stem}.mod{config.proteins.suffix}"
        log.info(
            "Already annotated %s -> %s; re-run with -f to recompute.",
            config.proteins.name, fasta_out.name,
        )
        return

    search_files, search_scalars = fingerprints["search"]
    computed: list[dict[str, HitResults]] = []
    run_orchestrated_step(
        manifest_dir, manifest, _key("search"),
        lambda: computed.append(_search_and_cache(config, cache)),
        step_dir=work_dir / "search",
        force=force,
        output_file=cache,
        input_files=search_files,
        input_scalars=search_scalars,
    )
    # A skipped search left its hits in the cache.
    results = computed[0] if computed else json.loads(cache.read_text())
    homology_res, hmmscan_res = results["homology"], results["pfam"]

    fasta_files, fasta_scalars = fingerprints["annotate_fasta"]
    run_orchestrated_step(
        manifest_dir, manifest, _key("annotate_fasta"),
        annotate_fasta, config.proteins, homology_res, hmmscan_res, mode.format_hit, work_dir,
        step_dir=work_dir / "annotate_fasta",
        force=force,
        input_files=fasta_files,
        input_scalars=fasta_scalars,
    )

    if config.gff3_path:
        gff3_files, gff3_scalars = fingerprints["annotate_gff3"]
        run_orchestrated_step(
            manifest_dir, manifest, _key("annotate_gff3"),
            annotate_gff3, config.gff3_path, homology_res, hmmscan_res,
            mode.apply_hits, work_dir,
            step_dir=work_dir / "annotate_gff3",
            force=force,
            input_files=gff3_files,
            input_scalars=gff3_scalars,
        )
