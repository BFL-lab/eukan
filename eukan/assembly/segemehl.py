"""segemehl read mapping — splice-agnostic alternative to STAR.

Unlike STAR, segemehl does not enforce canonical GT-AG splice sites, so it
captures non-canonical introns (e.g. the dominant CG-AG introns of diplonemids
such as *Hemistasia*) that STAR would miss or misplace. segemehl has no native
splice-junction table, so we derive a STAR-format ``SJ.out.tab`` from the BAM's
N-CIGAR junctions (:func:`align_hints.sj_table_from_bam`) and reuse the shared
post-alignment processing verbatim. Downstream steps (GeneMark, AUGUSTUS, and
transcript assembly) therefore see the identical contract STAR produces —
including the ``splice_site_summary.json`` that lets AUGUSTUS allow the
non-canonical splice sites via ``--allow_hinted_splicesites``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from eukan.assembly.align_hints import generate_rnaseq_hints, sj_table_from_bam
from eukan.exceptions import ExternalToolError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_piped
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_BAM = "segemehl_Aligned.sortedByCoord.out.bam"
_UNSORTED_BAM = "segemehl_unsorted.bam"
_SJ = "segemehl_SJ.out.tab"
_INDEX = "segemehl.idx"
# segemehl emits split-read byproducts as <base>.sngl.bed / .mult.bed / .trns.txt;
# pin <base> to a known path so they can be cleaned up afterwards.
_SPLITS_BASE = "segemehl_splits"
_SPLITS_SUFFIXES = (".sngl.bed", ".mult.bed", ".trns.txt")

# The coordinate-sort needs transient room for its temp spill (~= the sorted
# size) plus the final BAM, both landing on the run-dir filesystem — roughly
# twice the unsorted BAM. Warn before sorting when free space is below this.
_SORT_DISK_HEADROOM = 2
# samtools' multithreaded BGZF writer reports a failed write on a full disk
# with a misleading errno (e.g. "...: Illegal seek"); match these stderr
# fragments so we can translate the cryptic failure into a clear out-of-space
# message instead of letting the user chase a phantom seek bug.
_DISK_FULL_MARKERS = ("illegal seek", "no space", "failed writing", "write failed")


def _bam_is_complete(path: Path) -> bool:
    """True if *path* is a non-empty BAM with a valid BGZF EOF (``quickcheck``).

    Lets a re-run reuse an unsorted BAM from a previous attempt that failed
    only in the downstream sort (segemehl mapping is the multi-hour step, so we
    never want to redo it when its output survived). A truncated BAM from an
    interrupted or killed segemehl fails quickcheck and is re-mapped.
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        run_cmd(["samtools", "quickcheck", str(path)], cwd=path.parent)
    except ExternalToolError:
        return False
    return True


def _write_unmapped_fasta(unsorted_bam: Path, out_fasta: Path) -> int:
    """Extract unmapped records from *unsorted_bam* to *out_fasta* (FASTA).

    The segemehl transcript BAM keeps unmapped queries before the ``-F 4`` filter
    drops them; pulling them out here (pysam, no samtools dependency) preserves the
    transcripts that failed to map for inspection, matching what the STARlong path
    saves via ``--outReadsUnmapped``. Returns the number written.
    """
    import pysam

    n = 0
    bam = pysam.AlignmentFile(str(unsorted_bam), "rb")
    try:
        with open(out_fasta, "w") as fh:
            for read in bam:
                if not read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                seq = read.query_sequence
                if not read.query_name or not seq:
                    continue
                fh.write(f">{read.query_name}\n{seq}\n")
                n += 1
    finally:
        bam.close()
    return n


def _coordinate_sort_and_filter(
    unsorted: Path, out_bam: str, wd: Path, num_cpu: int
) -> None:
    """Drop unmapped reads (``-F 4``) and coordinate-sort into *out_bam*.

    Mirrors the validated recipe's ``samtools sort | samtools view -bh -F 4``;
    filtering before the sort is cheaper and yields the same result. On a tight
    disk this is the step that fails — samtools needs roughly twice the unsorted
    BAM in transient temp + output space — so we warn up front and, on a
    write/seek failure, re-raise with an explicit out-of-space hint.
    """
    if unsorted.exists():
        src = unsorted.stat().st_size
        free = shutil.disk_usage(wd).free
        if free < _SORT_DISK_HEADROOM * src:
            log.warning(
                "Low disk for coordinate-sort: %.1f GB free vs ~%.1f GB likely "
                "needed (temp spill + output for a %.1f GB BAM). If it fails, "
                "free space and re-run — the unsorted BAM is reused, so segemehl "
                "will not re-map.",
                free / 1e9, _SORT_DISK_HEADROOM * src / 1e9, src / 1e9,
            )
    try:
        run_piped(
            ["samtools", "view", "-bh", "-F", "4", str(unsorted)],
            ["samtools", "sort", "-@", str(num_cpu), "-o", out_bam, "-"],
            cwd=wd,
        )
    except ExternalToolError as exc:
        snippet = (exc.stderr_snippet or "").lower()
        if any(marker in snippet for marker in _DISK_FULL_MARKERS):
            raise ExternalToolError(
                "samtools sort failed writing the coordinate-sorted BAM, which "
                "on a full filesystem samtools reports as 'Illegal seek'. The "
                "sort needs roughly twice the unsorted BAM in transient temp + "
                "output space.",
                tool=exc.tool, returncode=exc.returncode, cmd=exc.cmd,
                stderr_snippet=exc.stderr_snippet,
                hint=(
                    "Free space on the run-dir filesystem (or run the assemble "
                    "step on a larger disk), then re-run: segemehl reuses the "
                    "existing unsorted BAM and skips the slow re-mapping."
                ),
            ) from exc
        raise


def map_reads_segemehl(config: AssemblyConfig) -> None:
    """Map RNA-seq reads to the genome using segemehl in split/spliced mode."""
    wd = config.work_dir
    log.info("Running segemehl read mapping...")

    index = wd / _INDEX
    unsorted = wd / _UNSORTED_BAM

    # segemehl mapping is the multi-hour step. If a previous run already
    # produced a complete unsorted BAM (and only the downstream sort failed —
    # e.g. on a full disk), reuse it rather than re-mapping. quickcheck rejects
    # a truncated BAM from an interrupted segemehl, which is then re-mapped.
    if _bam_is_complete(unsorted):
        log.info(
            "Reusing %s from a previous run; skipping segemehl mapping "
            "(delete it to force a full re-map).",
            unsorted.name,
        )
    else:
        if not index.exists():
            run_cmd(["segemehl.x", "-x", str(index), "-d", str(config.genome)], cwd=wd)

        # segemehl loads the whole suffix-array index + genome into memory, so
        # it is RAM-hungry on large genomes. Write the BAM straight to disk
        # with `-b -o` instead of piping into samtools: run_cmd then checks
        # segemehl's *own* exit status, so an OOM kill, a mate-count mismatch,
        # or any other crash surfaces as a clear segemehl error rather than a
        # misleading downstream "samtools: failed to read header".
        #
        # Flags mirror the validated `segemehl-H0.bam` recipe: `-H 0` reports
        # all (not just best-scoring) alignments, capturing multi-mapping
        # spliced reads; `-S <base>` enables split/spliced mapping and pins the
        # byproduct BED files to a known location for cleanup. segemehl reads
        # gzipped FASTQ natively, so `.fastq.gz` inputs need no decompression.
        run_cmd(
            [
                "segemehl.x",
                "-H", "0",
                "-i", str(index),
                "-d", str(config.genome),
                *config.reads_args_segemehl,
                "-S", str(wd / _SPLITS_BASE),
                "-t", str(config.num_cpu),
                "-b",
                "-o", str(unsorted),
            ],
            cwd=wd,
        )

    # Mapping is done; the index and split-read BEDs are never read again
    # (junctions come from the BAM). Delete them now, *before* the sort, so
    # their several GB don't compete with the sort's temp spill + output on a
    # tight disk — a `-H 0` run plus sort scratch can need well over 10 GB.
    index.unlink(missing_ok=True)
    for suffix in _SPLITS_SUFFIXES:
        (wd / f"{_SPLITS_BASE}{suffix}").unlink(missing_ok=True)

    _coordinate_sort_and_filter(unsorted, _BAM, wd, config.num_cpu)
    run_cmd(["samtools", "index", _BAM], cwd=wd)

    # Derive a STAR-format SJ.out.tab from the BAM, then reuse STAR's
    # post-processing so the downstream hints / splice summary are identical.
    sj = sj_table_from_bam(
        wd / _BAM, config.genome, wd,
        min_intron=config.min_intron_len,
        max_intron=config.max_intron_len,
        out_name=_SJ,
    )
    generate_rnaseq_hints(
        sj, wd / _BAM, config.genome, wd,
        diagnose=config.diagnose_softclips, source_label="segemehl",
    )

    # The unsorted BAM is no longer needed once the sorted BAM + hints exist.
    unsorted.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Transcript -> genome mapping config (consumed by star.map_transcripts)
# ---------------------------------------------------------------------------
# Trinity transcripts are mapped to the genome SPLICED (bounded intron + Local
# soft-clip). Both Trinity modes emit transcript-coordinate FASTAs (genome-guided
# only clusters reads per locus), so both are mapped here. The dispatcher
# (eukan.assembly.star.map_transcripts) uses STARlong by default, with segemehl
# `-S` (`-H 1`, see map_one_transcript_set_segemehl) as the fallback when STARlong
# fails/under-maps — and as the *primary* mapper when non-canonical splicing is
# extensive (or --aligner segemehl). These (query FASTA, output BAM) pairs and the
# jaccard-sibling resolver are shared across both paths and are the single source
# of truth for the mapped-transcript stems (see assembly.tracks).
_TRANSCRIPT_SETS: tuple[tuple[str, str], ...] = (
    ("trinity-denovo.fasta", "trinity-denovo.genome.bam"),
    ("trinity-gg.fasta", "trinity-gg.genome.bam"),
)
_GENOME_BAM_SUFFIX = ".genome.bam"


def _resolve_query(wd: Path, query_name: str) -> Path:
    """The transcript FASTA to map: the jaccard-clipped sibling if it exists.

    The jaccard step (``assembly/jaccard.py``) rewrites each de novo / genome-
    guided FASTA into a ``<stem>.jaccard.fasta``; prefer it so fused contigs are
    split before consolidation. Falls back to the original when clipping is off.
    """
    clipped = wd / query_name.replace(".fasta", ".jaccard.fasta")
    return clipped if clipped.exists() and clipped.stat().st_size > 0 else wd / query_name


_TX_INDEX = "segemehl_tx.idx"
_TX_SPLITS_BASE = "segemehl_tx_splits"


def map_one_transcript_set_segemehl(config: AssemblyConfig, query: Path, out_bam: str) -> None:
    """Spliced fallback for transcript->genome mapping (used when STARlong fails).

    segemehl ``-S`` natively splits long transcripts at introns, recovering the
    splice structure STAR can miss on long queries. ``-H 1`` (report only the best
    alignment) bounds the split-DP memory that ``-H 0`` blew past on this box; we
    extract the unmapped transcripts (for inspection — mirrors the STARlong primary
    path), then drop unmapped reads and coordinate-sort/index like the read path.
    """
    wd = config.work_dir
    final = wd / out_bam
    if _bam_is_complete(final):
        log.info("Reusing %s; skipping segemehl transcript mapping.", final.name)
        return

    index = wd / _TX_INDEX
    if not index.exists():
        run_cmd(["segemehl.x", "-x", str(index), "-d", str(config.genome)], cwd=wd)
    unsorted = wd / f"{out_bam}.unsorted.bam"
    log.warning(
        "Mapping %s with segemehl -S (-H 1) — memory-hungry; watch RSS on a "
        "low-memory box.", query.name,
    )
    run_cmd(
        [
            "segemehl.x",
            "-H", "1",
            "-i", str(index),
            "-d", str(config.genome),
            "-q", str(query),
            "-S", str(wd / _TX_SPLITS_BASE),
            "-t", str(config.num_cpu),
            "-b",
            "-o", str(unsorted),
        ],
        cwd=wd,
    )
    index.unlink(missing_ok=True)
    for suffix in _SPLITS_SUFFIXES:
        (wd / f"{_TX_SPLITS_BASE}{suffix}").unlink(missing_ok=True)
    stem = out_bam[: -len(_GENOME_BAM_SUFFIX)]
    _write_unmapped_fasta(unsorted, wd / f"{stem}.unmapped_transcripts.fasta")
    _coordinate_sort_and_filter(unsorted, out_bam, wd, config.num_cpu)
    run_cmd(["samtools", "index", out_bam], cwd=wd)
    unsorted.unlink(missing_ok=True)
