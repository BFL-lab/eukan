"""STAR read mapping and hint generation."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

from eukan.assembly.align_hints import generate_rnaseq_hints
from eukan.exceptions import ExternalToolError
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# Non-canonical-splice verdict (softclip_diagnostic_summary.json →
# verdict.non_canonical_splice.call) at which transcript→genome mapping switches
# from STARlong to splice-agnostic segemehl: STARlong's canonical-tuned seeding
# under-maps transcripts whose introns are mostly non-canonical.
_SEGEMEHL_NON_CANONICAL_CALL = "EXTENSIVE"

# STAR resource/index tuning for typical small eukaryotic genomes.
_STAR_SA_INDEX_NBASES = "3"                 # suffix-array pre-index size (small genomes)
_STAR_GENOME_GENERATE_RAM = "40317074816"   # byte cap for --limitGenomeGenerateRAM (~37 GiB)
_STAR_BAM_SORT_RAM = "27643756136"          # byte cap for --limitBAMsortRAM (~26 GiB)
# --outSJfilterIntronMaxVsReadN: max intron length allowed for junctions
# supported by 1, 2, 3+ reads respectively.
_STAR_SJ_FILTER_INTRON_MAX_VS_READN = ("100", "300", "500")

_BAM = "STAR_Aligned.sortedByCoord.out.bam"
_SJ = "STAR_SJ.out.tab"


def _is_gzipped(path: Path) -> bool:
    """Check if a file is gzip-compressed by reading the magic bytes."""
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def map_reads(config: AssemblyConfig) -> None:
    """Map RNA-seq reads to the genome using STAR."""
    wd = config.work_dir
    log.info("Running STAR read mapping...")

    index_dir = wd / "build-index"

    # Build genome index
    if not index_dir.exists():
        index_dir.mkdir()
        run_cmd(
            [
                "STAR",
                "--genomeSAindexNbases", _STAR_SA_INDEX_NBASES,
                "--limitGenomeGenerateRAM", _STAR_GENOME_GENERATE_RAM,
                "--runThreadN", str(config.num_cpu),
                "--runMode", "genomeGenerate",
                "--genomeDir", str(index_dir),
                "--genomeFastaFiles", str(config.genome),
            ],
            cwd=wd,
        )

    # Detect compressed input
    reads = config.reads_args_star
    zcat_args = []
    first_read_file = Path(reads[0])
    if first_read_file.suffix in (".gz", ".gzip") or _is_gzipped(first_read_file):
        zcat_args = ["--readFilesCommand", "zcat"]

    quality_args = ["--outQSconversionAdd", "-31"] if config.phred_quality == 64 else []

    max_intron_args = (
        ["--alignIntronMax", str(config.max_intron_len)]
        if config.max_intron_len
        else []
    )

    star_cmd = [
        "STAR",
        "--runThreadN", str(config.num_cpu),
        "--genomeDir", str(index_dir),
        "--alignEndsType", config.align_mode,
        "--readFilesIn", *reads,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outSJfilterIntronMaxVsReadN", *_STAR_SJ_FILTER_INTRON_MAX_VS_READN,
        "--alignIntronMin", str(config.min_intron_len),
        *max_intron_args,
        "--outFileNamePrefix", "STAR_",
        "--outSAMattributes", "All",
        "--outSAMattrIHstart", "0",
        "--outSAMstrandField", "intronMotif",
        "--limitBAMsortRAM", _STAR_BAM_SORT_RAM,
        *zcat_args,
        *quality_args,
    ]

    try:
        run_cmd(star_cmd, cwd=wd)
    except ExternalToolError:
        log.warning("STAR failed, falling back to STARlong")
        star_long_cmd = ["STARlong", *star_cmd[1:]]
        run_cmd(star_long_cmd, cwd=wd)

    # Report mapping rate from STAR log
    _log_mapping_rate(wd)

    # Splice junctions → hints + splice summary + diagnostic (shared with segemehl).
    generate_rnaseq_hints(
        wd / _SJ, wd / _BAM, config.genome, wd,
        diagnose=config.diagnose_softclips, source_label="STAR",
    )

    # Cleanup build index (large)
    shutil.rmtree(index_dir, ignore_errors=True)


def map_reads_auto(config: AssemblyConfig) -> None:
    """STAR read mapping, escalating to a segemehl re-map on extensive non-canonical splicing.

    STAR runs first — it is fast and it produces the soft-clip / splice diagnostic.
    If that diagnostic calls non-canonical splicing ``EXTENSIVE``, the reads are
    re-mapped with splice-agnostic segemehl, and that BAM becomes the one StringTie,
    the RNA-seq hints, and SL read-side detection consume (resolved via
    ``config.aligner_bam``): STAR's canonical-biased alignment otherwise mis-places
    reads across the non-canonical junctions, so the genome-guided assembly built on
    it is unreliable. ``--aligner star`` skips the escalation; ``--aligner segemehl``
    maps with segemehl from the start.
    """
    from eukan.assembly.segemehl import _BAM as _SEGEMEHL_READ_BAM
    from eukan.assembly.segemehl import map_reads_segemehl

    map_reads(config)

    seg_bam = config.work_dir / _SEGEMEHL_READ_BAM
    if _non_canonical_call(config.work_dir) == _SEGEMEHL_NON_CANONICAL_CALL:
        log.warning(
            "Non-canonical splicing EXTENSIVE — re-mapping the reads with segemehl "
            "(splice-agnostic) so genome-guided assembly and hints are not biased by "
            "STAR's canonical-splice alignment. Pass --aligner star to skip this."
        )
        map_reads_segemehl(config)
    elif seg_bam.exists():
        # A prior escalation's segemehl BAM is stale now that this run is
        # canonical-dominant; drop it so config.aligner_bam falls back to STAR.
        seg_bam.unlink(missing_ok=True)
        (config.work_dir / f"{_SEGEMEHL_READ_BAM}.bai").unlink(missing_ok=True)


def _non_canonical_call(work_dir: Path) -> str | None:
    """The ``non_canonical_splice`` verdict from the soft-clip diagnostic, if any.

    Reads ``softclip_diagnostic_summary.json`` (written after read mapping when
    ``--diagnose-softclips`` is on). Returns the call string
    (``EXTENSIVE`` / ``MODERATE`` / ``ABSENT``) or ``None`` when the file is
    absent or unreadable.
    """
    path = work_dir / Artifact.SOFTCLIP_DIAGNOSTIC.value
    try:
        data = json.loads(path.read_text())
        return data["verdict"]["non_canonical_splice"]["call"]
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _prefer_segemehl_for_transcripts(config: AssemblyConfig) -> bool:
    """True when transcript→genome mapping should use segemehl rather than STARlong.

    segemehl is chosen when the read aligner is explicitly segemehl, or when the
    soft-clip diagnostic called non-canonical splicing ``EXTENSIVE`` — both signal
    a splice landscape STARlong's canonical seeding handles poorly.
    """
    if config.aligner == "segemehl":
        return True
    return _non_canonical_call(config.work_dir) == _SEGEMEHL_NON_CANONICAL_CALL


def map_transcripts(config: AssemblyConfig) -> None:
    """Map de novo transcripts to the genome, routing on the splice landscape.

    Uses splice-agnostic segemehl ``-S`` when non-canonical splicing is extensive
    (or ``--aligner segemehl`` was chosen), else STARlong (with segemehl as a
    per-set fallback when STARlong errors or maps nothing). The output BAMs and
    downstream contract are identical either way.
    """
    if _prefer_segemehl_for_transcripts(config):
        log.info(
            "Mapping de novo transcripts with segemehl -S (non-canonical splicing "
            "extensive or --aligner segemehl); STARlong skipped."
        )
        _map_transcripts_segemehl(config)
    else:
        map_transcripts_star(config)
    _finalize_transcript_diagnostics(config)


def _finalize_transcript_diagnostics(config: AssemblyConfig) -> None:
    """Log unmapped-transcript counts and (when diagnosing) characterize poly-A.

    Path-agnostic: runs after either mapper over the files already on disk, so it
    is safe on a resumed/reused run. The unmapped FASTA is reported unconditionally
    (completeness); the poly-A characterization of the de novo transcript→genome BAM
    and the unmapped set is gated by ``diagnose_softclips`` (it is a soft-clip
    diagnostic) and written to ``polyA_diagnostic.json``, separate from the SL verdict.
    """
    from eukan.assembly.jaccard import _genome_stats
    from eukan.assembly.polya import (
        characterize_polya_bam,
        scan_fasta_polya,
        stats_to_dict,
        write_polya_section,
    )
    from eukan.assembly.segemehl import _GENOME_BAM_SUFFIX, _TRANSCRIPT_SETS, _resolve_query

    wd = config.work_dir
    for query_name, out_bam in _TRANSCRIPT_SETS:
        genome_bam = wd / out_bam
        if not genome_bam.exists():
            continue
        stem = out_bam[: -len(_GENOME_BAM_SUFFIX)]
        unmapped = wd / f"{stem}.unmapped_transcripts.fasta"

        # The unmapped FASTA is written only when the mapping actually runs; on a
        # resumed run with a complete BAM (or a run dir predating this feature) it
        # may be absent. Distinguish that from a genuine zero so the log/JSON never
        # claims "everything mapped" when the count is simply unavailable.
        if unmapped.exists():
            n_unmapped, n_unmapped_polya = scan_fasta_polya(unmapped)
            query = _resolve_query(wd, query_name)
            n_input = _genome_stats(query)[1] if query.exists() else 0
            pct = 100.0 * n_unmapped / n_input if n_input else 0.0
            log.info(
                "Unmapped de novo transcripts (%s): %d of %d (%.2f%%)%s",
                stem, n_unmapped, n_input, pct,
                f" -> {unmapped.name}" if n_unmapped else "",
            )
        else:
            log.info(
                "Unmapped de novo transcript FASTA not present for %s (reused BAM or "
                "older run dir); unmapped count unavailable this run.", stem,
            )

        if not config.diagnose_softclips:
            continue
        tx_stats = characterize_polya_bam(genome_bam, "transcripts")
        write_polya_section(wd, "transcripts", stats_to_dict(tx_stats))
        if unmapped.exists():
            write_polya_section(
                wd, "unmapped_transcripts",
                {"n_seqs": n_unmapped, "n_with_polyA_tail": n_unmapped_polya},
            )
        log.info(
            "Poly-A in de novo transcript mapping (%s): %d poly-A 3' soft-clips of %d "
            "(%.3f%%).",
            stem, tx_stats.n_polya, tx_stats.n_clips_examined, tx_stats.polya_pct_of_clips,
        )


def _map_transcripts_segemehl(config: AssemblyConfig) -> None:
    """segemehl ``-S`` transcript→genome mapping for every de novo set present.

    The segemehl primary path (vs. the STARlong-failure fallback): no STAR index
    is built. Per-set resume is handled inside
    :func:`segemehl.map_one_transcript_set_segemehl` (a complete BAM is reused).
    """
    from eukan.assembly.segemehl import (
        _TRANSCRIPT_SETS,
        _resolve_query,
        map_one_transcript_set_segemehl,
    )

    wd = config.work_dir
    sets = [
        (query, out_bam)
        for query_name, out_bam in _TRANSCRIPT_SETS
        if (query := _resolve_query(wd, query_name)).exists() and query.stat().st_size > 0
    ]
    if not sets:
        log.warning("No assembled transcripts found to map; skipping.")
        return
    for query, out_bam in sets:
        log.info("segemehl mapping transcripts %s -> %s ...", query.name, out_bam)
        map_one_transcript_set_segemehl(config, query, out_bam)


def map_transcripts_star(config: AssemblyConfig) -> None:
    """Map de novo assembled transcripts to the genome SPLICED with STARlong.

    Mapped with ``--alignIntronMax <max_intron_len>`` so cis-introns appear as
    ``N`` gaps — the splice structure :mod:`eukan.assembly.strand_correction` reads
    to homology-correct transcript strand — and ``--alignEndsType Local`` so the
    spliced leader still SOFT-CLIPS at the trans-splice acceptor (the SL-RNA locus
    is distal, not bridgeable within the bounded intron, so it is not split off),
    the signal SL detection (:mod:`eukan.assembly.sl_acceptors`) keys on. STARlong
    is the long-read STAR build, appropriate for long transcript queries; on failure
    or a zero-map result it falls back to segemehl ``-S`` (``-H 1``, memory-bounded).

    Produces one coordinate-sorted, indexed ``<stem>.genome.bam`` per de novo
    assembly present, with unmapped transcripts saved to
    ``<stem>.unmapped_transcripts.fasta``. Per-query resume: a BAM passing
    ``samtools quickcheck`` is left untouched.
    """
    from eukan.assembly.jaccard import _chr_bin_nbits, _genome_stats, _sa_index_nbases
    from eukan.assembly.segemehl import _GENOME_BAM_SUFFIX, _TRANSCRIPT_SETS, _resolve_query

    wd = config.work_dir
    sets = [
        (query, out_bam)
        for query_name, out_bam in _TRANSCRIPT_SETS
        if (query := _resolve_query(wd, query_name)).exists() and query.stat().st_size > 0
    ]
    if not sets:
        log.warning("No assembled transcripts found to map; skipping.")
        return

    index_dir = wd / "star_tx_index"
    total_len, n_seqs = _genome_stats(config.genome)
    shutil.rmtree(index_dir, ignore_errors=True)
    index_dir.mkdir()
    run_cmd(
        [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeDir", str(index_dir),
            "--genomeFastaFiles", str(config.genome),
            "--genomeSAindexNbases", _sa_index_nbases(total_len),
            "--genomeChrBinNbits", _chr_bin_nbits(total_len, n_seqs),
            "--limitGenomeGenerateRAM", _STAR_GENOME_GENERATE_RAM,
            "--runThreadN", str(config.num_cpu),
        ],
        cwd=wd,
    )
    try:
        for query, out_bam in sets:
            log.info("STAR mapping transcripts %s -> %s ...", query.name, out_bam)
            _star_map_one_transcript_set(config, index_dir, query, out_bam, _GENOME_BAM_SUFFIX)
    finally:
        # The index is regenerable; drop it so it doesn't linger in the run dir.
        shutil.rmtree(index_dir, ignore_errors=True)


def _count_mapped(bam_path: Path) -> int:
    """Mapped-read count from the BAM index (BAM must already be indexed)."""
    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return bam.mapped


def _count_fasta_records(path: Path) -> int:
    """Number of sequences in a (optionally gzipped) FASTA."""
    opener = gzip.open if _is_gzipped(path) else open
    n = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _star_input_reads(log_final: Path) -> int | None:
    """Parse 'Number of input reads' from a STAR ``Log.final.out``; None if missing."""
    try:
        for line in log_final.read_text().splitlines():
            if "Number of input reads" in line:
                return int(line.split("|")[-1].strip())
    except (OSError, ValueError):
        pass
    return None


def _star_map_one_transcript_set(
    config: AssemblyConfig, index_dir: Path, query: Path, out_bam: str, bam_suffix: str
) -> None:
    """Spliced-map one transcript FASTA → sorted, indexed *out_bam* (STARlong).

    Falls back to segemehl ``-S`` when STARlong errors or maps nothing — STARlong's
    short-read-tuned seeding can fail on long, multi-intron transcripts.
    """
    from eukan.assembly.segemehl import _bam_is_complete, map_one_transcript_set_segemehl

    wd = config.work_dir
    final = wd / out_bam
    if _bam_is_complete(final):
        log.info("Reusing %s; skipping STAR transcript mapping.", final.name)
        return

    stem = out_bam[: -len(bam_suffix)]
    prefix = f"startx_{stem}_"
    zcat_args = ["--readFilesCommand", "zcat"] if _is_gzipped(query) else []
    intron_max = str(config.max_intron_len) if config.max_intron_len else "1"
    starlong_cmd = [
        "STARlong",
        "--runThreadN", str(config.num_cpu),
        "--genomeDir", str(index_dir),
        "--readFilesIn", str(query),
        "--alignEndsType", "Local",        # soft-clip the spliced leader at the acceptor
        "--alignIntronMax", intron_max,    # spliced: recover cis-introns for strand inference
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outReadsUnmapped", "Fastx",
        "--outFileNamePrefix", prefix,
        "--limitBAMsortRAM", _STAR_BAM_SORT_RAM,
        *zcat_args,
    ]
    try:
        run_cmd(starlong_cmd, cwd=wd)
    except ExternalToolError:
        log.warning("STARlong failed mapping %s; falling back to segemehl -S.", query.name)
        map_one_transcript_set_segemehl(config, query, out_bam)
        return

    (wd / f"{prefix}Aligned.sortedByCoord.out.bam").rename(final)
    run_cmd(["samtools", "index", out_bam], cwd=wd)
    unmapped = wd / f"{prefix}Unmapped.out.mate1"
    if unmapped.exists():
        unmapped.rename(wd / f"{stem}.unmapped_transcripts.fasta")

    # STARlong is sometimes a plain short-read STAR build (no -DLONG_READS — e.g. a
    # bioconda SIMD variant that's a byte-identical copy of STAR), which mis-parses a
    # multi-record transcript FASTA into a single truncated read, then "maps" that one
    # read into a near-empty BAM. Detect that the run read far fewer reads than the FASTA
    # holds (parser breakage), as well as a zero-mapped result, and fall back to
    # segemehl -S, which maps long transcripts directly.
    n_input = _count_fasta_records(query)
    n_read = _star_input_reads(wd / f"{prefix}Log.final.out")
    n_mapped = _count_mapped(final)
    parser_broke = n_read is not None and n_input > 0 and n_read < n_input // 2
    if parser_broke or n_mapped == 0:
        log.warning(
            "STARlong read %s of %d transcripts in %s and mapped %d; "
            "falling back to segemehl -S.",
            f"only {n_read}" if n_read is not None else "an unknown number",
            n_input, query.name, n_mapped,
        )
        final.unlink(missing_ok=True)
        (wd / f"{out_bam}.bai").unlink(missing_ok=True)
        map_one_transcript_set_segemehl(config, query, out_bam)


def _log_mapping_rate(wd: Path) -> None:
    """Parse STAR Log.final.out and log the overall mapping rate."""
    log_file = wd / "STAR_Log.final.out"
    if not log_file.exists():
        return

    unique_pct = 0.0
    multi_pct = 0.0
    for line in log_file.read_text().splitlines():
        if "Uniquely mapped reads %" in line:
            unique_pct = float(line.split("|")[-1].strip().rstrip("%"))
        elif "% of reads mapped to multiple loci" in line:
            multi_pct = float(line.split("|")[-1].strip().rstrip("%"))

    total_pct = unique_pct + multi_pct
    if total_pct < 75:
        log.warning(
            "Detected low read mapping rate: %.1f%% (%.1f%% unique, %.1f%% multi)",
            total_pct, unique_pct, multi_pct,
        )
    else:
        log.info(
            "Read mapping rate: %.1f%% (%.1f%% unique, %.1f%% multi)",
            total_pct, unique_pct, multi_pct,
        )
