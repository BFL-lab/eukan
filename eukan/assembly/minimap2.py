"""minimap2 read/transcript mapping — the single aligner for the assembly pipeline.

minimap2 handles all three mapping roles with one binary:

* short-read RNA-seq → genome, spliced (``-x splice:sr``);
* full-length Trinity transcripts → genome, spliced (``-x splice:hq``);
* ungapped read → transcript mapping for the jaccard fusion-clip (``-x sr``).

It has no native splice-junction table, so a STAR-format ``SJ.out.tab`` is derived
from the BAM's N-CIGAR junctions (:func:`align_hints.sj_table_from_bam`) and the
shared post-alignment processing is reused verbatim — downstream steps (GeneMark,
AUGUSTUS, transcript assembly, SL detection, strand correction) see the identical
contract STAR/segemehl produced, including the ``splice_site_summary.json`` that
lets AUGUSTUS allow non-canonical splice sites.

When the soft-clip diagnostic reports extensive non-canonical splicing, mapping
escalates to the non-canonical flags (``-J 0 -C 3 --splice-flank=no``), which drop
the canonical GT-AG bias so non-canonical introns (e.g. the dominant CG-AG introns
of diplonemids such as *Hemistasia*) are captured rather than misplaced. The
escalation is gated by ``--non-canonical`` (``auto`` = on the EXTENSIVE verdict,
``force`` = always, ``off`` = never).
"""

from __future__ import annotations

import json
from pathlib import Path

from eukan.assembly.align_hints import generate_rnaseq_hints, sj_table_from_bam
from eukan.assembly.bam_utils import (
    _GENOME_BAM_SUFFIX,
    _TRANSCRIPT_SETS,
    _bam_is_complete,
    _coordinate_sort_and_filter,
    _resolve_query,
    _write_unmapped_fasta,
)
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_piped
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_BAM = "minimap2_Aligned.sortedByCoord.out.bam"
_SJ = "minimap2_SJ.out.tab"

# Non-canonical-splice verdict (softclip_diagnostic_summary.json →
# verdict.non_canonical_splice.call) at which mapping escalates to the
# non-canonical minimap2 flags.
_NON_CANONICAL_CALL = "EXTENSIVE"

# Layered onto the splice preset when non-canonical splicing is extensive:
# -J 0 selects minimap2's original splice model (required for -C to take effect),
# -C 3 lowers the non-canonical (non-GT-AG) penalty below the preset default, and
# --splice-flank=no drops the human/mouse donor/acceptor flanking bias — together
# capturing non-canonical introns a canonical-tuned aligner would misplace.
_NON_CANONICAL_FLAGS = ["-J", "0", "-C", "3", "--splice-flank=no"]


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


def _use_non_canonical(config: AssemblyConfig) -> bool:
    """Whether to layer the non-canonical flags onto the splice preset.

    ``force`` always; ``off`` never; ``auto`` only when the soft-clip diagnostic
    called non-canonical splicing ``EXTENSIVE`` (that verdict is written by the
    read-mapping step, so it is available when transcript mapping runs later).
    """
    if config.non_canonical == "force":
        return True
    if config.non_canonical == "off":
        return False
    return _non_canonical_call(config.work_dir) == _NON_CANONICAL_CALL


def _log_mapping_rate(bam: Path, wd: Path) -> None:
    """Log the read mapping rate from the coordinate-sorted, indexed BAM.

    Replaces STAR's ``Log.final.out`` parse: minimap2 emits no such log, so the
    rate is read from the BAM index statistics (``--secondary=no`` means one
    record per read, so mapped / (mapped + unmapped) is the read mapping rate).
    """
    import pysam

    try:
        with pysam.AlignmentFile(str(bam), "rb") as af:
            mapped = af.mapped
            unmapped = af.unmapped
    except (OSError, ValueError):
        return
    total = mapped + unmapped
    if total == 0:
        return
    pct = 100.0 * mapped / total
    if pct < 75:
        log.warning(
            "Detected low read mapping rate: %.1f%% (%d of %d records mapped)",
            pct, mapped, total,
        )
    else:
        log.info(
            "Read mapping rate: %.1f%% (%d of %d records mapped)", pct, mapped, total
        )


def _emit_sj_and_hints(config: AssemblyConfig) -> None:
    """Derive a STAR-format SJ.out.tab from the read BAM and generate hints.

    minimap2 has no native junction table, so the SJ table is synthesized from
    the BAM's N-CIGAR junctions and the shared post-processing (splice summary,
    intron/coverage hints, soft-clip diagnostic) is reused — identical to the
    contract STAR/segemehl produced.
    """
    wd = config.work_dir
    sj = sj_table_from_bam(
        wd / _BAM, config.genome, wd,
        min_intron=config.min_intron_len,
        max_intron=config.max_intron_len,
        out_name=_SJ,
    )
    generate_rnaseq_hints(
        sj, wd / _BAM, config.genome, wd,
        diagnose=config.diagnose_softclips, source_label="minimap2",
    )


def _map_reads_once(config: AssemblyConfig, *, non_canonical: bool) -> None:
    """Map reads to the genome (splice:sr) into the sorted, indexed read BAM."""
    wd = config.work_dir
    cmd = ["minimap2", "-a", "-x", "splice:sr", "-t", str(config.num_cpu)]
    if config.max_intron_len:
        cmd += ["-G", str(config.max_intron_len)]
    cmd.append("--secondary=no")
    if non_canonical:
        cmd += _NON_CANONICAL_FLAGS
    cmd += [str(config.genome), *config.reads_args_minimap2]
    run_piped(
        cmd,
        ["samtools", "sort", "-@", str(config.num_cpu), "-o", _BAM, "-"],
        cwd=wd,
    )
    run_cmd(["samtools", "index", _BAM], cwd=wd)


def map_reads_minimap2(config: AssemblyConfig) -> None:
    """Map RNA-seq reads to the genome with minimap2 (splice:sr).

    Maps canonically first, generates the SJ table + hints + soft-clip diagnostic,
    then — in ``--non-canonical auto`` mode — re-maps with the non-canonical flags
    when the diagnostic calls non-canonical splicing ``EXTENSIVE``, so genome-guided
    assembly and hints are not biased by the canonical-splice alignment. The re-map
    overwrites the same BAM and rebuilds the hints (the verdict is idempotent, so it
    stands). ``--non-canonical force`` maps non-canonically from the start; ``off``
    never escalates.
    """
    wd = config.work_dir
    log.info("Running minimap2 read mapping (splice:sr)...")
    _map_reads_once(config, non_canonical=config.non_canonical == "force")
    _log_mapping_rate(wd / _BAM, wd)
    _emit_sj_and_hints(config)

    if config.non_canonical == "auto" and _non_canonical_call(wd) == _NON_CANONICAL_CALL:
        log.warning(
            "Non-canonical splicing EXTENSIVE — re-mapping the reads with the "
            "non-canonical minimap2 flags (-J 0 -C 3 --splice-flank=no) so "
            "genome-guided assembly and hints capture the non-canonical introns. "
            "Pass --non-canonical off to skip this."
        )
        _map_reads_once(config, non_canonical=True)
        _log_mapping_rate(wd / _BAM, wd)
        _emit_sj_and_hints(config)


def map_reads_to_transcripts(config: AssemblyConfig, fasta: Path, tag: str) -> Path:
    """Map paired reads to *fasta* ungapped (``-x sr``); return a sorted+indexed BAM.

    The direct replacement for STAR's ``--alignEndsType EndToEnd --alignIntronMax 1``
    read→transcript mapping used by the jaccard fusion-clip: transcripts have no
    introns, so ``-x sr`` (non-spliced short-read) is the right preset. minimap2
    indexes the transcript FASTA on the fly, so no index build/teardown is needed.
    """
    wd = config.work_dir
    prefix = f"jaccard_{tag}_"
    bam = wd / f"{prefix}Aligned.sortedByCoord.out.bam"
    cmd = [
        "minimap2", "-a", "-x", "sr", "-t", str(config.num_cpu), "--secondary=no",
        str(fasta), *config.reads_args_minimap2,
    ]
    run_piped(
        cmd,
        ["samtools", "sort", "-@", str(config.num_cpu), "-o", bam.name, "-"],
        cwd=wd,
    )
    run_cmd(["samtools", "index", bam.name], cwd=wd)
    return bam


def _map_one_transcript_set(
    config: AssemblyConfig, query: Path, out_bam: str, *, non_canonical: bool
) -> None:
    """Spliced-map one transcript FASTA → sorted, indexed *out_bam* (splice:hq).

    ``-x splice:hq`` handles full-length cDNA natively (no STARlong truncation),
    emitting cis-introns as ``N`` gaps — the splice structure strand_correction
    reads — and soft-clipping the spliced leader at the trans-splice acceptor (the
    SL-RNA locus is distal, not colinear, so it stays a soft-clip) — the signal SL
    detection keys on. Unmapped queries are pulled out before the ``-F 4`` filter.
    Per-query resume: a BAM passing ``samtools quickcheck`` is left untouched.
    """
    wd = config.work_dir
    final = wd / out_bam
    if _bam_is_complete(final):
        log.info("Reusing %s; skipping minimap2 transcript mapping.", final.name)
        return

    unsorted = wd / f"{out_bam}.unsorted.bam"
    cmd = ["minimap2", "-a", "-x", "splice:hq", "-t", str(config.num_cpu)]
    if config.max_intron_len:
        cmd += ["-G", str(config.max_intron_len)]
    cmd.append("--secondary=no")
    if non_canonical:
        cmd += _NON_CANONICAL_FLAGS
    cmd += [str(config.genome), str(query)]
    run_piped(cmd, ["samtools", "view", "-b", "-o", unsorted.name, "-"], cwd=wd)

    stem = out_bam[: -len(_GENOME_BAM_SUFFIX)]
    _write_unmapped_fasta(unsorted, wd / f"{stem}.unmapped_transcripts.fasta")
    _coordinate_sort_and_filter(unsorted, out_bam, wd, config.num_cpu)
    run_cmd(["samtools", "index", out_bam], cwd=wd)
    unsorted.unlink(missing_ok=True)


def map_transcripts_minimap2(config: AssemblyConfig) -> None:
    """Map de novo + genome-guided Trinity transcripts to the genome SPLICED.

    Uses ``-x splice:hq`` for every transcript set present, layering the
    non-canonical flags when the splice landscape is extensively non-canonical
    (or ``--non-canonical force``). Produces one coordinate-sorted, indexed
    ``<stem>.genome.bam`` per set, with unmapped transcripts saved to
    ``<stem>.unmapped_transcripts.fasta``.
    """
    wd = config.work_dir
    sets = [
        (query, out_bam)
        for query_name, out_bam in _TRANSCRIPT_SETS
        if (query := _resolve_query(wd, query_name)).exists() and query.stat().st_size > 0
    ]
    if not sets:
        log.warning("No assembled transcripts found to map; skipping.")
    else:
        non_canonical = _use_non_canonical(config)
        if non_canonical:
            log.info(
                "Mapping transcripts with the non-canonical minimap2 flags "
                "(non-canonical splicing extensive or --non-canonical force)."
            )
        for query, out_bam in sets:
            log.info("minimap2 mapping transcripts %s -> %s ...", query.name, out_bam)
            _map_one_transcript_set(config, query, out_bam, non_canonical=non_canonical)
    _finalize_transcript_diagnostics(config)


def _finalize_transcript_diagnostics(config: AssemblyConfig) -> None:
    """Log unmapped-transcript counts and (when diagnosing) characterize poly-A.

    Path-agnostic: runs after mapping over the files already on disk, so it is
    safe on a resumed/reused run. The unmapped FASTA is reported unconditionally
    (completeness); the poly-A characterization of the de novo transcript→genome
    BAM and the unmapped set is gated by ``diagnose_softclips`` and written to
    ``polyA_diagnostic.json``, separate from the SL verdict.
    """
    from eukan.assembly.jaccard import _genome_stats
    from eukan.assembly.polya import (
        characterize_polya_bam,
        scan_fasta_polya,
        stats_to_dict,
        write_polya_section,
    )

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
