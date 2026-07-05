"""Tests for BAM N-CIGAR → STAR-format SJ.out.tab derivation and aligner wiring.

``sj_table_from_bam`` synthesizes a STAR-format junction table for aligners
without a native one (minimap2 today); once it does, the existing
``analyze_splice_sites`` / ``write_intron_hints`` are reused verbatim. These
tests build a tiny in-memory BAM + FASTA (pysam + BioPython only — no aligner
binaries needed) and assert the table, the splice summary, and the intron hints.
Also covers the (now single-aligner) minimap2 step wiring and the shared
coordinate-sort disk-full translation.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from eukan.assembly.align_hints import (
    analyze_splice_sites,
    sj_table_from_bam,
    write_intron_hints,
)
from eukan.assembly.pipeline import _steps_for, force_steps_from_run_flags
from eukan.exceptions import ExternalToolError
from eukan.settings import AssemblyConfig

# (start, end) 0-based half-open intron, n_unique (NH=1), n_multi (NH=2)
_JUNCTIONS = [
    (100, 160, 5, 2),   # GT-AG canonical  → motif 1, strand 1
    (300, 360, 8, 0),   # CG-AG non-canon  → motif 0, strand 0
    (500, 560, 4, 0),   # CT-AC rev-canon  → motif 2, strand 2
]
_FLANK = 10


def _write_fasta(path: Path, contigs: list[tuple[str, str]]) -> Path:
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n{seq}\n")
    return path


def _write_bam(path: Path, contigs, reads) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": n, "LN": ln} for n, ln in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["query_name"]
            a.query_sequence = r["query_sequence"]
            a.flag = r["flag"]
            a.reference_id = r["reference_id"]
            a.reference_start = r["reference_start"]
            a.mapping_quality = r["mapping_quality"]
            a.cigartuples = r["cigartuples"]
            a.set_tag("NH", r["nh"])
            out.write(a)
    return path


def _make_genome(length: int = 1000) -> str:
    seq = ["C"] * length
    def put(pos: int, bases: str) -> None:
        for i, ch in enumerate(bases):
            seq[pos + i] = ch
    put(100, "GT"); put(158, "AG")   # (100,160) GT-AG  # noqa: E702
    put(300, "CG"); put(358, "AG")   # (300,360) CG-AG  # noqa: E702
    put(500, "CT"); put(558, "AC")   # (500,560) CT-AC  # noqa: E702
    return "".join(seq)


def _reads():
    reads = []
    n = 0
    for start, end, n_uniq, n_multi in _JUNCTIONS:
        for nh, count in ((1, n_uniq), (2, n_multi)):
            for _ in range(count):
                reads.append(dict(
                    query_name=f"r{n}", query_sequence="A" * (2 * _FLANK),
                    flag=0, reference_id=0, reference_start=start - _FLANK,
                    mapping_quality=60,
                    cigartuples=[(0, _FLANK), (3, end - start), (0, _FLANK)],
                    nh=nh,
                ))
                n += 1
    # Sub-min-intron N-op (len 10) — must be length-filtered out.
    reads.append(dict(
        query_name=f"r{n}", query_sequence="A" * (2 * _FLANK),
        flag=0, reference_id=0, reference_start=690, mapping_quality=60,
        cigartuples=[(0, _FLANK), (3, 10), (0, _FLANK)], nh=1,
    ))
    return reads


def _build(tmp_path):
    genome = _write_fasta(tmp_path / "g.fa", [("chr1", _make_genome())])
    bam = _write_bam(tmp_path / "s.bam", [("chr1", 1000)], _reads())
    return bam, genome


def _sj_rows(sj_path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in sj_path.read_text().splitlines()]


def test_sj_table_from_bam(tmp_path):
    bam, genome = _build(tmp_path)
    sj = sj_table_from_bam(
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="mm2_SJ.out.tab",
    )
    rows = _sj_rows(sj)
    # Sorted by (chrom,start,end); the len-10 junction is filtered out.
    assert rows == [
        # chrom start end strand motif annotated n_unique n_multi overhang
        ["chr1", "101", "160", "1", "1", "0", "5", "2", "0"],  # GT-AG
        ["chr1", "301", "360", "0", "0", "0", "8", "0", "0"],  # CG-AG (non-canon)
        ["chr1", "501", "560", "2", "2", "0", "4", "0", "0"],  # CT-AC (rev-canon)
    ]


def test_sj_feeds_splice_summary(tmp_path):
    import json
    bam, genome = _build(tmp_path)
    sj = sj_table_from_bam(
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="mm2_SJ.out.tab",
    )
    analyze_splice_sites(sj, genome, tmp_path)
    summary = json.loads((tmp_path / "splice_site_summary.json").read_text())
    # Canonical via motif code; non-canonical CG-AG re-extracted from the genome.
    assert summary["GT-AG"] == {"count": 1, "unique_reads": 5}
    assert summary["CG-AG"] == {"count": 1, "unique_reads": 8}
    assert summary["CT-AC"] == {"count": 1, "unique_reads": 4}


def test_sj_feeds_intron_hints(tmp_path):
    bam, genome = _build(tmp_path)
    sj = sj_table_from_bam(
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="mm2_SJ.out.tab",
    )
    write_intron_hints(sj, tmp_path, "minimap2")
    lines = (tmp_path / "hints_introns.gff").read_text().splitlines()
    # score/mult = unique + multi; strand from SJ col4 (0→".").
    assert lines[0] == "chr1\tminimap2\tintron\t101\t160\t7\t+\t.\tmult=7;pri=4;src=E"
    assert lines[1] == "chr1\tminimap2\tintron\t301\t360\t8\t.\t.\tmult=8;pri=4;src=E"
    assert lines[2] == "chr1\tminimap2\tintron\t501\t560\t4\t-\t.\tmult=4;pri=4;src=E"


# ------------------------------------------------------------------
# Aligner-selection wiring (minimap2 is the sole aligner)
# ------------------------------------------------------------------

_CHAIN = [
    "minimap2", "trinity", "jaccard", "map_transcripts", "strand_correct",
    "defuse", "max_intron_split", "sl_detect", "sl_cut", "combinr",
]


def test_steps_for_lists_minimap2_first():
    steps = _steps_for()
    assert [s.name for s in steps] == _CHAIN
    assert steps[0].output == "minimap2_Aligned.sortedByCoord.out.bam"


def test_force_steps_from_run_minimap2_cascades():
    # --run-minimap2 cascades to the consumers of the read BAM (Trinity reads it
    # for genome-guided assembly; sl_detect reads it for SL detection), so the
    # transitive closure pulls in the whole Trinity track chain.
    assert force_steps_from_run_flags(run_minimap2=True) == [f"assembly/{s}" for s in _CHAIN]
    # --force re-runs the whole chain.
    assert force_steps_from_run_flags(force=True) == [f"assembly/{s}" for s in _CHAIN]


def test_assembly_config_aligner_fields(tmp_path):
    cfg = AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path,
        num_cpu=1, left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq",
    )
    assert cfg.aligner_bam == "minimap2_Aligned.sortedByCoord.out.bam"
    assert cfg.reads_args_minimap2 == [str(tmp_path / "l.fq"), str(tmp_path / "r.fq")]
    assert cfg.non_canonical == "auto"  # default


def test_non_canonical_field(tmp_path):
    cfg = AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path,
        num_cpu=1, non_canonical="force",
    )
    assert cfg.non_canonical == "force"


# ------------------------------------------------------------------
# Shared coordinate-sort disk-full translation (bam_utils)
# ------------------------------------------------------------------


def test_coordinate_sort_translates_disk_full(tmp_path, monkeypatch):
    """A samtools 'Illegal seek' (a full-disk write failure) is re-raised as a
    clear out-of-space error whose hint points at freeing space and notes that
    the unsorted BAM is reused on re-run."""
    from eukan.assembly import bam_utils

    (tmp_path / "unsorted.bam").write_bytes(b"BAMDATA")

    def boom(c1, c2, **kw):
        raise ExternalToolError(
            "samtools failed (exit 1)", tool="samtools", returncode=1, cmd=c2,
            stderr_snippet=(
                "[bam_sort_core] merging from 10 files...\n"
                'samtools sort: failed writing to "out.bam": Illegal seek'
            ),
        )

    monkeypatch.setattr(bam_utils, "run_piped", boom)

    with pytest.raises(ExternalToolError) as exc_info:
        bam_utils._coordinate_sort_and_filter(
            tmp_path / "unsorted.bam", "out.bam", tmp_path, 4
        )

    # Original samtools stderr preserved, plus an actionable disk hint.
    assert "Illegal seek" in exc_info.value.stderr_snippet
    assert exc_info.value.hint is not None
    assert "re-map" in exc_info.value.hint.lower() or "reuse" in exc_info.value.hint.lower()
