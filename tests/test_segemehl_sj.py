"""Tests for the segemehl path: BAM N-CIGAR → STAR-format SJ.out.tab → reuse.

`sj_table_from_bam` is the only segemehl-specific logic; once it produces a
STAR-format junction table, the existing `analyze_splice_sites` /
`write_intron_hints` are reused verbatim. These tests build a tiny in-memory
BAM + FASTA (pysam + BioPython only — no segemehl/samtools binaries needed) and
assert the table, the splice summary, and the intron hints. Also covers the
aligner-selection wiring.
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
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="seg_SJ.out.tab",
    )
    rows = _sj_rows(sj)
    # Sorted by (chrom,start,end); the len-10 junction is filtered out.
    assert rows == [
        # chrom start end strand motif annotated n_unique n_multi overhang
        ["chr1", "101", "160", "1", "1", "0", "5", "2", "0"],  # GT-AG
        ["chr1", "301", "360", "0", "0", "0", "8", "0", "0"],  # CG-AG (non-canon)
        ["chr1", "501", "560", "2", "2", "0", "4", "0", "0"],  # CT-AC (rev-canon)
    ]


def test_segemehl_sj_feeds_splice_summary(tmp_path):
    import json
    bam, genome = _build(tmp_path)
    sj = sj_table_from_bam(
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="seg_SJ.out.tab",
    )
    analyze_splice_sites(sj, genome, tmp_path)
    summary = json.loads((tmp_path / "splice_site_summary.json").read_text())
    # Canonical via motif code; non-canonical CG-AG re-extracted from the genome.
    assert summary["GT-AG"] == {"count": 1, "unique_reads": 5}
    assert summary["CG-AG"] == {"count": 1, "unique_reads": 8}
    assert summary["CT-AC"] == {"count": 1, "unique_reads": 4}


def test_segemehl_sj_feeds_intron_hints(tmp_path):
    bam, genome = _build(tmp_path)
    sj = sj_table_from_bam(
        bam, genome, tmp_path, min_intron=20, max_intron=2000, out_name="seg_SJ.out.tab",
    )
    write_intron_hints(sj, tmp_path, "segemehl")
    lines = (tmp_path / "hints_introns.gff").read_text().splitlines()
    # score/mult = unique + multi; strand from SJ col4 (0→".").
    assert lines[0] == "chr1\tsegemehl\tintron\t101\t160\t7\t+\t.\tmult=7;pri=4;src=E"
    assert lines[1] == "chr1\tsegemehl\tintron\t301\t360\t8\t.\t.\tmult=8;pri=4;src=E"
    assert lines[2] == "chr1\tsegemehl\tintron\t501\t560\t4\t-\t.\tmult=4;pri=4;src=E"


# ------------------------------------------------------------------
# Aligner-selection wiring
# ------------------------------------------------------------------


def test_steps_for_selects_aligner():
    star = _steps_for("star")
    seg = _steps_for("segemehl")
    # Trinity (both de novo + genome-guided) replaces the old rnaspades +
    # stringtie steps; both Trinity FASTAs are jaccard-clipped then mapped.
    assert [s.name for s in star] == [
        "star", "trinity", "jaccard", "map_transcripts", "strand_correct", "defuse",
        "max_intron_split", "sl_detect", "sl_cut", "combinr",
    ]
    assert [s.name for s in seg] == [
        "segemehl", "trinity", "jaccard", "map_transcripts", "strand_correct", "defuse",
        "max_intron_split", "sl_detect", "sl_cut", "combinr",
    ]
    assert seg[0].output == "segemehl_Aligned.sortedByCoord.out.bam"


def test_force_steps_respects_active_aligner():
    # --run-<aligner> cascades to the consumers of the read BAM (Trinity reads it
    # for genome-guided assembly; sl_detect reads it for SL detection), so the
    # transitive closure pulls in the whole Trinity track chain; the active
    # aligner's name is the one in the returned chain.
    assert force_steps_from_run_flags(aligner="star", run_star=True) == [
        "assembly/star", "assembly/trinity", "assembly/jaccard",
        "assembly/map_transcripts", "assembly/strand_correct", "assembly/defuse",
        "assembly/max_intron_split", "assembly/sl_detect", "assembly/sl_cut",
        "assembly/combinr",
    ]
    assert force_steps_from_run_flags(aligner="segemehl", run_segemehl=True) == [
        "assembly/segemehl", "assembly/trinity", "assembly/jaccard",
        "assembly/map_transcripts", "assembly/strand_correct", "assembly/defuse",
        "assembly/max_intron_split", "assembly/sl_detect", "assembly/sl_cut",
        "assembly/combinr",
    ]
    # --force re-runs the active aligner's whole chain.
    assert force_steps_from_run_flags(aligner="segemehl", force=True) == [
        "assembly/segemehl", "assembly/trinity",
        "assembly/jaccard", "assembly/map_transcripts",
        "assembly/strand_correct", "assembly/defuse", "assembly/max_intron_split",
        "assembly/sl_detect", "assembly/sl_cut", "assembly/combinr",
    ]


def test_assembly_config_aligner_fields(tmp_path):
    cfg = AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path,
        num_cpu=1, aligner="segemehl",
        left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq",
    )
    assert cfg.aligner_bam == "segemehl_Aligned.sortedByCoord.out.bam"
    assert cfg.reads_args_segemehl == [
        "-q", str(tmp_path / "l.fq"), "-p", str(tmp_path / "r.fq"),
    ]
    star_cfg = cfg.model_copy(update={"aligner": "star"})
    assert star_cfg.aligner_bam == "STAR_Aligned.sortedByCoord.out.bam"


def test_steps_for_auto_uses_escalating_read_step():
    """`auto` keeps the `star` step name but uses the escalating read mapper."""
    from eukan.assembly.star import map_reads_auto

    steps = _steps_for("auto")
    assert steps[0].name == "star"
    assert steps[0].fn is map_reads_auto


def test_aligner_bam_auto_resolves_to_segemehl_once_escalated(tmp_path):
    """In auto mode, downstream reads STAR's BAM until a segemehl re-map appears."""
    cfg = AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path, num_cpu=1,
    )
    assert cfg.aligner == "auto"
    assert cfg.aligner_bam == "STAR_Aligned.sortedByCoord.out.bam"
    (tmp_path / "segemehl_Aligned.sortedByCoord.out.bam").write_text("escalated")
    assert cfg.aligner_bam == "segemehl_Aligned.sortedByCoord.out.bam"


def test_map_reads_segemehl_matches_h0_recipe(tmp_path, monkeypatch):
    """segemehl writes BAM to a file with the validated `-H 0` recipe, then a
    separate stage filters unmapped reads (`-F 4`) and coordinate-sorts.

    Writing the BAM via run_cmd (not piping segemehl's stdout) means
    segemehl's own exit status surfaces on failure — an OOM kill or a
    mate-count mismatch shows as a segemehl error, not a downstream
    "samtools: failed to read header".
    """
    from eukan.assembly import segemehl as seg_mod

    cmds: list[list[str]] = []
    pipes: list[tuple[list[str], list[str]]] = []
    monkeypatch.setattr(seg_mod, "run_cmd", lambda cmd, **kw: cmds.append(list(cmd)))
    monkeypatch.setattr(
        seg_mod, "run_piped", lambda c1, c2, **kw: pipes.append((list(c1), list(c2)))
    )
    monkeypatch.setattr(seg_mod, "sj_table_from_bam", lambda *a, **k: tmp_path / "sj.tab")
    monkeypatch.setattr(seg_mod, "generate_rnaseq_hints", lambda *a, **k: None)

    cfg = AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path,
        num_cpu=4, aligner="segemehl",
        left_reads=tmp_path / "l.fq.gz", right_reads=tmp_path / "r.fq.gz",
    )
    seg_mod.map_reads_segemehl(cfg)

    # Mapping invocation: reports all hits (-H 0), split mode (-S), writes a
    # BAM to a file (-b -o) rather than to stdout.
    seg_cmd = next(c for c in cmds if c and c[0] == "segemehl.x" and "-i" in c)
    assert seg_cmd[seg_cmd.index("-H") + 1] == "0"
    assert "-b" in seg_cmd
    assert seg_cmd[seg_cmd.index("-o") + 1].endswith(".bam")
    assert "-S" in seg_cmd

    # Post-process: drop unmapped reads (-F 4), then coordinate-sort, as one
    # piped stage; the sorted BAM is then indexed.
    view, sort = pipes[0]
    assert view[:2] == ["samtools", "view"]
    assert view[view.index("-F") + 1] == "4"
    assert sort[:2] == ["samtools", "sort"]
    assert any(c[:2] == ["samtools", "index"] for c in cmds)


def _segemehl_cfg(tmp_path):
    return AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, manifest_dir=tmp_path,
        num_cpu=4, aligner="segemehl",
        left_reads=tmp_path / "l.fq.gz", right_reads=tmp_path / "r.fq.gz",
    )


def test_map_reads_segemehl_resumes_and_frees_index_splits(tmp_path, monkeypatch):
    """A complete unsorted BAM from a prior run is reused (no re-mapping), and
    the index + split-read BEDs are deleted *before* the sort so they don't eat
    the disk the sort needs. This is the recovery path for a run that mapped
    successfully (multi-hour) but failed downstream in the sort.
    """
    from eukan.assembly import segemehl as seg_mod

    cmds: list[list[str]] = []
    pipes: list[tuple[list[str], list[str]]] = []
    monkeypatch.setattr(seg_mod, "run_cmd", lambda cmd, **kw: cmds.append(list(cmd)))
    monkeypatch.setattr(
        seg_mod, "run_piped", lambda c1, c2, **kw: pipes.append((list(c1), list(c2)))
    )
    monkeypatch.setattr(seg_mod, "sj_table_from_bam", lambda *a, **k: tmp_path / "sj.tab")
    monkeypatch.setattr(seg_mod, "generate_rnaseq_hints", lambda *a, **k: None)
    monkeypatch.setattr(seg_mod, "_bam_is_complete", lambda p: True)

    # Leftovers from the prior (failed-at-sort) run.
    (tmp_path / "segemehl_unsorted.bam").write_bytes(b"BAMDATA")
    (tmp_path / "segemehl.idx").write_bytes(b"INDEX")
    for suffix in (".sngl.bed", ".mult.bed", ".trns.txt"):
        (tmp_path / f"segemehl_splits{suffix}").write_bytes(b"SPLIT")

    seg_mod.map_reads_segemehl(_segemehl_cfg(tmp_path))

    # No segemehl invocation at all — neither index build nor (re-)mapping.
    assert not any(c and c[0] == "segemehl.x" for c in cmds)
    # Index + split BEDs were reclaimed before the sort.
    assert not (tmp_path / "segemehl.idx").exists()
    for suffix in (".sngl.bed", ".mult.bed", ".trns.txt"):
        assert not (tmp_path / f"segemehl_splits{suffix}").exists()
    # The sort + filter + index still run.
    assert pipes and pipes[0][0][:2] == ["samtools", "view"]
    assert any(c[:2] == ["samtools", "index"] for c in cmds)


def test_map_reads_segemehl_translates_disk_full(tmp_path, monkeypatch):
    """A samtools 'Illegal seek' (a full-disk write failure during the merge)
    is re-raised as a clear out-of-space error whose hint points at freeing
    space and notes that the unsorted BAM is reused on re-run.
    """
    from eukan.assembly import segemehl as seg_mod

    monkeypatch.setattr(seg_mod, "run_cmd", lambda cmd, **kw: None)
    monkeypatch.setattr(seg_mod, "sj_table_from_bam", lambda *a, **k: tmp_path / "sj.tab")
    monkeypatch.setattr(seg_mod, "generate_rnaseq_hints", lambda *a, **k: None)
    monkeypatch.setattr(seg_mod, "_bam_is_complete", lambda p: True)
    (tmp_path / "segemehl_unsorted.bam").write_bytes(b"BAMDATA")

    def boom(c1, c2, **kw):
        raise ExternalToolError(
            "samtools failed (exit 1)", tool="samtools", returncode=1, cmd=c2,
            stderr_snippet=(
                "[bam_sort_core] merging from 10 files...\n"
                'samtools sort: failed writing to "out.bam": Illegal seek'
            ),
        )

    monkeypatch.setattr(seg_mod, "run_piped", boom)

    with pytest.raises(ExternalToolError) as exc_info:
        seg_mod.map_reads_segemehl(_segemehl_cfg(tmp_path))

    # Original samtools stderr preserved, plus an actionable disk hint.
    assert "Illegal seek" in exc_info.value.stderr_snippet
    assert exc_info.value.hint is not None
    assert "re-map" in exc_info.value.hint.lower() or "reuse" in exc_info.value.hint.lower()
