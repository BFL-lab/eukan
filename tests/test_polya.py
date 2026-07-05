"""Unit tests for eukan.assembly.polya (poly-A / poly-T characterization)."""

from __future__ import annotations

import json
from pathlib import Path

import pysam

from eukan.assembly.bam_utils import _write_unmapped_fasta
from eukan.assembly.polya import (
    POLYA_DIAGNOSTIC,
    PolyAStats,
    characterize_polya_bam,
    classify_clip,
    scan_fasta_polya,
    stats_to_dict,
    tally_clip,
    write_polya_section,
)
from eukan.infra.artifacts import Artifact
from eukan.settings import AssemblyConfig


def _write_bam(path, contigs, reads):
    """Minimal BAM writer (mirrors tests/test_bam_diagnostic._write_bam)."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": n, "LN": ln} for n, ln in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["query_name"]
            a.query_sequence = r["query_sequence"]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["query_sequence"]))
            a.flag = r["flag"]
            a.reference_id = r["reference_id"]
            a.reference_start = r["reference_start"]
            a.mapping_quality = r.get("mapping_quality", 60)
            a.cigartuples = r["cigartuples"]
            out.write(a)
    return path


# --- classify_clip ---------------------------------------------------------


def test_classify_clip_polya_3p():
    assert classify_clip("3p", "A" * 20) == "polyA"
    assert classify_clip("3p", "A" * 17 + "GCT") == "polyA"  # 17/20 = 0.85 >= 0.8


def test_classify_clip_polyt_5p():
    assert classify_clip("5p", "T" * 20) == "polyT"
    assert classify_clip("5p", "T" * 17 + "GCA") == "polyT"


def test_classify_clip_wrong_side_is_none():
    # A-rich but at the mRNA 5' end is NOT a poly-A tail; T-rich at 3' is neither.
    assert classify_clip("5p", "A" * 20) is None
    assert classify_clip("3p", "T" * 20) is None


def test_classify_clip_length_and_fraction_thresholds():
    assert classify_clip("3p", "A" * 7) is None             # below min_len 8
    assert classify_clip("3p", "A" * 8) == "polyA"           # at min_len
    assert classify_clip("3p", "A" * 15 + "GGGGG") is None   # 15/20 = 0.75 < 0.8
    # tunable thresholds
    assert classify_clip("3p", "A" * 6, min_len=6) == "polyA"
    assert classify_clip("3p", "A" * 15 + "GGGGG", min_frac=0.7) == "polyA"


def test_classify_clip_mixed_is_none():
    assert classify_clip("3p", "ACGT" * 5) is None


# --- tally_clip / stats ----------------------------------------------------


def test_tally_clip_accumulates():
    s = PolyAStats()
    tally_clip(s, "3p", "A" * 20, "c1")
    tally_clip(s, "3p", "A" * 12, "c2")
    tally_clip(s, "5p", "T" * 15, "c1")
    tally_clip(s, "3p", "ACGT" * 5, "c1")  # examined but not poly-A
    tally_clip(s, "3p", "A" * 5, "c1")     # below min_len -> ignored entirely
    assert s.n_clips_examined == 4
    assert s.n_polya == 2
    assert s.n_polyt == 1
    assert s.polya_len_max == 20
    assert s.polya_mean_len == 16.0  # (20 + 12) / 2
    assert s.contigs_with_polya == {"c1", "c2"}
    assert round(s.polya_pct_of_clips, 1) == 50.0


def test_stats_to_dict_shape():
    s = PolyAStats()
    tally_clip(s, "3p", "A" * 20, "c1")
    d = stats_to_dict(s)
    assert d["n_polyA_3p"] == 1 and d["n_softclips_examined"] == 1
    assert set(d) == {
        "n_softclips_examined", "n_polyA_3p", "n_polyT_5p",
        "polyA_pct_of_softclips", "polyA_mean_len", "polyA_max_len",
        "n_contigs_with_polyA",
    }


# --- characterize_polya_bam (own BAM pass) ---------------------------------


def test_characterize_polya_bam_forward_and_reverse(tmp_path):
    bam = tmp_path / "tx.bam"
    reads = [
        # fwd, 3' soft-clip of 20 A's -> poly-A tail
        dict(query_name="t1", query_sequence="C" * 10 + "A" * 20, flag=0,
             reference_id=0, reference_start=5, cigartuples=[(0, 10), (4, 20)]),
        # fwd, 3' soft-clip that is NOT poly-A (examined but not counted)
        dict(query_name="t2", query_sequence="C" * 10 + "ACGT" * 5, flag=0,
             reference_id=0, reference_start=5, cigartuples=[(0, 10), (4, 20)]),
        # reverse read whose mRNA 3' tail is poly-A: stored RC -> LEADING 20S of T's;
        # _extract_clips re-orients to a "3p" A-rich clip.
        dict(query_name="t3", query_sequence="T" * 20 + "G" * 10, flag=16,
             reference_id=0, reference_start=5, cigartuples=[(4, 20), (0, 10)]),
    ]
    _write_bam(bam, [("c1", 100)], reads)
    stats = characterize_polya_bam(bam, "transcripts")
    assert stats.label == "transcripts"
    assert stats.n_clips_examined == 3
    assert stats.n_polya == 2          # t1 (fwd) + t3 (rev, re-oriented)
    assert stats.n_polyt == 0
    assert stats.contigs_with_polya == {"c1"}


# --- scan_fasta_polya ------------------------------------------------------


def test_scan_fasta_polya(tmp_path):
    fa = tmp_path / "u.fasta"
    records = {
        "a": "ACGT" * 20 + "A" * 16,  # 3' poly-A tail (sense)
        "b": "ACGT" * 20,             # no tail (ends ...ACGT)
        "c": "GGGG" + "A" * 8,        # short, trailing 8 A's
        "d": "T" * 12 + "ACGT" * 20,  # 5' poly-T (antisense poly-A on unstranded data)
        "e": "ACG",                   # below min_len -> skipped, never counted
    }
    fa.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    n, n_polya = scan_fasta_polya(fa)
    assert n == 5
    assert n_polya == 3  # a (3'A), c (3'A), d (5'T antisense)


def test_has_section(tmp_path):
    from eukan.assembly.polya import has_section

    assert has_section(tmp_path, "reads") is False  # no file yet
    write_polya_section(tmp_path, "reads", {"n_polyA_3p": 1})
    assert has_section(tmp_path, "reads") is True
    assert has_section(tmp_path, "transcripts") is False
    # corrupt JSON -> treated as absent (never crashes a producer)
    (tmp_path / POLYA_DIAGNOSTIC).write_text("{not json")
    assert has_section(tmp_path, "reads") is False


# --- write_polya_section (JSON merge) --------------------------------------


def test_write_polya_section_merges(tmp_path):
    p1 = write_polya_section(tmp_path, "reads", {"n_polyA_3p": 5})
    assert p1.name == POLYA_DIAGNOSTIC
    write_polya_section(tmp_path, "transcripts", {"n_polyA_3p": 0})
    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert data == {"reads": {"n_polyA_3p": 5}, "transcripts": {"n_polyA_3p": 0}}
    # re-writing a section overwrites just that key
    write_polya_section(tmp_path, "reads", {"n_polyA_3p": 9})
    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert data["reads"] == {"n_polyA_3p": 9} and data["transcripts"] == {"n_polyA_3p": 0}


# --- unmapped-transcript capture -------------------------------------------


def test_write_unmapped_fasta_extracts_only_unmapped(tmp_path):
    bam = tmp_path / "unsorted.bam"
    reads = [
        # mapped transcript -> excluded
        dict(query_name="tx_M", query_sequence="ACGT" * 10, flag=0,
             reference_id=0, reference_start=5, cigartuples=[(0, 40)]),
        # unmapped transcript -> written
        dict(query_name="tx_U", query_sequence="TTTTAAAACCCCGGGG", flag=4,
             reference_id=-1, reference_start=-1, cigartuples=None),
    ]
    _write_bam(bam, [("c1", 100)], reads)
    out = tmp_path / "trinity-denovo.unmapped_transcripts.fasta"
    n = _write_unmapped_fasta(bam, out)
    assert n == 1
    text = out.read_text()
    assert ">tx_U" in text and "TTTTAAAACCCCGGGG" in text
    assert "tx_M" not in text


def test_map_one_transcript_set_captures_unmapped(tmp_path, monkeypatch):
    # Drive the REAL minimap2._map_one_transcript_set body: the unmapped FASTA must be
    # extracted from the unsorted BAM BEFORE the -F 4 sort/filter drops the records.
    from eukan.assembly import minimap2

    def fake_piped(cmd1, cmd2, **kw):
        # emulate `minimap2 ... | samtools view -b -o <unsorted> -`
        out = Path(cmd2[cmd2.index("-o") + 1])
        _write_bam(tmp_path / out.name, [("c1", 100)], [
            dict(query_name="tx_M", query_sequence="ACGT" * 10, flag=0,
                 reference_id=0, reference_start=5, cigartuples=[(0, 40)]),
            dict(query_name="tx_U", query_sequence="GGGGCCCCAAAATTTT", flag=4,
                 reference_id=-1, reference_start=-1, cigartuples=None),
        ])
        return ""

    monkeypatch.setattr(minimap2, "run_piped", fake_piped)
    monkeypatch.setattr(minimap2, "run_cmd", lambda *a, **k: None)
    monkeypatch.setattr(minimap2, "_coordinate_sort_and_filter", lambda *a, **k: None)
    monkeypatch.setattr(minimap2, "_bam_is_complete", lambda _p: False)

    (tmp_path / "genome.fa").write_text(">c1\n" + "ACGT" * 100 + "\n")
    query = tmp_path / "trinity-denovo.fasta"
    query.write_text(">tx_M\nACGT\n>tx_U\nGGGG\n")
    config = AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=1)
    minimap2._map_one_transcript_set(
        config, query, "trinity-denovo.genome.bam", non_canonical=False
    )

    fa = tmp_path / "trinity-denovo.unmapped_transcripts.fasta"
    assert fa.exists()
    text = fa.read_text()
    assert ">tx_U" in text and "GGGGCCCCAAAATTTT" in text
    assert "tx_M" not in text


# --- align_hints integration (the "reads" section) -------------------------


def _reads_bam_with_polya(tmp_path):
    """A small read BAM whose only soft-clips are 3' poly-A tails + an N-only genome."""
    (tmp_path / "g.fa").write_text(">c1\n" + "N" * 400 + "\n")
    reads = [
        dict(query_name=f"r{i}", query_sequence="ACGTACGTAC" + "A" * 20, flag=0,
             reference_id=0, reference_start=20 + i * 30, mapping_quality=60,
             cigartuples=[(0, 10), (4, 20)])
        for i in range(3)
    ]
    bam = _write_bam(tmp_path / "reads.bam", [("c1", 400)], reads)
    return bam, tmp_path / "g.fa"


def test_run_softclip_diagnostic_writes_reads_polya_section(tmp_path):
    from eukan.assembly.align_hints import run_softclip_diagnostic

    bam, genome = _reads_bam_with_polya(tmp_path)
    run_softclip_diagnostic(bam, genome, tmp_path)

    assert (tmp_path / Artifact.SOFTCLIP_DIAGNOSTIC.value).exists()
    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert data["reads"]["n_polyA_3p"] == 3


def test_run_softclip_diagnostic_backfills_reads_section_on_resume(tmp_path):
    from eukan.assembly.align_hints import run_softclip_diagnostic

    bam, genome = _reads_bam_with_polya(tmp_path)
    run_softclip_diagnostic(bam, genome, tmp_path)            # cold: writes both
    (tmp_path / POLYA_DIAGNOSTIC).unlink()                    # simulate pre-feature run dir
    run_softclip_diagnostic(bam, genome, tmp_path)            # summary exists -> backfill

    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert data["reads"]["n_polyA_3p"] == 3  # re-created via the poly-A-only pass


# --- minimap2._finalize_transcript_diagnostics -----------------------------


def _finalize_setup(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "ACGT" * 200 + "\n")
    (tmp_path / "trinity-denovo.fasta").write_text(">t1\n" + "ACGT" * 50 + "\n")
    _write_bam(tmp_path / "trinity-denovo.genome.bam", [("c1", 800)], [
        dict(query_name="t1", query_sequence="C" * 10 + "A" * 20, flag=0,
             reference_id=0, reference_start=5, cigartuples=[(0, 10), (4, 20)]),
    ])


def test_finalize_transcript_diagnostics_writes_sections(tmp_path):
    from eukan.assembly.minimap2 import _finalize_transcript_diagnostics

    _finalize_setup(tmp_path)
    (tmp_path / "trinity-denovo.unmapped_transcripts.fasta").write_text(
        ">u1\n" + "ACGT" * 20 + "A" * 16 + "\n"
    )
    config = AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=1)
    _finalize_transcript_diagnostics(config)

    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert data["transcripts"]["n_polyA_3p"] == 1
    assert data["unmapped_transcripts"] == {"n_seqs": 1, "n_with_polyA_tail": 1}


def test_finalize_transcript_diagnostics_respects_diagnose_off(tmp_path):
    from eukan.assembly.minimap2 import _finalize_transcript_diagnostics

    _finalize_setup(tmp_path)
    (tmp_path / "trinity-denovo.unmapped_transcripts.fasta").write_text(">u1\nACGT\n")
    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=1,
        diagnose_softclips=False,
    )
    _finalize_transcript_diagnostics(config)
    assert not (tmp_path / POLYA_DIAGNOSTIC).exists()  # no poly-A work when diagnosing off


def test_finalize_transcript_diagnostics_handles_absent_unmapped_fasta(tmp_path):
    from eukan.assembly.minimap2 import _finalize_transcript_diagnostics

    _finalize_setup(tmp_path)  # NO unmapped FASTA on disk (e.g. reused BAM)
    config = AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=1)
    _finalize_transcript_diagnostics(config)

    data = json.loads((tmp_path / POLYA_DIAGNOSTIC).read_text())
    assert "transcripts" in data
    # No misleading zeroed "unmapped" section when the count is simply unavailable.
    assert "unmapped_transcripts" not in data
