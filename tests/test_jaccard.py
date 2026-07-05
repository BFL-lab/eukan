"""Unit tests for eukan.assembly.jaccard (Trinity-style jaccard clipping)."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from eukan.assembly import jaccard as j
from eukan.assembly.jaccard import (
    Trough,
    _candidate_troughs,
    _group_and_pick_best,
    _partition_exons,
    _require_hills,
    clip_gff3,
    coverage_array,
    find_clip_points,
    iter_fragment_spans,
    jaccard_array,
    run_jaccard,
    split_fasta_record,
)
from eukan.settings import AssemblyConfig

# A 120 bp non-palindromic genome (seeded) so reverse-complement is distinct —
# essential for the minus-strand split test to actually exercise the sign.
GENOME = "".join(random.Random(42).choice("ACGT") for _ in range(120))


def _read_fasta(path) -> dict[str, str]:
    recs, name, seq = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                recs[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name is not None:
        recs[name] = "".join(seq)
    return recs


# --- jaccard math ----------------------------------------------------------


def test_jaccard_single_edge_fragments():
    # One fragment touches the left window edge, one the right, neither spans:
    # n_both=0, n_single=2 -> jaccard = 1/3.
    jac = jaccard_array([(1, 50), (110, 160)], 300)
    # window_lend = mid - 49; mid=60 -> window_lend=11, window_rend=110.
    assert jac[60] == round(1 / 3, 4)


def test_jaccard_no_fragments_yields_no_clip():
    # With no fragments the pseudocount makes empty regions read as high jaccard
    # (1.0), never a trough — so no false split on a contig with no read support.
    jac = jaccard_array([], 100)
    assert len(jac) == 101
    assert find_clip_points([], 100) == []


def test_jaccard_window_length_boundary():
    # A fragment exactly W long spans exactly one junction position fully.
    W = j._WINDOW
    jac = jaccard_array([(1, W)], 200)
    # at window_lend=1 (mid=1+(W-1)//2) this single frag spans both edges.
    mid = 1 + (W - 1) // 2
    assert jac[mid] == round((1 + 1) / (0 + 1 + 1), 4)  # n_both=1,n_single=0 -> 2/2...


def test_clean_break_produces_trough_and_one_clip():
    # Two well-covered clusters separated by a gap (> W): read pairs do not
    # bridge the gap, so jaccard troughs there and exactly one clip is called.
    left = [(i, i + 150) for i in range(1, 80)]       # covers ~1..230
    right = [(i, i + 150) for i in range(400, 480)]   # covers ~400..630
    frags = left + right
    L = 700
    jac = jaccard_array(frags, L)
    assert min(jac[231:400]) <= j._MAX_TROUGH_VAL  # a real trough in the gap
    clips = find_clip_points(frags, L)
    assert len(clips) == 1
    assert 230 < clips[0] < 400  # the clip lands inside the break


def test_well_supported_contig_has_no_clip():
    # A single dense cluster, every position bridged -> no trough, no clip.
    frags = [(i, i + 150) for i in range(1, 200)]
    assert find_clip_points(frags, 350) == []


# --- coverage-adaptive greediness ------------------------------------------


def _two_clusters(start_a, start_b, n=6, step=16, ins=150):
    """Two low-depth fragment clusters separated by a gap (a faint fusion)."""
    a = [(start_a + k * step, start_a + k * step + ins) for k in range(n)]
    b = [(start_b + k * step, start_b + k * step + ins) for k in range(n)]
    return a + b


def test_lowcov_fusion_split_only_when_greedy():
    # Two ~6x clusters with a gap: the junction jaccard troughs to ~0.2, above the
    # fixed 0.05 floor, so the Trinity-faithful pass (greed=0) misses it. The
    # coverage-adaptive gate (greed>0) widens the floor at this depth and splits it.
    frags = _two_clusters(1, 430)
    L = 700
    assert find_clip_points(frags, L, greed=0.0) == []
    clips = find_clip_points(frags, L, greed=1.5)
    assert len(clips) == 1
    assert 230 < clips[0] < 430  # the clip lands inside the gap


def test_greed_does_not_overclip_dense_or_change_highcov_break():
    # A dense single contig is never split, even when greedy (bridging holds).
    dense = [(i, i + 150) for i in range(1, 200)]
    assert find_clip_points(dense, 350, greed=1.5) == []
    # A high-coverage clean break is called identically with or without greed:
    # the achievable floor is far below 0.05 there, so the strict gate stands.
    hi = [(i, i + 150) for i in range(1, 80)] + [(i, i + 150) for i in range(400, 480)]
    assert find_clip_points(hi, 700, greed=0.0) == find_clip_points(hi, 700, greed=1.5)


# --- tunable detection thresholds ------------------------------------------


def test_find_clip_points_honors_max_trough():
    # A faint two-cluster fusion troughs to ~0.2: missed at the strict 0.05 floor
    # (greed=0), but a directly-raised max_trough gate catches it with NO coverage
    # adaptation — reproducing what greediness achieves at this depth.
    frags = _two_clusters(1, 430)
    L = 700
    assert find_clip_points(frags, L, greed=0.0) == []
    raised = find_clip_points(frags, L, greed=0.0, max_trough=0.30)
    assert raised == find_clip_points(frags, L, greed=1.5)
    assert len(raised) == 1


def test_find_clip_points_honors_min_delta():
    # A confident high-coverage clean break is normally one clip; requiring an
    # impossibly tall flanking hill (> the 1.0 jaccard ceiling) suppresses it.
    frags = [(i, i + 150) for i in range(1, 80)] + [(i, i + 150) for i in range(400, 480)]
    L = 700
    assert len(find_clip_points(frags, L)) == 1
    assert find_clip_points(frags, L, min_delta=1.5) == []


def _dip_jac(length: int, center: int, val: float) -> list[float]:
    arr = [1.0] * (length + 1)
    arr[center] = val
    return arr


def test_candidate_troughs_adaptive_threshold():
    L, center, dip = 600, 300, 0.12
    jac = _dip_jac(L, center, dip)
    # Low flank coverage: achievable floor rises above the dip -> candidate.
    low = _candidate_troughs(jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=[3] * (L + 1), greed=1.5)
    assert [t.pos for t in low] == [center]
    # High flank coverage: achievable floor stays below 0.05, the dip (0.12) is
    # NOT a clean fusion at that depth -> rejected (strict gate holds).
    high = _candidate_troughs(jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=[50] * (L + 1), greed=1.5)
    assert high == []
    # greed=0 (or no cov) reproduces the fixed-floor behaviour exactly.
    assert _candidate_troughs(jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=[3] * (L + 1), greed=0.0) == []
    assert _candidate_troughs(jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL) == []


def test_candidate_troughs_no_cov0_phantom_at_start():
    # The first scanned center (100) must use a real left-flank coverage, not the
    # cov[0]=0 sentinel; with uniform coverage its gate then matches an interior
    # center, so a dip the interior rejects is rejected at the start too.
    L = 600
    cov = [10] * (L + 1)
    start = _candidate_troughs(
        _dip_jac(L, 100, 0.10), j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5
    )
    mid = _candidate_troughs(
        _dip_jac(L, 300, 0.10), j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5
    )
    assert [t.pos for t in start] == [t.pos for t in mid] == []


def test_candidate_troughs_adaptive_cap():
    L, center = 600, 300
    cov = [1] * (L + 1)  # achievable 1/3 -> *1.5 = 0.5, capped at _MAX_ADAPTIVE_TROUGH
    below = _candidate_troughs(
        _dip_jac(L, center, j._MAX_ADAPTIVE_TROUGH - 0.02),
        j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5,
    )
    assert [t.pos for t in below] == [center]
    above = _candidate_troughs(
        _dip_jac(L, center, j._MAX_ADAPTIVE_TROUGH + 0.05),
        j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5,
    )
    assert above == []


def test_candidate_troughs_adaptive_cap_is_tunable():
    # The adaptive cap is a parameter: a 0.20 dip at very low depth is a candidate
    # under the default 0.30 cap, but a tightened cap (0.10) rejects it — the gate
    # is never allowed to relax that far, so low-coverage clipping stays stringent.
    L, center = 600, 300
    cov = [1] * (L + 1)  # achievable 1/3 -> *1.5 = 0.5, so the cap binds
    jac = _dip_jac(L, center, 0.20)
    assert [t.pos for t in _candidate_troughs(
        jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5)] == [center]
    assert _candidate_troughs(
        jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL, cov=cov, greed=1.5,
        max_adaptive_trough=0.10) == []


# --- coverage --------------------------------------------------------------


def test_coverage_array_inclusive():
    cov = coverage_array([(1, 3), (2, 5)], 6)
    assert cov == [0, 1, 2, 2, 1, 1, 0]  # index 0 unused


# --- trough detection helpers ----------------------------------------------


def _flat_jac(length: int, dips: list[int]) -> list[float]:
    arr = [0.0] + [1.0] * length
    for d in dips:
        arr[d] = 0.0
    return arr


def test_require_hills_needs_both_sides():
    jac = _flat_jac(600, [300])
    kept = _require_hills([Trough(300, 0.0, 1.0)], jac, j._TROUGH_WIN, j._MIN_JACCARD_DELTA)
    assert len(kept) == 1

    # No hill on the left (everything below the trough+delta there) -> dropped.
    one_sided = [0.0] * 301 + [1.0] * 300  # index 1..300 = 0.0, 301..600 = 1.0
    dropped = _require_hills([Trough(300, 0.0, 1.0)], one_sided, j._TROUGH_WIN, j._MIN_JACCARD_DELTA)
    assert dropped == []


def test_group_picks_single_deepest_within_window():
    jac = _flat_jac(700, [300, 450])  # two dips 150 bp apart (<= trough_win)
    troughs = _candidate_troughs(jac, j._TROUGH_WIN, j._MAX_TROUGH_VAL)
    best = _group_and_pick_best(troughs, j._TROUGH_WIN)
    assert len(best) == 1  # merged into one group


# --- fasta splitting -------------------------------------------------------


def test_split_fasta_record():
    segs = split_fasta_record("N" * 400, [100, 250], 25)
    assert [len(s) for s in segs] == [100, 150, 150]


def test_split_fasta_record_drops_short_pieces():
    segs = split_fasta_record("N" * 400, [10, 250], 25)  # first piece 10 bp < 25
    assert [len(s) for s in segs] == [240, 150]  # [11..250]=240, [251..400]=150


# --- GFF3 / GTF split path -------------------------------------------------


def test_partition_exons_plus_strand():
    # exons (1,20) and (31,50); spliced length 40; clip after spliced base 25.
    segs = _partition_exons([(1, 20), (31, 50)], [25], "+", 40)
    assert segs == [[(1, 20), (31, 35)], [(36, 50)]]


def test_partition_exons_minus_strand():
    # 5'->3' order is descending genomic for '-': exon (31,50) then (1,20).
    segs = _partition_exons([(31, 50), (1, 20)], [25], "-", 40)
    # segment1 = first 25 spliced bases from the high-genomic 5' end.
    assert segs == [[(31, 50), (16, 20)], [(1, 15)]]


# Two exons (1..40) + (51..90); spliced length 80 so a clip at 45 leaves both
# pieces above the 25 bp floor.
_GFF_TEMPLATE = (
    "##gff-version 3\n"
    "chr1\ttest\tgene\t1\t90\t.\t{s}\t.\tID=g1\n"
    "chr1\ttest\tmRNA\t1\t90\t.\t{s}\t.\tID=t1;Parent=g1\n"
    "chr1\ttest\texon\t1\t40\t.\t{s}\t.\tID=t1.e1;Parent=t1\n"
    "chr1\ttest\texon\t51\t90\t.\t{s}\t.\tID=t1.e2;Parent=t1\n"
)


def _orig_spliced(gff, genome):
    from eukan.gff import create_gff_db
    from eukan.gff.io import iter_assembled_sequences

    return dict(
        (m.id, s)
        for m, s in iter_assembled_sequences(
            create_gff_db(str(gff)), genome, child_featuretype="exon"
        )
    )["t1"]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_clip_gff3_splits_and_reextracts(tmp_path, strand):
    (tmp_path / "genome.fa").write_text(f">chr1\n{GENOME}\n")
    gff = tmp_path / "tx.gff3"
    gff.write_text(_GFF_TEMPLATE.format(s=strand))

    orig = _orig_spliced(gff, tmp_path / "genome.fa")
    assert len(orig) == 80

    out_gff = tmp_path / "clipped.gff3"
    out_fa = tmp_path / "clipped.fasta"
    clip_gff3(gff, tmp_path / "genome.fa", {"t1": [45]}, out_gff, out_fa)

    recs = _read_fasta(out_fa)
    assert set(recs) == {"t1.j1", "t1.j2"}
    assert recs["t1.j1"] == orig[:45]      # guards the minus-strand sign
    assert recs["t1.j2"] == orig[45:]
    assert recs["t1.j1"] + recs["t1.j2"] == orig


def test_split_models_at_clips_counts_and_passthrough(tmp_path):
    gff = tmp_path / "tx.gff3"
    gff.write_text(_GFF_TEMPLATE.format(s="+"))
    # spliced length 80; a clip at 45 splits into two >=25 bp pieces.
    out, n_split = j._split_models_at_clips(gff, {"t1": [45]}, 25)
    assert n_split == 1
    assert sorted(m.tid for m in out) == ["t1.j1", "t1.j2"]
    # no clip for this id -> passthrough, n_split 0.
    out2, n2 = j._split_models_at_clips(gff, {}, 25)
    assert n2 == 0 and [m.tid for m in out2] == ["t1"]


def test_write_spliced_fasta_matches_iter_assembled(tmp_path):
    # _write_spliced_fasta must use the SAME orientation as iter_assembled_sequences
    # (ascending-genomic concat, RC on '-'), so clip coords line up across paths.
    from eukan.assembly.jaccard import _parse_transcript_models, _write_transcript_models_gff3
    from eukan.gff import create_gff_db
    from eukan.gff.io import iter_assembled_sequences

    (tmp_path / "genome.fa").write_text(f">chr1\n{GENOME}\n")
    gtf = tmp_path / "stringtie.gtf"
    gtf.write_text(
        'chr1\tStringTie\texon\t1\t40\t.\t+\t.\ttranscript_id "p";\n'
        'chr1\tStringTie\texon\t51\t90\t.\t+\t.\ttranscript_id "p";\n'
        'chr1\tStringTie\texon\t1\t40\t.\t-\t.\ttranscript_id "m";\n'
        'chr1\tStringTie\texon\t51\t90\t.\t-\t.\ttranscript_id "m";\n'
    )
    models = _parse_transcript_models(gtf)
    mine_path = tmp_path / "spliced.fasta"
    j._write_spliced_fasta(models, tmp_path / "genome.fa", mine_path)
    mine = _read_fasta(mine_path)

    norm = tmp_path / "norm.gff3"
    _write_transcript_models_gff3(models, norm)
    ref = {
        m.id: s
        for m, s in iter_assembled_sequences(
            create_gff_db(str(norm)), tmp_path / "genome.fa", child_featuretype="exon"
        )
    }
    assert mine == ref
    assert mine["m"] != mine["p"]  # minus strand really is reverse-complemented


def test_resolve_stringtie_models_prefers_clipped(tmp_path):
    from eukan.assembly.jaccard import STRINGTIE_JACCARD_GFF3, resolve_stringtie_models

    (tmp_path / "stringtie.gtf").write_text("x")
    assert resolve_stringtie_models(tmp_path).name == "stringtie.gtf"
    # empty clipped file is ignored (size 0); non-empty one is preferred.
    (tmp_path / STRINGTIE_JACCARD_GFF3).write_text("")
    assert resolve_stringtie_models(tmp_path).name == "stringtie.gtf"
    (tmp_path / STRINGTIE_JACCARD_GFF3).write_text("##gff-version 3\n")
    assert resolve_stringtie_models(tmp_path).name == STRINGTIE_JACCARD_GFF3


def _header_bam(path, refs):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": ln} for n, ln in refs.items()]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass


def test_clip_stringtie_gtf_splits_fused_model(tmp_path, monkeypatch):
    (tmp_path / "genome.fa").write_text(f">chr1\n{GENOME}\n")
    gtf = tmp_path / "stringtie.gtf"
    gtf.write_text(
        'chr1\tStringTie\ttranscript\t1\t90\t.\t+\t.\ttranscript_id "STRG.1.1";\n'
        'chr1\tStringTie\texon\t1\t40\t.\t+\t.\ttranscript_id "STRG.1.1";\n'
        'chr1\tStringTie\texon\t51\t90\t.\t+\t.\ttranscript_id "STRG.1.1";\n'
    )
    bam = tmp_path / "jaccard_stringtie_Aligned.sortedByCoord.out.bam"
    _header_bam(bam, {"STRG.1.1": 80})  # spliced length of the transcript

    monkeypatch.setattr(j, "map_reads_to_transcripts", lambda cfg, fa, tag: bam)
    monkeypatch.setattr(j, "iter_fragment_spans", lambda b: iter([("STRG.1.1", [(1, 80)])]))
    # decouple from the jaccard math (tested elsewhere): clip the 80 bp tx at 45.
    monkeypatch.setattr(j, "find_clip_points", lambda spans, length, **kw: [45])

    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path,
        left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq", num_cpu=1,
    )
    n_in, n_split = j._clip_stringtie_gtf(config)
    assert (n_in, n_split) == (1, 1)
    out = (tmp_path / "stringtie.jaccard.gff3").read_text()
    assert out.count("\tmRNA\t") == 2  # split into two models
    assert "STRG.1.1.j1" in out and "STRG.1.1.j2" in out
    assert not (tmp_path / "stringtie.spliced.fasta").exists()  # intermediate cleaned


def test_clip_threads_config_knobs_into_find_clip_points(tmp_path, monkeypatch):
    # The config's jaccard detection knobs reach find_clip_points unchanged, so a
    # user-supplied stringency setting actually governs clip detection.
    (tmp_path / "genome.fa").write_text(f">chr1\n{GENOME}\n")
    gtf = tmp_path / "stringtie.gtf"
    gtf.write_text(
        'chr1\tStringTie\texon\t1\t40\t.\t+\t.\ttranscript_id "STRG.1.1";\n'
        'chr1\tStringTie\texon\t51\t90\t.\t+\t.\ttranscript_id "STRG.1.1";\n'
    )
    bam = tmp_path / "jaccard_stringtie_Aligned.sortedByCoord.out.bam"
    _header_bam(bam, {"STRG.1.1": 80})
    monkeypatch.setattr(j, "map_reads_to_transcripts", lambda cfg, fa, tag: bam)
    monkeypatch.setattr(j, "iter_fragment_spans", lambda b: iter([("STRG.1.1", [(1, 80)])]))
    captured: dict[str, float] = {}

    def _capture(spans, length, **kw):
        captured.update(kw)
        return []

    monkeypatch.setattr(j, "find_clip_points", _capture)

    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path,
        left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq", num_cpu=1,
        jaccard_greediness=2.0, jaccard_max_trough=0.08,
        jaccard_min_delta=0.5, jaccard_max_adaptive_trough=0.4,
    )
    j._clip_stringtie_gtf(config)
    assert captured == {
        "greed": 2.0, "max_trough": 0.08,
        "min_delta": 0.5, "max_adaptive_trough": 0.4,
    }


def test_clip_one_fasta_threads_config_knobs(tmp_path, monkeypatch):
    # The de novo FASTA path (the primary jaccard target — rnaSPAdes contigs) threads
    # the same config knobs into find_clip_points as the StringTie path, so reverting
    # it to the old greed-only call would be caught here too.
    src = tmp_path / "rnaspades.fasta"
    src.write_text(">c1\n" + "ACGT" * 50 + "\n")
    bam = tmp_path / "rnaspades.jaccard_Aligned.sortedByCoord.out.bam"
    _header_bam(bam, {"c1": 200})
    monkeypatch.setattr(j, "map_reads_to_transcripts", lambda cfg, fa, tag: bam)
    monkeypatch.setattr(j, "iter_fragment_spans", lambda b: iter([("c1", [(1, 200)])]))
    captured: dict[str, float] = {}

    def _capture(spans, length, **kw):
        captured.update(kw)
        return []

    monkeypatch.setattr(j, "find_clip_points", _capture)

    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path,
        left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq", num_cpu=1,
        jaccard_greediness=2.0, jaccard_max_trough=0.08,
        jaccard_min_delta=0.5, jaccard_max_adaptive_trough=0.4,
    )
    j._clip_one_fasta(config, src, tmp_path / "rnaspades.jaccard.fasta")
    assert captured == {
        "greed": 2.0, "max_trough": 0.08,
        "min_delta": 0.5, "max_adaptive_trough": 0.4,
    }


def test_assembly_config_jaccard_knob_defaults():
    # Defaults preserve the Trinity-faithful behaviour (no behavioural change unless
    # a knob is set), matching the jaccard.py module constants.
    cfg = AssemblyConfig(genome=Path("g.fa"))
    assert cfg.jaccard_greediness == 1.5
    assert cfg.jaccard_max_trough == j._MAX_TROUGH_VAL == 0.05
    assert cfg.jaccard_min_delta == j._MIN_JACCARD_DELTA == 0.35
    assert cfg.jaccard_max_adaptive_trough == j._MAX_ADAPTIVE_TROUGH == 0.30


def test_run_jaccard_clips_stringtie_when_present(tmp_path, monkeypatch):
    monkeypatch.setattr(j, "_clip_one_fasta", lambda cfg, src, out: (1, 1, 0))
    called: list[str] = []
    monkeypatch.setattr(j, "_clip_stringtie_gtf", lambda cfg: called.append("st") or (3, 1))

    (tmp_path / "rnaspades.fasta").write_text(">c1\nACGTACGT\n")
    (tmp_path / "stringtie.gtf").write_text("chr1\tx\texon\t1\t5\t.\t+\t.\ttranscript_id \"a\";\n")
    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path,
        left_reads=tmp_path / "l.fq", right_reads=tmp_path / "r.fq", num_cpu=1,
    )
    run_jaccard(config)
    assert called == ["st"]  # StringTie GTF clipped alongside the de novo FASTA


def test_clip_gff3_passthrough_when_no_clip(tmp_path):
    (tmp_path / "genome.fa").write_text(f">chr1\n{GENOME}\n")
    gff = tmp_path / "tx.gff3"
    gff.write_text(_GFF_TEMPLATE.format(s="+"))
    out_gff = tmp_path / "out.gff3"
    out_fa = tmp_path / "out.fasta"
    clip_gff3(gff, tmp_path / "genome.fa", {}, out_gff, out_fa)  # no clips
    assert set(_read_fasta(out_fa)) == {"t1"}


# --- run_jaccard wiring ----------------------------------------------------


def test_run_jaccard_noop_on_single_end(tmp_path, monkeypatch, caplog):
    def _boom(*a, **k):  # pragma: no cover - must not be called
        raise AssertionError("STAR must not run on single-end input")

    monkeypatch.setattr(j, "_clip_one_fasta", _boom)
    (tmp_path / "rnaspades.fasta").write_text(">c1\nACGTACGT\n")
    config = AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path,
        single_reads=tmp_path / "reads.fq", num_cpu=1,
    )
    with caplog.at_level("WARNING"):
        run_jaccard(config)
    assert "paired reads" in caplog.text
    # No .jaccard.fasta written.
    assert not (tmp_path / "rnaspades.jaccard.fasta").exists()


# --- fragment extraction from a BAM ----------------------------------------


def _write_bam(path, reads, ref="tx1", ref_len=800):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ref, "LN": ref_len}]}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for name, flag, start, length in reads:
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = name
            seg.flag = flag
            seg.reference_id = 0
            seg.reference_start = start
            seg.mapping_quality = 60
            seg.cigartuples = [(0, length)]  # length M
            out.write(seg)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))


def test_iter_fragment_spans_filters_pairs(tmp_path):
    bam = tmp_path / "tx.bam"
    # flags: 99 = paired,proper,first,mate-reverse ; 147 = paired,proper,second,reverse
    _write_bam(bam, [
        ("good", 99, 100, 50), ("good", 147, 250, 50),     # insert 200 -> kept (101,300)
        ("short", 99, 100, 50), ("short", 147, 100, 50),   # insert 50  -> dropped
        ("long", 99, 100, 50), ("long", 147, 650, 50),     # insert 600 -> dropped
        ("sec", 99 | 0x100, 120, 50),                       # secondary -> skipped
    ])
    spans = dict(iter_fragment_spans(bam))
    assert spans["tx1"] == [(101, 300)]
