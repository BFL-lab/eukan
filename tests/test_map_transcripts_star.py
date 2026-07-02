"""Unit tests for transcript->genome mapping (STARlong-spliced command construction)."""

from __future__ import annotations

import json
from pathlib import Path

from eukan.assembly import segemehl, star
from eukan.infra.artifacts import Artifact
from eukan.settings import AssemblyConfig


def _config(tmp_path, **kw):
    (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 100 + "\n")
    return AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=4, **kw)


def _mock(monkeypatch, tmp_path, *, mapped: int = 1):
    """Record run_cmd calls; materialize the mapper's sorted BAM so the rename succeeds.

    *mapped* stands in for the post-mapping ``_count_mapped`` so the segemehl fallback
    fires only when we want it to (mapped == 0).
    """
    cmds: list[list[str]] = []

    def fake(cmd, **kw):
        cmds.append(cmd)
        if cmd[0] in ("STAR", "STARlong") and "--outFileNamePrefix" in cmd:
            prefix = cmd[cmd.index("--outFileNamePrefix") + 1]
            (tmp_path / f"{prefix}Aligned.sortedByCoord.out.bam").write_text("bam")

    monkeypatch.setattr(star, "run_cmd", fake)
    monkeypatch.setattr(star, "_count_mapped", lambda _p: mapped)
    return cmds


def _map_cmds(cmds):
    return [c for c in cmds if c and c[0] in ("STAR", "STARlong") and "--readFilesIn" in c]


def test_maps_both_trinity_tracks(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    (tmp_path / "trinity-gg.fasta").write_text(">t\nACGTACGT\n")
    cmds = _mock(monkeypatch, tmp_path)

    star.map_transcripts_star(_config(tmp_path))

    assert sum(1 for c in cmds if "genomeGenerate" in c) == 1  # one genome index
    queries = {Path(c[c.index("--readFilesIn") + 1]).name for c in _map_cmds(cmds)}
    assert queries == {"trinity-denovo.fasta", "trinity-gg.fasta"}
    assert (tmp_path / "trinity-denovo.genome.bam").exists()
    assert (tmp_path / "trinity-gg.genome.bam").exists()


def test_starlong_flags_spliced_local(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    cmds = _mock(monkeypatch, tmp_path)

    star.map_transcripts_star(_config(tmp_path, max_intron_len=5000))

    (cmd,) = _map_cmds(cmds)
    assert cmd[0] == "STARlong"                                 # long-read STAR build
    assert cmd[cmd.index("--alignEndsType") + 1] == "Local"     # soft-clip the SL
    assert cmd[cmd.index("--alignIntronMax") + 1] == "5000"     # spliced, bounded intron
    assert cmd[cmd.index("--outReadsUnmapped") + 1] == "Fastx"
    assert ["samtools", "index", "trinity-denovo.genome.bam"] in cmds


def test_prefers_jaccard_clipped_query(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGT\n")
    (tmp_path / "trinity-denovo.jaccard.fasta").write_text(">t\nACGTACGT\n")
    cmds = _mock(monkeypatch, tmp_path)

    star.map_transcripts_star(_config(tmp_path))

    (cmd,) = _map_cmds(cmds)
    assert Path(cmd[cmd.index("--readFilesIn") + 1]).name == "trinity-denovo.jaccard.fasta"
    assert (tmp_path / "trinity-denovo.genome.bam").exists()  # output keyed to the stem


def test_unmapped_captured(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGT\n")

    def fake(cmd, **kw):
        if cmd[0] in ("STAR", "STARlong") and "--outFileNamePrefix" in cmd:
            p = cmd[cmd.index("--outFileNamePrefix") + 1]
            (tmp_path / f"{p}Aligned.sortedByCoord.out.bam").write_text("bam")
            (tmp_path / f"{p}Unmapped.out.mate1").write_text(">u\nACGT\n")

    monkeypatch.setattr(star, "run_cmd", fake)
    monkeypatch.setattr(star, "_count_mapped", lambda _p: 1)
    star.map_transcripts_star(_config(tmp_path))

    assert (tmp_path / "trinity-denovo.unmapped_transcripts.fasta").exists()


def test_segemehl_fallback_on_zero_map(tmp_path, monkeypatch):
    """STARlong mapping nothing falls back to segemehl -S for that transcript set."""
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    _mock(monkeypatch, tmp_path, mapped=0)

    called: list[str] = []
    monkeypatch.setattr(
        segemehl, "map_one_transcript_set_segemehl",
        lambda config, query, out_bam: called.append(out_bam),
    )

    star.map_transcripts_star(_config(tmp_path))

    assert called == ["trinity-denovo.genome.bam"]


def test_segemehl_fallback_on_parser_breakage(tmp_path, monkeypatch):
    """A short-read STAR build standing in for STARlong mis-parses a multi-record FASTA
    into a single read; that triggers the segemehl fallback even though the one read maps.
    """
    # 10 transcripts, but STARlong's Log.final.out will report only 1 input read.
    (tmp_path / "trinity-gg.fasta").write_text(
        "".join(f">t{i}\nACGTACGTACGT\n" for i in range(10))
    )

    def fake(cmd, **kw):
        if cmd[0] in ("STAR", "STARlong") and "--outFileNamePrefix" in cmd:
            p = cmd[cmd.index("--outFileNamePrefix") + 1]
            (tmp_path / f"{p}Aligned.sortedByCoord.out.bam").write_text("bam")
            (tmp_path / f"{p}Log.final.out").write_text(
                "                          Number of input reads |\t1\n"
            )

    monkeypatch.setattr(star, "run_cmd", fake)
    monkeypatch.setattr(star, "_count_mapped", lambda _p: 1)  # the single read "mapped"

    called: list[str] = []
    monkeypatch.setattr(
        segemehl, "map_one_transcript_set_segemehl",
        lambda config, query, out_bam: called.append(out_bam),
    )

    star.map_transcripts_star(_config(tmp_path))

    assert called == ["trinity-gg.genome.bam"]


def test_no_assemblies_is_noop(tmp_path, monkeypatch):
    cmds = _mock(monkeypatch, tmp_path)
    star.map_transcripts_star(_config(tmp_path))
    assert cmds == []


def _write_verdict(tmp_path, call):
    (tmp_path / Artifact.SOFTCLIP_DIAGNOSTIC.value).write_text(
        json.dumps({"verdict": {"non_canonical_splice": {"call": call}}})
    )


def _record_dispatch(monkeypatch):
    chosen: list[str] = []
    monkeypatch.setattr(star, "map_transcripts_star", lambda _c: chosen.append("starlong"))
    monkeypatch.setattr(star, "_map_transcripts_segemehl", lambda _c: chosen.append("segemehl"))
    return chosen


class TestMapTranscriptsDispatch:
    """``map_transcripts`` routes to STARlong or segemehl on the splice landscape."""

    def test_default_uses_starlong(self, tmp_path, monkeypatch):
        chosen = _record_dispatch(monkeypatch)
        star.map_transcripts(_config(tmp_path))
        assert chosen == ["starlong"]

    def test_aligner_segemehl_uses_segemehl(self, tmp_path, monkeypatch):
        chosen = _record_dispatch(monkeypatch)
        star.map_transcripts(_config(tmp_path, aligner="segemehl"))
        assert chosen == ["segemehl"]

    def test_extensive_non_canonical_uses_segemehl(self, tmp_path, monkeypatch):
        _write_verdict(tmp_path, "EXTENSIVE")
        chosen = _record_dispatch(monkeypatch)
        star.map_transcripts(_config(tmp_path))
        assert chosen == ["segemehl"]

    def test_moderate_non_canonical_uses_starlong(self, tmp_path, monkeypatch):
        """Only EXTENSIVE routes to segemehl; MODERATE stays on STARlong."""
        _write_verdict(tmp_path, "MODERATE")
        chosen = _record_dispatch(monkeypatch)
        star.map_transcripts(_config(tmp_path))
        assert chosen == ["starlong"]

    def test_unreadable_verdict_uses_starlong(self, tmp_path, monkeypatch):
        (tmp_path / Artifact.SOFTCLIP_DIAGNOSTIC.value).write_text("{not json")
        chosen = _record_dispatch(monkeypatch)
        star.map_transcripts(_config(tmp_path))
        assert chosen == ["starlong"]


class TestReadMappingEscalation:
    """``map_reads_auto`` re-maps the reads with segemehl when non-canonical EXTENSIVE."""

    def _patch(self, monkeypatch, verdict=None):
        calls: list[str] = []

        def fake_star(config):
            calls.append("star")
            if verdict is not None:
                _write_verdict(config.work_dir, verdict)

        monkeypatch.setattr(star, "map_reads", fake_star)
        # map_reads_segemehl is imported lazily from the segemehl module, so patch
        # it there (the function-local import resolves the attribute at call time).
        monkeypatch.setattr(
            segemehl, "map_reads_segemehl", lambda config: calls.append("segemehl")
        )
        return calls

    def test_escalates_on_extensive(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="EXTENSIVE")
        star.map_reads_auto(_config(tmp_path))
        assert calls == ["star", "segemehl"]

    def test_no_escalation_on_moderate(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="MODERATE")
        star.map_reads_auto(_config(tmp_path))
        assert calls == ["star"]

    def test_no_escalation_without_verdict(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict=None)
        star.map_reads_auto(_config(tmp_path))
        assert calls == ["star"]

    def test_stale_segemehl_bam_dropped_when_not_extensive(self, tmp_path, monkeypatch):
        from eukan.assembly.segemehl import _BAM as seg

        (tmp_path / seg).write_text("stale")
        (tmp_path / f"{seg}.bai").write_text("stale")
        calls = self._patch(monkeypatch, verdict="MODERATE")
        star.map_reads_auto(_config(tmp_path))
        assert calls == ["star"]
        assert not (tmp_path / seg).exists()  # stale escalation BAM cleared
        assert not (tmp_path / f"{seg}.bai").exists()


class TestSegemehlPrimaryPath:
    """The segemehl-primary transcript path maps each set, no STAR index built."""

    def test_maps_each_set(self, tmp_path, monkeypatch):
        (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
        (tmp_path / "trinity-gg.fasta").write_text(">t\nACGTACGT\n")
        called: list[tuple[str, str]] = []
        monkeypatch.setattr(
            segemehl, "map_one_transcript_set_segemehl",
            lambda config, query, out_bam: called.append((query.name, out_bam)),
        )
        star._map_transcripts_segemehl(_config(tmp_path))
        assert called == [
            ("trinity-denovo.fasta", "trinity-denovo.genome.bam"),
            ("trinity-gg.fasta", "trinity-gg.genome.bam"),
        ]

    def test_prefers_jaccard_clipped_query(self, tmp_path, monkeypatch):
        (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGT\n")
        (tmp_path / "trinity-denovo.jaccard.fasta").write_text(">t\nACGTACGT\n")
        called: list[str] = []
        monkeypatch.setattr(
            segemehl, "map_one_transcript_set_segemehl",
            lambda config, query, out_bam: called.append(query.name),
        )
        star._map_transcripts_segemehl(_config(tmp_path))
        assert called == ["trinity-denovo.jaccard.fasta"]

    def test_noop_without_transcripts(self, tmp_path, monkeypatch):
        called: list[int] = []
        monkeypatch.setattr(
            segemehl, "map_one_transcript_set_segemehl",
            lambda config, query, out_bam: called.append(1),
        )
        star._map_transcripts_segemehl(_config(tmp_path))
        assert called == []
