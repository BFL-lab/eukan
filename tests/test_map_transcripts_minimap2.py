"""Unit tests for minimap2 read/transcript mapping (command construction + escalation)."""

from __future__ import annotations

import json
import os
from dataclasses import asdict

import pytest

from eukan.assembly import minimap2
from eukan.assembly.bam_diagnostic import NonCanonicalSpliceCall, TransSplicingCall, Verdict
from eukan.assembly.pipeline import STEPS
from eukan.exceptions import ExternalToolError
from eukan.infra.artifacts import Artifact
from eukan.infra.manifest import ASSEMBLY, step_key
from eukan.infra.pipeline import run_simple_pipeline
from eukan.settings import AssemblyConfig


def _config(tmp_path, **kw):
    (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 100 + "\n")
    return AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=4, **kw)


def _config_reads(tmp_path, **kw):
    (tmp_path / "R1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "R2.fq").write_text("@r\nACGT\n+\nIIII\n")
    return _config(
        tmp_path, left_reads=tmp_path / "R1.fq", right_reads=tmp_path / "R2.fq", **kw
    )


def _write_verdict(work_dir, call):
    verdict = Verdict(
        TransSplicingCall(
            call="ABSENT", top_non_trivial_cluster_key="", top_non_trivial_cluster_consensus="",
            top_non_trivial_cluster_n_loci=0, top_non_trivial_cluster_n_reads=0,
            sl_bucket_pct_of_consistent=0.0,
        ),
        NonCanonicalSpliceCall(call=call, canonical_pct=50.0, top_non_canonical_dinuc=""),
    )
    (work_dir / Artifact.SOFTCLIP_DIAGNOSTIC.value).write_text(
        json.dumps({"verdict": asdict(verdict)})
    )


# ---------------------------------------------------------------------------
# Transcript -> genome mapping (splice:hq)
# ---------------------------------------------------------------------------


def _mock_transcript_mapping(monkeypatch, tmp_path):
    """Record run_piped commands; stub unmapped-extraction, sort, and finalize so
    the mapping runs without a real minimap2/samtools/pysam."""
    piped: list[tuple[list[str], list[str]]] = []
    unmapped: list[str] = []

    def fake_piped(cmd1, cmd2, **kw):
        piped.append((cmd1, cmd2))
        # emulate `samtools view -b -o <unsorted> -`: materialize the unsorted BAM
        if "-o" in cmd2:
            (tmp_path / cmd2[cmd2.index("-o") + 1]).write_text("bam")
        return ""

    def fake_unmapped(unsorted_bam, out_fasta):
        unmapped.append(out_fasta.name)
        out_fasta.write_text(">u\nACGT\n")
        return 1

    monkeypatch.setattr(minimap2, "run_piped", fake_piped)
    monkeypatch.setattr(minimap2, "run_cmd", lambda *a, **k: None)
    monkeypatch.setattr(minimap2, "_write_unmapped_fasta", fake_unmapped)
    monkeypatch.setattr(minimap2, "_bam_is_complete", lambda _p: False)
    monkeypatch.setattr(minimap2, "_finalize_transcript_diagnostics", lambda _c: None)
    return piped, unmapped


def _minimap2_cmds(piped):
    return [c1 for c1, _ in piped if c1 and c1[0] == "minimap2"]


def test_failed_transcript_sort_reports_samtools_error(tmp_path, monkeypatch):
    """A failed sort of a transcript BAM surfaces samtools' own error. The
    segemehl-era translation re-raised it as out-of-space with a hint that a
    re-run reuses the unsorted BAM and skips the mapper; minimap2 always re-maps."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "samtools").write_text(
        "#!/bin/sh\n"
        'if [ "$1" = sort ]; then\n'
        '  echo \'samtools sort: failed writing to "o.bam": Illegal seek\' >&2; exit 1\n'
        "fi\n"
    )
    (bindir / "samtools").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    real_piped = minimap2.run_piped

    def piped(cmd1, cmd2, **kw):
        if cmd1[0] != "minimap2":
            return real_piped(cmd1, cmd2, **kw)
        (tmp_path / cmd2[cmd2.index("-o") + 1]).write_text("bam")  # the unsorted BAM
        return ""

    monkeypatch.setattr(minimap2, "run_piped", piped)
    monkeypatch.setattr(minimap2, "_write_unmapped_fasta", lambda u, o: 0)
    monkeypatch.setattr(minimap2, "_bam_is_complete", lambda _p: False)
    query = tmp_path / "trinity-denovo.fasta"
    query.write_text(">t\nACGT\n")

    with pytest.raises(ExternalToolError) as exc_info:
        minimap2._map_one_transcript_set(
            _config(tmp_path), query, "trinity-denovo.genome.bam", non_canonical=False,
        )
    assert "Illegal seek" in exc_info.value.stderr_snippet
    assert "reuse" not in (exc_info.value.hint or "").lower()


def test_maps_both_trinity_tracks(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    (tmp_path / "trinity-gg.fasta").write_text(">t\nACGTACGT\n")
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)

    minimap2.map_transcripts_minimap2(_config(tmp_path))

    queries = {c[-1].rsplit("/", 1)[-1] for c in _minimap2_cmds(piped)}
    assert queries == {"trinity-denovo.fasta", "trinity-gg.fasta"}
    assert (tmp_path / "trinity-denovo.genome.bam").exists()
    assert (tmp_path / "trinity-gg.genome.bam").exists()


def test_splice_hq_flags(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)

    minimap2.map_transcripts_minimap2(_config(tmp_path, max_intron_len=5000))

    (cmd,) = _minimap2_cmds(piped)
    assert cmd[0] == "minimap2"
    assert cmd[cmd.index("-x") + 1] == "splice:hq"      # full-length cDNA, spliced
    assert cmd[cmd.index("-G") + 1] == "5000"           # native intron bound
    assert "--secondary=no" in cmd
    assert "-a" in cmd
    assert "-J" not in cmd                              # canonical by default


def test_prefers_jaccard_clipped_query(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGT\n")
    (tmp_path / "trinity-denovo.jaccard.fasta").write_text(">t\nACGTACGT\n")
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)

    minimap2.map_transcripts_minimap2(_config(tmp_path))

    (cmd,) = _minimap2_cmds(piped)
    assert cmd[-1].rsplit("/", 1)[-1] == "trinity-denovo.jaccard.fasta"
    assert (tmp_path / "trinity-denovo.genome.bam").exists()  # output keyed to the stem


def test_unmapped_captured(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGT\n")
    _, unmapped = _mock_transcript_mapping(monkeypatch, tmp_path)

    minimap2.map_transcripts_minimap2(_config(tmp_path))

    assert unmapped == ["trinity-denovo.unmapped_transcripts.fasta"]
    assert (tmp_path / "trinity-denovo.unmapped_transcripts.fasta").exists()


def test_no_assemblies_is_noop(tmp_path, monkeypatch):
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)
    minimap2.map_transcripts_minimap2(_config(tmp_path))
    assert _minimap2_cmds(piped) == []


def test_resume_skips_complete_bam(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)
    monkeypatch.setattr(minimap2, "_bam_is_complete", lambda _p: True)

    minimap2.map_transcripts_minimap2(_config(tmp_path))

    assert _minimap2_cmds(piped) == []  # a complete BAM is reused, not re-mapped


# ---------------------------------------------------------------------------
# Non-canonical selection (--non-canonical auto/force/off x verdict)
# ---------------------------------------------------------------------------


class TestUseNonCanonical:
    def test_auto_without_verdict_is_canonical(self, tmp_path):
        assert minimap2._use_non_canonical(_config(tmp_path)) is False

    def test_auto_extensive_is_non_canonical(self, tmp_path):
        _write_verdict(tmp_path, "EXTENSIVE")
        assert minimap2._use_non_canonical(_config(tmp_path)) is True

    def test_auto_moderate_is_canonical(self, tmp_path):
        _write_verdict(tmp_path, "MODERATE")
        assert minimap2._use_non_canonical(_config(tmp_path)) is False

    def test_force_is_non_canonical_without_verdict(self, tmp_path):
        assert minimap2._use_non_canonical(_config(tmp_path, non_canonical="force")) is True

    def test_off_ignores_extensive_verdict(self, tmp_path):
        _write_verdict(tmp_path, "EXTENSIVE")
        assert minimap2._use_non_canonical(_config(tmp_path, non_canonical="off")) is False


def test_transcript_mapping_applies_nc_flags_when_extensive(tmp_path, monkeypatch):
    (tmp_path / "trinity-denovo.fasta").write_text(">t\nACGTACGT\n")
    _write_verdict(tmp_path, "EXTENSIVE")
    piped, _ = _mock_transcript_mapping(monkeypatch, tmp_path)

    minimap2.map_transcripts_minimap2(_config(tmp_path))

    (cmd,) = _minimap2_cmds(piped)
    assert cmd[cmd.index("-J") + 1] == "0"
    assert cmd[cmd.index("-C") + 1] == "3"
    assert "--splice-flank=no" in cmd


# ---------------------------------------------------------------------------
# Read mapping escalation (map_reads_minimap2 canonical -> non-canonical)
# ---------------------------------------------------------------------------


class TestReadMappingEscalation:
    """``map_reads_minimap2`` re-maps with the non-canonical flags on an EXTENSIVE
    verdict in ``auto`` mode; ``force`` maps non-canonically from the start; ``off``
    never escalates."""

    def _patch(self, monkeypatch, verdict=None):
        calls: list[bool] = []
        monkeypatch.setattr(
            minimap2, "_map_reads_once",
            lambda config, *, non_canonical: calls.append(non_canonical),
        )
        monkeypatch.setattr(minimap2, "_log_mapping_rate", lambda *a, **k: None)

        def fake_emit(config):
            if verdict is not None:
                _write_verdict(config.work_dir, verdict)

        monkeypatch.setattr(minimap2, "_emit_sj_and_hints", fake_emit)
        return calls

    def test_escalates_on_extensive(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="EXTENSIVE")
        minimap2.map_reads_minimap2(_config_reads(tmp_path))
        assert calls == [False, True]  # canonical, then non-canonical re-map

    def test_no_escalation_on_moderate(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="MODERATE")
        minimap2.map_reads_minimap2(_config_reads(tmp_path))
        assert calls == [False]

    def test_no_escalation_without_verdict(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict=None)
        minimap2.map_reads_minimap2(_config_reads(tmp_path))
        assert calls == [False]

    def test_force_maps_non_canonical_from_start(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="EXTENSIVE")
        minimap2.map_reads_minimap2(_config_reads(tmp_path, non_canonical="force"))
        assert calls == [True]  # forced from the start, no second pass

    def test_off_never_escalates(self, tmp_path, monkeypatch):
        calls = self._patch(monkeypatch, verdict="EXTENSIVE")
        minimap2.map_reads_minimap2(_config_reads(tmp_path, non_canonical="off"))
        assert calls == [False]


def test_rerun_of_read_mapping_rebuilds_softclip_verdict(tmp_path, monkeypatch):
    """A re-run of the completed minimap2 step clears the soft-clip summary before
    mapping, so the verdict is rebuilt from the new BAM; the non-canonical re-map
    inside that run still reuses the verdict its first pass wrote."""
    monkeypatch.setattr(
        minimap2, "_map_reads_once",
        lambda config, *, non_canonical: (config.work_dir / minimap2._BAM).write_text("bam"),
    )
    monkeypatch.setattr(minimap2, "_log_mapping_rate", lambda *a, **k: None)
    summary_present: list[bool] = []

    def fake_emit(config):
        # run_softclip_diagnostic walks the BAM only when no summary exists
        summary = config.work_dir / Artifact.SOFTCLIP_DIAGNOSTIC.value
        summary_present.append(summary.exists())
        if not summary.exists():
            _write_verdict(config.work_dir, "EXTENSIVE")

    monkeypatch.setattr(minimap2, "_emit_sj_and_hints", fake_emit)
    cfg = _config_reads(tmp_path)
    steps = [s for s in STEPS if s.name == "minimap2"]
    run_simple_pipeline(ASSEMBLY, steps, cfg)
    polya = tmp_path / Artifact.POLYA_DIAGNOSTIC.value
    polya.write_text(json.dumps({"transcripts": {"trinity-denovo": {}}}))
    summary_present.clear()

    run_simple_pipeline(ASSEMBLY, steps, cfg, force_steps=[step_key(ASSEMBLY, "minimap2")])

    assert summary_present == [False, True]  # cleared before mapping, reused by the re-map
    # the transcript sections map_transcripts wrote are not the read mapping's to drop
    assert "transcripts" in json.loads(polya.read_text())


def test_read_mapping_splice_sr_flags(tmp_path, monkeypatch):
    piped: list[tuple[list[str], list[str]]] = []
    monkeypatch.setattr(
        minimap2, "run_piped", lambda c1, c2, **k: piped.append((c1, c2)) or ""
    )
    monkeypatch.setattr(minimap2, "run_cmd", lambda *a, **k: None)
    monkeypatch.setattr(minimap2, "_log_mapping_rate", lambda *a, **k: None)
    monkeypatch.setattr(minimap2, "_emit_sj_and_hints", lambda *a, **k: None)

    minimap2.map_reads_minimap2(_config_reads(tmp_path, max_intron_len=5000))

    (cmd, sort_cmd) = next((c1, c2) for c1, c2 in piped if c1[0] == "minimap2")
    assert cmd[cmd.index("-x") + 1] == "splice:sr"      # short-read spliced preset
    assert cmd[cmd.index("-G") + 1] == "5000"
    assert "--secondary=no" in cmd
    assert "-J" not in cmd                              # canonical by default (auto, no verdict)
    assert cmd[-2:] == [str(tmp_path / "R1.fq"), str(tmp_path / "R2.fq")]  # reads positional
    assert sort_cmd[:2] == ["samtools", "sort"]
