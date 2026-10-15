"""Optional assembly outputs left by an earlier run must not shadow fresh models.

map_transcripts maps ``<track>.jaccard.fasta`` and max_intron_split reads
``<track>.genome.defuse.gff3`` whenever they exist, so run_assembly removes them
when their step does not run this time.
"""

from __future__ import annotations

import logging

from eukan.assembly import pipeline
from eukan.assembly.bam_utils import _resolve_query
from eukan.assembly.tracks import resolve_model_source
from eukan.settings import AssemblyConfig

_JACCARD = ("trinity-denovo.jaccard.fasta", "trinity-gg.jaccard.fasta")
_DEFUSE = ("trinity-denovo.genome.defuse.gff3", "trinity-gg.genome.defuse.gff3")


def _run_assembly(tmp_path, monkeypatch, *, paired=True, **kw):
    """Run run_assembly over leftover jaccard/defuse files with the driver stubbed;
    return the leftovers still present when the driver starts and the skipped steps."""
    for name in (*_JACCARD, *_DEFUSE, "trinity-denovo.fasta",
                 "trinity-denovo.genome.stranded.gff3"):
        (tmp_path / name).write_text("x\n")
    reads = (
        {"left_reads": tmp_path / "R1.fq", "right_reads": tmp_path / "R2.fq"}
        if paired else {"single_reads": tmp_path / "R.fq"}
    )
    cfg = AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=1,
                         **reads, **kw)
    seen: dict[str, set[str]] = {}

    def fake_driver(name, steps, config, *, force_steps=None, skip=None):
        seen["present"] = {f for f in (*_JACCARD, *_DEFUSE) if (tmp_path / f).exists()}
        seen["skipped"] = {s.name for s in steps if skip is not None and skip(s)}

    monkeypatch.setattr(pipeline, "run_simple_pipeline", fake_driver)
    pipeline.run_assembly(cfg)
    return seen["present"], seen["skipped"]


def test_disabled_steps_lose_their_leftovers(tmp_path, monkeypatch, caplog):
    with caplog.at_level(logging.INFO, logger="eukan.assembly.pipeline"):
        present, skipped = _run_assembly(
            tmp_path, monkeypatch, jaccard_clip=False, defuse=False,
        )
    assert skipped == {"jaccard", "defuse"}
    assert present == set()
    removed = [r.getMessage() for r in caplog.records if r.getMessage().startswith("Removing")]
    assert len(removed) == len(_JACCARD) + len(_DEFUSE)  # one INFO line per file
    # downstream now falls back to the fresh, unclipped / un-de-fused models
    assert _resolve_query(tmp_path, "trinity-denovo.fasta").name == "trinity-denovo.fasta"
    assert resolve_model_source(tmp_path, "trinity-denovo.genome").name == (
        "trinity-denovo.genome.stranded.gff3"
    )


def test_defuse_without_uniprot_loses_its_leftovers(tmp_path, monkeypatch):
    present, skipped = _run_assembly(tmp_path, monkeypatch, defuse=True)
    assert skipped == {"defuse"}
    assert present == set(_JACCARD)


def test_single_end_jaccard_no_op_loses_its_leftovers(tmp_path, monkeypatch):
    present, _ = _run_assembly(
        tmp_path, monkeypatch, paired=False, defuse=True, uniprot_db=tmp_path / "db.dmnd",
    )
    assert present == set(_DEFUSE)


def test_enabled_steps_keep_their_outputs(tmp_path, monkeypatch):
    present, skipped = _run_assembly(
        tmp_path, monkeypatch, defuse=True, uniprot_db=tmp_path / "db.dmnd",
    )
    assert skipped == set()
    assert present == {*_JACCARD, *_DEFUSE}
