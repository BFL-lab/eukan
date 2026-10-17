"""Tests for the eukan.repeats package — sort, GFF→hints, config."""

from Bio import SeqIO

from eukan.repeats.masker import gff_to_hints
from eukan.repeats.modeler import _salvage_partial_library, sort_and_uppercase
from eukan.repeats.pipeline import force_steps_from_run_flags
from eukan.settings import RepeatsConfig


class TestSortAndUppercase:
    def test_sorts_by_length_descending(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">short\nacgt\n>long\nacgtacgtacgt\n>medium\nacgtacgt\n")
        dst = tmp_path / "out.fa"

        sort_and_uppercase(src, dst)

        records = list(SeqIO.parse(str(dst), "fasta"))
        assert [r.id for r in records] == ["long", "medium", "short"]

    def test_uppercases_sequence(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">a\nacgtacgt\n")
        dst = tmp_path / "out.fa"

        sort_and_uppercase(src, dst)

        records = list(SeqIO.parse(str(dst), "fasta"))
        assert str(records[0].seq) == "ACGTACGT"


class TestGffToHints:
    def test_rewrites_columns(self, tmp_path):
        repeat_gff = tmp_path / "in.gff"
        repeat_gff.write_text(
            "##gff-version 3\n"
            "\n"
            "# a comment\n"
            "chr1\tRepeatMasker\tdispersed_repeat\t100\t200\t12.3\t+\t.\t"
            'Target "Motif:foo" 1 100\n'
            "chr1\tRepeatMasker\tsimilarity\t300\t400\t1.2\t-\t.\tTarget=bar\n"
        )
        hints = tmp_path / "hints.gff"

        gff_to_hints(repeat_gff, hints)

        rows = [
            line.split("\t") for line in hints.read_text().splitlines() if line.strip()
        ]
        assert len(rows) == 2
        for row in rows:
            assert row[2] == "nonexonpart"
            assert row[8] == "src=RM"
        assert rows[0][:2] == ["chr1", "RepeatMasker"]
        assert rows[0][3:5] == ["100", "200"]
        assert rows[1][3:5] == ["300", "400"]

    def test_skips_short_lines(self, tmp_path):
        repeat_gff = tmp_path / "in.gff"
        repeat_gff.write_text(
            "chr1\tRepeatMasker\tdispersed_repeat\t100\t200\t12.3\t+\t.\tT=foo\n"
            "broken\tline\n"
        )
        hints = tmp_path / "hints.gff"

        gff_to_hints(repeat_gff, hints)

        assert len(hints.read_text().splitlines()) == 1


class TestRepeatsConfig:
    def test_defaults(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = RepeatsConfig(genome=genome)

        assert config.engine == "rmblast"
        assert config.lib is None
        assert config.name == "genome"
        assert config.num_cpu >= 1

    def test_lib_override(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        lib = tmp_path / "families.fa"
        lib.touch()
        config = RepeatsConfig(genome=genome, lib=lib)
        assert config.lib == lib


class TestSalvagePartialLibrary:
    """``_salvage_partial_library`` recovers usable FASTA from crashed runs."""

    def test_no_rm_dir_returns_none(self, tmp_path):
        assert _salvage_partial_library(tmp_path) is None

    def test_prefers_cumulative_consensi(self, tmp_path):
        rm = tmp_path / "RM_123.SomeStamp"
        (rm / "round-1").mkdir(parents=True)
        (rm / "consensi.fa").write_text(">cumulative\nACGT\n")
        (rm / "round-1" / "consensi.fa").write_text(">round-1-only\nAAAA\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == rm / "consensi.fa"

    def test_falls_back_to_latest_round_classified(self, tmp_path):
        rm = tmp_path / "RM_123.SomeStamp"
        # Cumulative is empty → skip; round-2 has classified, round-1 unclassified
        (rm / "round-1").mkdir(parents=True)
        (rm / "round-2").mkdir(parents=True)
        (rm / "consensi.fa").write_text("")  # empty -> skipped
        (rm / "round-1" / "consensi.fa").write_text(">r1\nAAAA\n")
        (rm / "round-2" / "consensi.fa.classified").write_text(">r2#LINE/L1\nCCCC\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == rm / "round-2" / "consensi.fa.classified"

    def test_round_consensi_unclassified_when_classified_missing(self, tmp_path):
        rm = tmp_path / "RM_123.SomeStamp"
        (rm / "round-1").mkdir(parents=True)
        (rm / "round-1" / "consensi.fa").write_text(">r1\nAAAA\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == rm / "round-1" / "consensi.fa"

    def test_refined_cons_is_last_resort(self, tmp_path):
        rm = tmp_path / "RM_123.SomeStamp"
        (rm / "round-1").mkdir(parents=True)
        (rm / "round-1" / "refined-cons.fa").write_text(">refiner\nGGGG\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == rm / "round-1" / "refined-cons.fa"

    def test_picks_latest_rm_dir(self, tmp_path):
        old = tmp_path / "RM_111.A"
        new = tmp_path / "RM_999.Z"
        for d in (old, new):
            d.mkdir()
        old.joinpath("consensi.fa").write_text(">old\nACGT\n")
        new.joinpath("consensi.fa").write_text(">new\nACGT\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == new / "consensi.fa"

    def test_skips_empty_files(self, tmp_path):
        rm = tmp_path / "RM_123.A"
        (rm / "round-1").mkdir(parents=True)
        (rm / "consensi.fa").write_text("")
        (rm / "round-1" / "consensi.fa.classified").write_text("")
        (rm / "round-1" / "consensi.fa").write_text(">x\nA\n")

        salvaged = _salvage_partial_library(tmp_path)
        assert salvaged == rm / "round-1" / "consensi.fa"


class TestForceStepsFromRunFlags:
    """Mirrors annotation/assembly: ``--run-X`` and ``--force`` translate
    to a list of full ``repeats/<step>`` manifest keys to force.
    """

    def test_no_flags_returns_empty(self):
        assert force_steps_from_run_flags() == []

    def test_explicit_run_flag_forces_that_step(self):
        assert force_steps_from_run_flags(run_masker=True) == ["repeats/masker"]

    def test_run_modeler_also_forces_masker(self):
        # The masker reads the library the modeler rebuilds.
        assert force_steps_from_run_flags(run_modeler=True) == [
            "repeats/modeler", "repeats/masker",
        ]

    def test_force_alone_returns_all_keys(self):
        assert force_steps_from_run_flags(force=True) == [
            "repeats/modeler", "repeats/masker",
        ]

    def test_returned_keys_are_prefixed(self):
        result = force_steps_from_run_flags(run_modeler=True, run_masker=True)
        assert all(k.startswith("repeats/") for k in result)


class TestRunRepeatsResume:
    """Through run_repeats: ``--run-modeler`` re-runs the masker on the rebuilt
    library, the masker re-runs when its library or the genome changed, and a
    changed genome rebuilds the library."""

    @staticmethod
    def _stub(monkeypatch, ran):
        from eukan.repeats import pipeline

        def fake_modeler(config):
            ran.append("modeler")
            out = pipeline._modeler_output(config)
            out.write_text(f">family{len(ran)}\nACGT\n")  # RepeatModeler samples at random
            return out

        def fake_masker(config, families):
            ran.append("masker")
            masked = config.work_dir / f"{config.name}.masked.fasta"
            masked.write_text(">chr1\nacgt\n")
            return masked, config.work_dir / "hints_repeatmask.gff"

        monkeypatch.setattr(pipeline, "run_modeler", fake_modeler)
        monkeypatch.setattr(pipeline, "run_masker", fake_masker)
        return pipeline.run_repeats

    def test_run_modeler_reruns_masker(self, tmp_path, monkeypatch):
        ran: list[str] = []
        run_repeats = self._stub(monkeypatch, ran)
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        config = RepeatsConfig(genome=genome, work_dir=tmp_path)
        run_repeats(config)
        ran.clear()
        run_repeats(config, force_steps=force_steps_from_run_flags(run_modeler=True))
        assert ran == ["modeler", "masker"]

    def test_changed_library_or_genome_reruns_masker(self, tmp_path, monkeypatch):
        ran: list[str] = []
        run_repeats = self._stub(monkeypatch, ran)
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        lib = tmp_path / "families.fa"
        lib.write_text(">a\nACGT\n")
        config = RepeatsConfig(genome=genome, work_dir=tmp_path, lib=lib)
        run_repeats(config)
        run_repeats(config)
        assert ran == ["masker"]  # unchanged: skipped
        lib.write_text(">a\nACGTT\n")
        run_repeats(config)
        genome.write_text(">chr1\nACGTA\n")
        run_repeats(config)
        assert ran == ["masker"] * 3

    def test_changed_genome_rebuilds_the_library(self, tmp_path, monkeypatch):
        ran: list[str] = []
        run_repeats = self._stub(monkeypatch, ran)
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        config = RepeatsConfig(genome=genome, work_dir=tmp_path)
        run_repeats(config)
        run_repeats(config)
        assert ran == ["modeler", "masker"]  # unchanged: skipped
        genome.write_text(">chr1\nACGTA\n")
        run_repeats(config)
        assert ran == ["modeler", "masker"] * 2
