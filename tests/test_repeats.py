"""Tests for the eukan.repeats package — sort, GFF→hints, config."""

from Bio import SeqIO

from eukan.repeats.masker import gff_to_hints
from eukan.repeats.modeler import sort_and_uppercase
from eukan.repeats.orchestrator import steps_and_force_from_run_flags
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


class TestStepsAndForce:
    def test_default_runs_both(self):
        steps, force = steps_and_force_from_run_flags()
        assert steps == ["modeler", "masker"]
        assert force is False

    def test_explicit_run_flag_implies_force(self):
        steps, force = steps_and_force_from_run_flags(run_masker=True)
        assert steps == ["masker"]
        assert force is True

    def test_force_alone(self):
        steps, force = steps_and_force_from_run_flags(force=True)
        assert steps == ["modeler", "masker"]
        assert force is True
