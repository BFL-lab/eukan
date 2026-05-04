"""Tests for eukan.settings — pydantic-settings configuration."""

from pathlib import Path

from eukan.settings import AssemblyConfig, FunctionalConfig, Kingdom, PipelineConfig


class TestPipelineConfig:
    def test_defaults(self, tmp_path):
        """Should fill in defaults for optional fields."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = PipelineConfig(genome=genome, proteins=[genome])

        assert config.genetic_code == "11"
        assert config.weights == [2, 1, 3]
        assert config.kingdom is None
        assert config.num_cpu >= 1
        assert config.name == "genome"
        assert len(config.shortname) == 5

    def test_name_derived_from_genome(self, tmp_path):
        """Name should default to genome stem."""
        genome = tmp_path / "my_organism.fasta"
        genome.touch()
        config = PipelineConfig(genome=genome, proteins=[genome])
        assert config.name == "my_organism"

    def test_explicit_name_preserved(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = PipelineConfig(genome=genome, proteins=[genome], name="custom")
        assert config.name == "custom"

    def test_kingdom_enum(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = PipelineConfig(genome=genome, proteins=[genome], kingdom="fungus")
        assert config.kingdom == Kingdom.fungus
        assert config.is_fungus is True
        assert config.is_protist is False

    def test_has_transcripts(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()

        config = PipelineConfig(genome=genome, proteins=[genome])
        assert config.has_transcripts is False

        config2 = PipelineConfig(
            genome=genome,
            proteins=[genome],
            transcripts_fasta=genome,
            transcripts_gff=genome,
            rnaseq_hints=genome,
        )
        assert config2.has_transcripts is True

    def test_auto_discover_assembly_outputs(self, tmp_path):
        """Should auto-discover assembly outputs in work_dir."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        (tmp_path / "nr_transcripts.fasta").write_text(">seq\nACGT\n")
        (tmp_path / "nr_transcripts.gff3").write_text("##gff-version 3\n")
        (tmp_path / "hints_rnaseq.gff").write_text("# hints\n")

        config = PipelineConfig(genome=genome, proteins=[genome], work_dir=tmp_path)
        assert config.has_transcripts is True
        assert config.transcripts_fasta == tmp_path / "nr_transcripts.fasta"
        assert config.transcripts_gff == tmp_path / "nr_transcripts.gff3"
        assert config.rnaseq_hints == tmp_path / "hints_rnaseq.gff"

    def test_no_auto_discover_when_missing(self, tmp_path):
        """Should not set transcripts when assembly outputs are absent."""
        genome = tmp_path / "genome.fa"
        genome.touch()

        config = PipelineConfig(genome=genome, proteins=[genome], work_dir=tmp_path)
        assert config.has_transcripts is False
        assert config.transcripts_fasta is None

    def test_no_auto_discover_when_explicit(self, tmp_path):
        """Explicit transcript paths should not be overwritten by auto-discovery."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        custom = tmp_path / "custom.fasta"
        custom.touch()
        # Create assembly outputs that would be auto-discovered
        (tmp_path / "nr_transcripts.fasta").write_text(">seq\nACGT\n")
        (tmp_path / "nr_transcripts.gff3").write_text("##gff-version 3\n")
        (tmp_path / "hints_rnaseq.gff").write_text("# hints\n")

        config = PipelineConfig(
            genome=genome, proteins=[genome], work_dir=tmp_path,
            transcripts_fasta=custom, transcripts_gff=custom, rnaseq_hints=custom,
        )
        assert config.transcripts_fasta == custom

    def test_partial_assembly_outputs_warns(self, tmp_path, caplog):
        """Should warn when only some assembly outputs exist."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        (tmp_path / "nr_transcripts.fasta").write_text(">seq\nACGT\n")
        # Missing nr_transcripts.gff3 and hints_rnaseq.gff

        import logging
        with caplog.at_level(logging.WARNING, logger="eukan.settings"):
            config = PipelineConfig(genome=genome, proteins=[genome], work_dir=tmp_path)
        assert config.has_transcripts is False
        assert "Partial assembly outputs" in caplog.text

    def test_env_var_override(self, tmp_path, monkeypatch):
        """EUKAN_ env vars should override defaults."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        monkeypatch.setenv("EUKAN_GENETIC_CODE", "6")
        monkeypatch.setenv("EUKAN_NUM_CPU", "4")

        config = PipelineConfig(genome=genome, proteins=[genome])
        assert config.genetic_code == "6"
        assert config.num_cpu == 4

    def test_constructor_overrides_env(self, tmp_path, monkeypatch):
        """Explicit constructor args should beat env vars."""
        genome = tmp_path / "genome.fa"
        genome.touch()
        monkeypatch.setenv("EUKAN_GENETIC_CODE", "6")

        config = PipelineConfig(genome=genome, proteins=[genome], genetic_code="11")
        assert config.genetic_code == "11"


class TestAssemblyConfig:
    def test_defaults(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = AssemblyConfig(genome=genome)

        assert config.min_intron_len == 20
        assert config.max_intron_len == 5000
        assert config.align_mode == "Local"
        assert config.name == "genome"

    def test_reads_args_star(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        left = tmp_path / "left.fq"
        right = tmp_path / "right.fq"
        config = AssemblyConfig(genome=genome, left_reads=left, right_reads=right)
        assert config.reads_args_star == [str(left), str(right)]

    def test_reads_args_single(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        single = tmp_path / "reads.fq"
        config = AssemblyConfig(genome=genome, single_reads=single)
        assert config.reads_args_star == [str(single)]

    def test_memory_gb_default_uses_meminfo(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = AssemblyConfig(genome=genome)
        # Default factory always returns at least the 4 GiB floor.
        assert isinstance(config.memory_gb, int)
        assert config.memory_gb >= 4

    def test_memory_gb_explicit_override(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.touch()
        config = AssemblyConfig(genome=genome, memory_gb=12)
        assert config.memory_gb == 12


class TestDefaultTrinityMemoryGb:
    """Direct tests for ``_default_trinity_memory_gb`` against fixture meminfo."""

    def test_uses_mem_available_at_60_percent(self, tmp_path):
        from eukan.settings import _default_trinity_memory_gb

        # 32 GiB total, 16 GiB available -> 0.6 * 16 = 9.6 -> 9
        meminfo = tmp_path / "meminfo"
        meminfo.write_text(
            "MemTotal:       33554432 kB\n"
            "MemFree:         8388608 kB\n"
            "MemAvailable:   16777216 kB\n"
        )
        assert _default_trinity_memory_gb(str(meminfo)) == 9

    def test_falls_back_to_half_total_without_mem_available(self, tmp_path):
        from eukan.settings import _default_trinity_memory_gb

        # No MemAvailable line. 32 GiB total -> 16 GiB.
        meminfo = tmp_path / "meminfo"
        meminfo.write_text(
            "MemTotal:       33554432 kB\n"
            "MemFree:         8388608 kB\n"
        )
        assert _default_trinity_memory_gb(str(meminfo)) == 16

    def test_floor_at_4_gib_on_low_memory(self, tmp_path):
        from eukan.settings import _default_trinity_memory_gb

        # 4 GiB total, 1 GiB available -> 0.6 * 1 = 0.6 -> floored to 4.
        meminfo = tmp_path / "meminfo"
        meminfo.write_text(
            "MemTotal:        4194304 kB\n"
            "MemAvailable:    1048576 kB\n"
        )
        assert _default_trinity_memory_gb(str(meminfo)) == 4

    def test_returns_4_when_meminfo_missing(self, tmp_path):
        from eukan.settings import _default_trinity_memory_gb

        assert _default_trinity_memory_gb(str(tmp_path / "no-such-file")) == 4

    def test_returns_4_on_malformed_meminfo(self, tmp_path):
        from eukan.settings import _default_trinity_memory_gb

        meminfo = tmp_path / "meminfo"
        meminfo.write_text("MemAvailable:   not-a-number kB\n")
        assert _default_trinity_memory_gb(str(meminfo)) == 4


class TestFunctionalConfig:
    def test_default_db_paths(self, tmp_path):
        proteins = tmp_path / "proteins.faa"
        proteins.touch()
        config = FunctionalConfig(proteins=proteins)

        assert config.uniprot_db == Path("databases/uniprot_sprot.faa")
        assert config.pfam_db == Path("databases/Pfam-A.hmm")

    def test_env_var_override(self, tmp_path, monkeypatch):
        proteins = tmp_path / "proteins.faa"
        proteins.touch()
        monkeypatch.setenv("EUKAN_FUNC_UNIPROT_DB", "/custom/uniprot.faa")

        config = FunctionalConfig(proteins=proteins)
        assert config.uniprot_db == Path("/custom/uniprot.faa")
