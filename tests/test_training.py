"""Tests for eukan.annotation.training.build_training_set."""

from types import SimpleNamespace

from eukan.annotation import training


def _model(gene_id, start, end):
    return (
        f"chr1\tpred\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}\n"
        f"chr1\tpred\tmRNA\t{start}\t{end}\t.\t+\t.\tID=m_{gene_id};Parent={gene_id}\n"
        f"chr1\tpred\tCDS\t{start}\t{end}\t.\t+\t0\tID=c_{gene_id};Parent=m_{gene_id}\n"
    )


def test_flank_is_half_the_first_source_mean_gene_length(tmp_path, monkeypatch):
    gm = tmp_path / "genemark.gff3"
    gm.write_text(_model("g1", 101, 400) + _model("g2", 1001, 1900) + _model("g3", 3001, 3301))
    prot = tmp_path / "prot.gff3"
    prot.write_text(_model("p1", 101, 400) + _model("p3", 3001, 3301))
    calls = []
    monkeypatch.setattr(training, "run_cmd", lambda cmd, **kw: calls.append(cmd))
    config = SimpleNamespace(genome=tmp_path / "genome.fa")

    training.build_training_set(config, (gm, prot), tmp_path)

    # mean(end - start) = (299 + 899 + 300) / 3 = 499.33 -> flank int(249.67)
    assert calls == [
        ["gff2gbSmallDNA.pl", "training_set.gff3", str(tmp_path / "genome.fa"), "249", "genbank.gb"],
        ["randomSplit.pl", "genbank.gb", "1"],
    ]
    assert (tmp_path / "training_set.gff3").read_text().count("\tgene\t") == 2
