"""Tests for eukan.annotation.orf.create_transcriptome_orf_db."""

from eukan.annotation.orf import create_transcriptome_orf_db


class TestCreateTranscriptomeOrfDb:
    def test_exon_ids_repeated_across_transcripts(self, tmp_path):
        """Exon IDs that restart per transcript (exon1, exon2 in each) no
        longer stop the parse with gffutils' "Duplicate ID" error."""
        orf = "ATG" + "GCC" * 40 + "TAA"  # 126 nt
        locus = orf[:60] + "GT" + "C" * 46 + "AG" + orf[60:]
        pad = "C" * 100
        (tmp_path / "genome.fa").write_text(f">chr1\n{pad}{locus}{pad}{locus}{pad}\n")
        rows = []
        for tid, s in (("t1", 101), ("t2", 101 + len(locus) + 100)):
            rows.append(f"chr1\tx\texon\t{s}\t{s + 59}\t.\t+\t.\tID=exon1;Parent={tid}")
            rows.append(f"chr1\tx\texon\t{s + 110}\t{s + len(locus) - 1}\t.\t+\t.\tID=exon2;Parent={tid}")
        (tmp_path / "transcripts.gff3").write_text("\n".join(rows) + "\n")

        db = create_transcriptome_orf_db(tmp_path / "transcripts.gff3", tmp_path / "genome.fa")

        assert [g.id for g in db.features_of_type("gene", order_by="start")] == ["t1", "t2"]
        cds = [(f.start, f.end) for f in db.features_of_type("CDS", order_by="start")]
        assert cds == [(101, 160), (211, 276), (377, 436), (487, 552)]
