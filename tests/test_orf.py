"""Tests for eukan.annotation.orf — ORF finding and coordinate mapping."""

import pandas as pd

from eukan.annotation.orf import fetch_longest_orf, find_starts_stops

# ---------------------------------------------------------------------------
# find_starts_stops
# ---------------------------------------------------------------------------


class TestFindStartsStops:
    def test_basic_start_stop(self):
        """Should find ATG starts and TAA/TAG/TGA stops."""
        seqs = [("seq1", "+", "AAATGCCCCTAAGGG")]
        #                       ATG at pos 3, TAA at pos 10 (0-based: 2,9)
        df = find_starts_stops(seqs, genetic_code=1)

        starts = df[df["codon"] == "start"]
        stops = df[df["codon"] == "stop"]

        assert len(starts) >= 1
        assert len(stops) >= 1
        # ATG at position 3 (1-based)
        assert 3 in starts["pos"].values

    def test_empty_sequence(self):
        """Empty sequence should return empty DataFrame."""
        df = find_starts_stops([("seq1", "+", "")], genetic_code=1)
        assert len(df) == 0

    def test_no_starts(self):
        """Sequence with no start codons."""
        df = find_starts_stops([("seq1", "+", "CCCCCCCCCC")])
        starts = df[df["codon"] == "start"]
        assert len(starts) == 0

    def test_genetic_code_6_taa_not_stop(self):
        """Under code 6 (Ciliate), TAA/TAG encode glutamine; TGA is the only stop."""
        seqs = [("seq1", "+", "AAATGCCCCTAAGGG")]
        #                       ATG at pos 3, TAA at pos 10

        # Code 1: TAA is a stop codon
        df_code1 = find_starts_stops(seqs, genetic_code=1)
        stops_1 = df_code1[df_code1["codon"] == "stop"]
        assert any(stops_1["pos"] == 10)

        # Code 6: TAA is NOT a stop (only TGA is)
        df_code6 = find_starts_stops(seqs, genetic_code=6)
        stops_6 = df_code6[df_code6["codon"] == "stop"]
        assert 10 not in stops_6["pos"].values

    def test_overlapping_start_codons_all_found(self):
        """Code 11 starts overlap: ATT@3 and TTG@4 share bases in CCATTGCC.
        A non-overlapping scan resumed after ATT and lost the TTG."""
        df = find_starts_stops([("seq1", "+", "CCATTGCC")], genetic_code=11)
        starts = df[df["codon"] == "start"]
        assert sorted(starts["pos"]) == [3, 4]

    def test_overlapping_stop_codons_all_found(self):
        """Code 2 stops overlap: TAG@2 and AGA@3 share bases in CTAGAC."""
        df = find_starts_stops([("seq1", "+", "CTAGAC")], genetic_code=2)
        stops = df[df["codon"] == "stop"]
        assert sorted(stops["pos"]) == [2, 3]


# ---------------------------------------------------------------------------
# fetch_longest_orf
# ---------------------------------------------------------------------------


class TestFetchLongestORF:
    def test_basic_orf(self):
        """Should find the longest in-frame start-stop pair."""
        # Create a simple case: one start at pos 1, one stop at pos 100
        data = pd.DataFrame(
            {
                "id": ["s1"] * 4,
                "seqlen": [200] * 4,
                "pos": [1, 10, 100, 150],
                "strand": ["+"] * 4,
                "codon": ["start", "start", "stop", "stop"],
                "frame": [0, 0, 0, 0],
            }
        )
        result = fetch_longest_orf(data, min_orf_len=50, orf_len_frac=0.1)
        assert len(result) == 1
        assert result.iloc[0]["start"] == 1
        # stop adjusted +2 for full codon
        assert result.iloc[0]["stop"] == 102

    def test_min_length_filter(self):
        """ORFs below min_orf_len should be filtered out."""
        data = pd.DataFrame(
            {
                "id": ["s1", "s1"],
                "seqlen": [200, 200],
                "pos": [1, 20],
                "strand": ["+", "+"],
                "codon": ["start", "stop"],
                "frame": [0, 0],
            }
        )
        result = fetch_longest_orf(data, min_orf_len=50, orf_len_frac=0.0)
        assert len(result) == 0  # ORF is only 19 nt

    def test_empty_input(self):
        """Empty input should return empty DataFrame."""
        data = pd.DataFrame(
            columns=["id", "seqlen", "pos", "strand", "codon", "frame"]
        )
        result = fetch_longest_orf(data)
        assert len(result) == 0

    def test_fraction_filter(self):
        """ORFs below orf_len_frac should be filtered out."""
        data = pd.DataFrame(
            {
                "id": ["s1", "s1"],
                "seqlen": [1000, 1000],
                "pos": [1, 100],
                "strand": ["+", "+"],
                "codon": ["start", "stop"],
                "frame": [0, 0],
            }
        )
        # ORF is 99 nt / 1000 nt = 0.099, below default 0.50
        result = fetch_longest_orf(data, min_orf_len=10, orf_len_frac=0.50)
        assert len(result) == 0

    def test_longest_orf_starts_at_overlapped_start_codon(self):
        """Under code 11 the longest ORF starts at TTG@2 (frame 1), which
        overlaps ATT@1 (frame 0). With that start lost, the scan fell back
        to the in-frame ATG@20 and reported a shorter ORF."""
        seq = "ATTG" + "AAA" * 5 + "ATG" + "AAA" * 35 + "TAA"
        result = fetch_longest_orf(find_starts_stops([("t1", "+", seq)], genetic_code=11))
        assert len(result) == 1
        orf = result.iloc[0]
        assert (orf["start"], orf["stop"], orf["frame"]) == (2, len(seq), 1)
