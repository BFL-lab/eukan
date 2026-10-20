"""Tests for eukan.annotation.orf — ORF finding and coordinate mapping."""

from eukan.annotation.orf import find_longest_orfs


def _orf(seq, genetic_code=1, **kw):
    return find_longest_orfs([("t1", "+", seq)], genetic_code=genetic_code, **kw).get("t1")


# ---------------------------------------------------------------------------
# find_longest_orfs
# ---------------------------------------------------------------------------


class TestFindLongestOrfs:
    def test_basic_orf(self):
        """1-based start of the start codon to the last base of the stop codon."""
        seq = "CC" + "ATG" + "GCC" * 40 + "TAA" + "CC"
        assert _orf(seq) == (3, 128)

    def test_first_start_after_the_previous_stop(self):
        """An in-frame start downstream of the first one does not shorten the ORF."""
        seq = "ATG" + "GCC" * 10 + "ATG" + "GCC" * 30 + "TAG"
        assert _orf(seq) == (1, len(seq))

    def test_upstream_orf_wins_a_tie_in_frame(self):
        orf = "ATG" + "GCC" * 40 + "TAA"
        assert _orf(orf + orf, orf_len_frac=0.4) == (1, len(orf))

    def test_cross_frame_tie_goes_to_the_frame_whose_first_start_comes_first(self):
        """Frame 2 has a start at position 3 with no ORF behind it, ahead of
        frame 0's first start; both frames hold an equally long ORF."""
        orf = "ATG" + "GCC" * 40 + "TAA"
        seq = "GGATGTAA" + "C" * 22 + orf + "GG" + orf + "CC"
        assert _orf(seq, orf_len_frac=0.3) == (len(seq) - 2 - len(orf) + 1, len(seq) - 2)

    def test_min_length_is_strict(self):
        assert _orf("ATG" + "GCC" * 29 + "TAA", orf_len_frac=0.0) is None  # 90 nt
        assert _orf("ATG" + "GCC" * 30 + "TAA", orf_len_frac=0.0) == (1, 96)  # 93 nt

    def test_fraction_filter_is_strict(self):
        seq = "ATG" + "GCC" * 49 + "TAA" + "C" * 147  # 150 / 300 nt
        assert _orf(seq) is None
        assert _orf(seq[:-3]) == (1, 153)

    def test_no_orf(self):
        assert _orf("") is None
        assert _orf("C" * 300) is None
        assert _orf("ATG" + "GCC" * 100) is None  # no stop

    def test_genetic_code_6_taa_not_stop(self):
        """Under code 6 (Ciliate), TAA/TAG encode glutamine; TGA is the only stop."""
        seq = "ATG" + "GCC" * 40 + "TAA" + "GCC" * 40 + "TGA"
        assert _orf(seq, orf_len_frac=0.0) == (1, 126)
        assert _orf(seq, genetic_code=6, orf_len_frac=0.0) == (1, len(seq))

    def test_longest_orf_starts_at_overlapped_start_codon(self):
        """Under code 11 the longest ORF starts at TTG@2 (frame 1), which
        overlaps ATT@1 (frame 0). With that start lost, the scan fell back
        to the in-frame ATG@20 and reported a shorter ORF."""
        seq = "ATTG" + "AAA" * 5 + "ATG" + "AAA" * 35 + "TAA"
        assert _orf(seq, genetic_code=11) == (2, len(seq))

    def test_overlapping_stop_codons_all_found(self):
        """Code 2 stops overlap in TAGA: the in-frame AGA, one base after an
        out-of-frame TAG, ends the ORF; a resumed scan read through it."""
        seq = "CCC" + "ATG" + "GCC" * 39 + "GCT" + "AGA" + "GCC" * 10 + "TAA"
        assert _orf(seq, genetic_code=2, orf_len_frac=0.0) == (4, 129)

    def test_one_entry_per_transcript_with_an_orf(self):
        orf = "ATG" + "GCC" * 40 + "TAA"
        seqs = [("a", "+", orf), ("b", "-", "C" * 200), ("c", "-", "G" + orf)]
        assert find_longest_orfs(seqs) == {"a": (1, 126), "c": (2, 127)}
