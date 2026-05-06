"""Tests for eukan.annotation.consensus — PASA output resolution."""

from __future__ import annotations

import os
import time

from eukan.annotation.consensus import _resolve_consensus_path


class TestResolveConsensusPath:
    """``_resolve_consensus_path`` must return the *latest* PASA iteration.

    Regression: previously sorted file paths ascending and returned
    ``[0]`` — that's the *oldest* iteration. PASA suffixes can be
    variable-width process IDs, so even a "right index" lexicographic
    pick is fragile; we now sort by mtime instead.
    """

    def test_no_pasa_outputs_returns_evm_consensus(self, tmp_path):
        result = _resolve_consensus_path(tmp_path)
        assert result == tmp_path / "consensus_models.gff3"

    def test_picks_most_recently_modified_pasa_file(self, tmp_path):
        # Create three PASA iterations with monotonically-increasing mtime.
        old = tmp_path / "db.gene_structures_post_PASA_updates.99999.gff3"
        mid = tmp_path / "db.gene_structures_post_PASA_updates.55555.gff3"
        new = tmp_path / "db.gene_structures_post_PASA_updates.11111.gff3"
        for path in (old, mid, new):
            path.write_text("# placeholder\n")

        # Force ordered mtimes so the test doesn't depend on filesystem
        # tie-break order: old < mid < new.
        now = time.time()
        os.utime(old, (now - 200, now - 200))
        os.utime(mid, (now - 100, now - 100))
        os.utime(new, (now,       now))

        assert _resolve_consensus_path(tmp_path) == new

    def test_lex_sort_would_pick_wrong_file_with_variable_width_pids(self, tmp_path):
        """A 4-digit PID sorts lexicographically before a 5-digit one,
        so name-sorted ``[0]`` (or even ``[-1]``) can pick the wrong
        file. Mtime-based selection ignores filename ordering.
        """
        # Create the 5-digit (older) file first, then the 4-digit (newer).
        old_5digit = tmp_path / "db.gene_structures_post_PASA_updates.55555.gff3"
        new_4digit = tmp_path / "db.gene_structures_post_PASA_updates.9999.gff3"
        old_5digit.write_text("# old\n")
        new_4digit.write_text("# new\n")

        now = time.time()
        os.utime(old_5digit, (now - 100, now - 100))
        os.utime(new_4digit, (now,       now))

        # Lex sort would put 55555 last (after 9999), but the newer file
        # is 9999. The mtime-based selection picks 9999.
        assert _resolve_consensus_path(tmp_path) == new_4digit
