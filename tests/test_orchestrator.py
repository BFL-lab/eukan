"""Tests for eukan.annotation.orchestrator — CLI flag → manifest key translation."""

from __future__ import annotations

from eukan.annotation.orchestrator import force_steps_from_run_flags


class TestForceStepsFromRunFlags:
    """``--run-*`` CLI flags translate to the right manifest keys.

    The CLI surface is the dict of booleans accepted as kwargs here;
    the orchestrator surface is the list of ``annotation/<step>`` keys
    that get popped from the manifest before re-execution.
    """

    def test_no_flags_returns_empty(self):
        assert force_steps_from_run_flags() == []

    def test_run_genemark_groups_orf_finder(self):
        """--run-genemark forces both genemark and orf_finder (shared flag)."""
        result = force_steps_from_run_flags(run_genemark=True)
        assert set(result) == {"annotation/genemark", "annotation/orf_finder"}

    def test_run_snap_groups_codingquarry(self):
        """--run-snap forces both snap and codingquarry (shared flag)."""
        result = force_steps_from_run_flags(run_snap=True)
        assert set(result) == {"annotation/snap", "annotation/codingquarry"}

    def test_run_augustus_alone(self):
        assert force_steps_from_run_flags(run_augustus=True) == ["annotation/augustus"]

    def test_run_consensus_alone(self):
        assert force_steps_from_run_flags(run_consensus=True) == [
            "annotation/evm_consensus_models"
        ]

    def test_run_prot_align_default_picks_non_ssp(self):
        """Without spaln_ssp, --run-prot-align forces prot_align only."""
        result = force_steps_from_run_flags(run_prot_align=True)
        assert result == ["annotation/prot_align"]

    def test_run_prot_align_with_spsp_picks_ssp(self):
        """With spaln_ssp=True, --run-prot-align forces prot_align_ssp only."""
        result = force_steps_from_run_flags(run_prot_align=True, spaln_ssp=True)
        assert result == ["annotation/prot_align_ssp"]

    def test_spsp_alone_does_not_force(self):
        """spaln_ssp gates which prot-align step gets forced; alone it's a no-op."""
        assert force_steps_from_run_flags(spaln_ssp=True) == []

    def test_all_flags_together(self):
        """Every --run-* flag set forces every step exactly once."""
        result = force_steps_from_run_flags(
            spaln_ssp=False,
            run_genemark=True,
            run_prot_align=True,
            run_augustus=True,
            run_snap=True,
            run_consensus=True,
        )
        assert set(result) == {
            "annotation/genemark",
            "annotation/orf_finder",
            "annotation/prot_align",
            "annotation/augustus",
            "annotation/snap",
            "annotation/codingquarry",
            "annotation/evm_consensus_models",
        }
        assert "annotation/prot_align_ssp" not in result

    def test_all_flags_with_spsp(self):
        """Same as above but spaln_ssp swaps prot_align → prot_align_ssp."""
        result = force_steps_from_run_flags(
            spaln_ssp=True,
            run_genemark=True,
            run_prot_align=True,
            run_augustus=True,
            run_snap=True,
            run_consensus=True,
        )
        assert "annotation/prot_align_ssp" in result
        assert "annotation/prot_align" not in result

    def test_returned_keys_are_prefixed(self):
        """Every returned key carries the ``annotation/`` prefix."""
        result = force_steps_from_run_flags(
            run_genemark=True, run_augustus=True, run_snap=True,
        )
        assert all(k.startswith("annotation/") for k in result)
