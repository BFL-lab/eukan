"""Tests for eukan.{annotation,assembly}.orchestrator — CLI flag → step translation."""

from __future__ import annotations

from typing import ClassVar

from eukan.annotation.orchestrator import force_steps_from_run_flags
from eukan.assembly.orchestrator import steps_and_force_from_run_flags


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


class TestStepsAndForceFromRunFlags:
    """``--run-X`` on assemble narrows the step list AND forces re-run.

    Regression: the old CLI mapped ``--run-X`` to ``steps=[X]`` only,
    leaving ``force=False``, so a completed manifest entry silently
    no-opped despite the help text claiming "Force re-run".
    """

    _ALL: ClassVar[list[str]] = ["star", "trinity", "pasa"]

    def test_no_flags_runs_all_no_force(self):
        steps, force = steps_and_force_from_run_flags()
        assert steps == self._ALL
        assert force is False

    def test_force_alone_runs_all_with_force(self):
        steps, force = steps_and_force_from_run_flags(force=True)
        assert steps == self._ALL
        assert force is True

    def test_run_star_alone_forces_star_only(self):
        """The bug we just fixed: --run-star alone must force-rerun star."""
        steps, force = steps_and_force_from_run_flags(run_star=True)
        assert steps == ["star"]
        assert force is True

    def test_run_trinity_alone_forces_trinity_only(self):
        steps, force = steps_and_force_from_run_flags(run_trinity=True)
        assert steps == ["trinity"]
        assert force is True

    def test_run_pasa_alone_forces_pasa_only(self):
        steps, force = steps_and_force_from_run_flags(run_pasa=True)
        assert steps == ["pasa"]
        assert force is True

    def test_run_star_with_force_is_idempotent(self):
        """--run-star --force still scopes to star (force was already implied)."""
        steps, force = steps_and_force_from_run_flags(run_star=True, force=True)
        assert steps == ["star"]
        assert force is True

    def test_multiple_run_flags(self):
        """Two --run-X flags select both steps and force them."""
        steps, force = steps_and_force_from_run_flags(run_star=True, run_pasa=True)
        assert steps == ["star", "pasa"]
        assert force is True

    def test_step_order_is_pipeline_order(self):
        """Returned steps follow pipeline order (star → trinity → pasa) regardless of kwarg order."""
        steps, _ = steps_and_force_from_run_flags(run_pasa=True, run_trinity=True, run_star=True)
        assert steps == ["star", "trinity", "pasa"]
