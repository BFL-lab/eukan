"""Tests for the EVM commands.list line parser."""

from pathlib import Path

from eukan.annotation.evm import _parse_evm_command


class TestParseEvmCommand:
    def test_simple_command(self):
        result = _parse_evm_command("evidence_modeler.pl --genome g.fa")
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert argv == ["evidence_modeler.pl", "--genome", "g.fa"]
        assert cwd is None
        assert stdout is None
        assert stderr is None

    def test_stdout_redirect(self):
        result = _parse_evm_command("evidence_modeler.pl > consensus.out")
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert argv == ["evidence_modeler.pl"]
        assert stdout == "consensus.out"
        assert stderr is None

    def test_stderr_redirect(self):
        result = _parse_evm_command("evidence_modeler.pl 2> errs.log")
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert stderr == "errs.log"
        assert stdout is None

    def test_both_redirects(self):
        result = _parse_evm_command(
            "evidence_modeler.pl --weights w.txt > out.gff 2> err.log"
        )
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert argv == ["evidence_modeler.pl", "--weights", "w.txt"]
        assert stdout == "out.gff"
        assert stderr == "err.log"

    def test_cd_prefix(self):
        result = _parse_evm_command("cd /tmp/part1 && evidence_modeler.pl --x y")
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert cwd == Path("/tmp/part1")
        assert argv == ["evidence_modeler.pl", "--x", "y"]

    def test_cd_with_redirects(self):
        result = _parse_evm_command(
            "cd part && evidence_modeler.pl > out 2> err"
        )
        assert result is not None
        argv, cwd, stdout, stderr = result
        assert cwd == Path("part")
        assert argv == ["evidence_modeler.pl"]
        assert stdout == "out"
        assert stderr == "err"

    def test_pipe_returns_none(self):
        # Pipes must fall back to shell; we don't try to handle them.
        assert _parse_evm_command("a | b") is None

    def test_double_redirect_returns_none(self):
        assert _parse_evm_command("a > x > y") is None

    def test_empty_returns_none(self):
        assert _parse_evm_command("") is None

    def test_unbalanced_quotes_returns_none(self):
        # shlex.split raises on unclosed quotes
        assert _parse_evm_command('evidence_modeler.pl --foo "unterminated') is None

    def test_quoted_path_with_spaces(self):
        result = _parse_evm_command(
            'evidence_modeler.pl --weights "/tmp/with space/w.txt"'
        )
        assert result is not None
        argv, _, _, _ = result
        assert argv == ["evidence_modeler.pl", "--weights", "/tmp/with space/w.txt"]
