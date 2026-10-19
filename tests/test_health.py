"""Tests for eukan.infra.health."""

import importlib
import os
import subprocess
import sys
from pathlib import Path

import eukan
from eukan.cli.check import format_results
from eukan.infra import health
from eukan.infra.health import (
    PythonCheckResult,
    _check_module_imports,
    _crash_signal,
    check_tool,
    cpu_baseline,
    run_checks,
)
from eukan.infra.tools_registry import EnvVarSpec, Tool, load_tools


class TestCheckTool:
    def test_finds_python(self):
        """Python itself should always be found."""
        tool = Tool("Python", "python3", ("python3", "--version"), ("test",))
        result = check_tool(tool)
        assert result.found
        assert result.version_ok
        assert "Python" in result.version_output

    def test_missing_tool(self):
        """A nonexistent binary should fail cleanly."""
        tool = Tool("fake", "nonexistent_tool_xyz", ("nonexistent_tool_xyz",), ("test",))
        result = check_tool(tool)
        assert not result.found
        assert not result.version_ok

    def test_env_var_check(self):
        """Missing env var should be flagged."""
        tool = Tool(
            "test", "python3", ("python3", "--version"), ("test",),
            env_vars=(EnvVarSpec(var="NONEXISTENT_VAR_XYZ"),),
        )
        result = check_tool(tool)
        assert result.found
        assert not result.env_ok

    def test_non_utf8_version_output_is_ok(self, tmp_path):
        """A working tool whose version banner holds a non-UTF-8 byte (Debian's
        samtools prints a Latin-1 0xab) is not reported broken."""
        script = tmp_path / "fakesamtools.sh"
        script.write_text("printf 'samtools 1.21\\nbuilt in \\253BUILDPATH\\273\\n'\n")
        tool = Tool("fakesamtools", "sh", ("sh", str(script)), ("test",))
        result = check_tool(tool)
        assert result.found
        assert result.version_ok
        assert result.version_output == "samtools 1.21"

    @staticmethod
    def _fake_genemark(tmp_path, monkeypatch, check_install):
        """The registry's genemark entry, with a GeneMark dir on PATH whose
        gmes_petap.pl prints its usage banner and whose check_install.bash is
        *check_install*."""
        gm = tmp_path / "gmes_linux_64_4"
        gm.mkdir()
        petap = gm / "gmes_petap.pl"
        petap.write_text("#!/bin/sh\necho '# -------------------'\necho 'Usage: gmes_petap.pl'\n")
        petap.chmod(0o755)
        (gm / "check_install.bash").write_text(check_install)
        monkeypatch.setenv("PATH", f"{gm}{os.pathsep}{os.environ['PATH']}")
        (genemark,) = [t for t in load_tools() if t.name == "genemark"]
        return genemark

    def test_genemark_with_an_ended_license_fails(self, tmp_path, monkeypatch):
        """gmes_petap.pl runs without the key, so only GeneMark's own install
        check (gmhmme3) shows that the license period has ended."""
        genemark = self._fake_genemark(tmp_path, monkeypatch, (
            'echo "Checking GeneMark-ES installation"\n'
            'echo "GeneMark.hmm is set"\n'
            'echo "Your license period has ended" >&2\n'
            'echo "Error, unexpected condition was detected"\n'
            "exit 1\n"
        ))
        result = check_tool(genemark)
        assert result.problems == [
            "install check failed: Your license period has ended "
            "Error, unexpected condition was detected",
        ]

    def test_genemark_whose_install_check_passes_is_ok(self, tmp_path, monkeypatch):
        genemark = self._fake_genemark(
            tmp_path, monkeypatch, 'echo "All required components for GeneMark-ES were found"\n',
        )
        result = check_tool(genemark)
        assert result.ok
        assert result.version_output == "# -------------------"

    def test_direct_sigill_detected(self):
        """A binary killed by SIGILL should be flagged as broken."""
        tool = Tool(
            "selfkill", "python3",
            ("python3", "-c", "import os, signal; os.kill(os.getpid(), signal.SIGILL)"),
            ("test",),
        )
        result = check_tool(tool)
        assert result.found
        assert not result.version_ok
        assert result.crash_signal == "SIGILL"
        assert "SIGILL" in result.version_output

    def test_shell_encoded_sigill_detected(self):
        """A wrapper script that exits 128+4 (shell SIGILL) should be flagged."""
        tool = Tool(
            "wrappedkill", "bash",
            ("bash", "-c", "exit 132"),
            ("test",),
        )
        result = check_tool(tool)
        assert result.found
        assert not result.version_ok
        assert result.crash_signal == "SIGILL"

    def test_illegal_instruction_string_detected(self):
        """A wrapper that prints 'Illegal instruction' and exits 0 should still be flagged.

        Mirrors the bioconda STAR wrapper case: bash prints the signal name
        when the wrapped child dies, and the `for SIMD` loop falls through
        in a way that masks the encoded exit code. The string match is the
        only signal we have.
        """
        tool = Tool(
            "fakestar", "bash",
            ("bash", "-c", "echo 'Illegal instruction (core dumped)' >&2; exit 0"),
            ("test",),
        )
        result = check_tool(tool)
        assert result.found
        assert not result.version_ok
        assert result.crash_signal == "SIGILL"


class TestCrashSignal:
    """_crash_signal() handles both encodings: -N and 128+N."""

    def test_direct_signal_negative(self):
        assert _crash_signal(-4) == "SIGILL"
        assert _crash_signal(-11) == "SIGSEGV"

    def test_shell_encoded_signal(self):
        assert _crash_signal(132) == "SIGILL"
        assert _crash_signal(139) == "SIGSEGV"

    def test_normal_exit_codes_are_none(self):
        assert _crash_signal(0) is None
        assert _crash_signal(1) is None
        assert _crash_signal(127) is None  # command not found — not a signal
        assert _crash_signal(160) is None  # past the 128+31 cap


class TestCpuBaseline:
    """cpu_baseline() returns a recognized x86-64-vN level on Linux x86_64."""

    def test_returns_known_level_on_linux(self):
        result = cpu_baseline()
        if result is None:  # non-Linux or no /proc/cpuinfo
            return
        level, flags = result
        assert level in ("x86-64", "x86-64-v2", "x86-64-v3", "x86-64-v4")
        assert isinstance(flags, set)


class TestFormatSigill:
    """SIGILL crashes show CPU baseline + remediation hint."""

    def test_sigill_includes_baseline_hint(self):
        tool = Tool(
            "selfkill", "python3",
            ("python3", "-c", "import os, signal; os.kill(os.getpid(), signal.SIGILL)"),
            ("test",),
        )
        result = check_tool(tool)
        output = format_results([], [result])
        assert "SIGILL" in output
        # The CPU baseline footer should appear when any tool crashed.
        assert "CPU baseline:" in output


class TestToolRegistry:
    def test_loads_from_toml(self):
        """Should load tools from tools.toml."""
        tools = load_tools()
        assert len(tools) > 0
        names = [t.name for t in tools]
        assert "augustus" in names
        assert "samtools" in names

    def test_tool_fields(self):
        """Tools should have all required fields populated."""
        tools = load_tools()
        for tool in tools:
            assert tool.binary
            assert tool.version_cmd
            assert len(tool.required_by) > 0

    def test_minimap2_registered(self):
        """minimap2 (bioconda >= 2.29 for splice:sr) is registered for assemble."""
        mm2 = {t.name: t for t in load_tools()}["minimap2"]
        assert mm2.binary == "minimap2"
        assert mm2.conda_package == "minimap2"
        assert mm2.min_version == "2.29"
        assert "assemble" in mm2.required_by


class TestRunChecks:
    def test_filters_by_subcommand(self):
        """Should only check tools for the requested subcommand."""
        passed, failed, db_results, python_results = run_checks(["db-fetch"])
        all_tools = passed + failed
        for r in all_tools:
            assert "db-fetch" in r.tool.required_by

    def test_returns_db_results_for_func_annot(self, tmp_path):
        """Should include database checks when func-annot is in scope."""
        passed, failed, db_results, python_results = run_checks(["func-annot"], db_dir=tmp_path)
        assert len(db_results) > 0
        assert all(not ok for _, _, ok in db_results)

    def test_no_db_results_for_assemble(self):
        """Should not check databases for assemble-only."""
        passed, failed, db_results, python_results = run_checks(["assemble"])
        assert len(db_results) == 0


class TestPythonChecks:
    def test_module_check_covers_every_first_party_module(self):
        """Every module under eukan/ is imported, not a hand-kept list of 30."""
        root = Path(eukan.__file__).parent
        n_files = sum(1 for p in root.rglob("*.py") if "__pycache__" not in p.parts)
        result = _check_module_imports()
        assert result.ok
        assert result.detail == f"all {n_files} modules OK"

    def test_module_check_reports_a_module_that_fails_to_import(self, monkeypatch):
        """functional.kofam was not on the old list, so a broken one passed."""
        real = importlib.import_module

        def fake(name, *args, **kwargs):
            if name == "eukan.functional.kofam":
                raise ImportError("boom")
            return real(name, *args, **kwargs)

        monkeypatch.setattr(importlib, "import_module", fake)
        result = _check_module_imports()
        assert not result.ok
        assert "eukan.functional.kofam: boom" in result.detail

    def test_failed_probe_reports_why(self, monkeypatch):
        """A failed probe used to report an empty detail (a bare assert)."""
        monkeypatch.setattr(health, "load_tools", lambda: ())
        result = health._check_python_dep("tool registry", health._probe_tool_registry)
        assert not result.ok
        assert result.detail == "only 0 tools loaded from tools.toml"

    def test_probe_still_fails_under_python_O(self):
        """`python -O` strips asserts, which turned every probe into a pass."""
        code = (
            "from eukan.infra import health\n"
            "health.load_tools = lambda: ()\n"
            "print(health._check_python_dep('r', health._probe_tool_registry).ok)\n"
        )
        out = subprocess.run(
            [sys.executable, "-O", "-c", code], capture_output=True, text=True, check=True,
        ).stdout
        assert out.strip() == "False"

    def test_all_pass(self):
        """Python dep checks should pass in a working environment."""
        from eukan.infra.health import check_python_deps
        results = check_python_deps()
        assert len(results) > 0
        for r in results:
            assert r.ok, f"{r.name} failed: {r.detail}"

    def test_included_in_run_checks(self):
        """run_checks should return python_results."""
        passed, failed, db_results, python_results = run_checks(["assemble"])
        assert len(python_results) > 0


class TestFormatResults:
    def test_format_output(self):
        """Should produce readable output with counts."""
        tool = Tool("Python", "python3", ["python3", "--version"], ["test"])
        result = check_tool(tool)
        output = format_results([result], [])
        assert "1 tools OK" in output
        assert "Python" in output

    def test_format_with_db_results(self):
        """Should include database section when results provided."""
        db_ok = [("uniprot", "uniprot_sprot.faa OK (md5:abc...)", True)]
        db_fail = [("pfam", "Pfam-A.hmm not found", False)]
        output = format_results([], [], db_ok + db_fail)
        assert "1 databases OK" in output
        assert "1 databases MISSING" in output
        assert "eukan db-fetch" in output

    def test_format_with_python_results(self):
        """Should include Python section when results provided."""
        py_ok = [PythonCheckResult("pyhmmer", True, "works")]
        py_fail = [PythonCheckResult("missing_lib", False, "not installed")]
        output = format_results([], [], python_results=py_ok + py_fail)
        assert "1 Python checks OK" in output
        assert "1 Python checks FAILED" in output

    def test_install_hint_shown(self):
        """Missing tool with install_hint should show the hint."""
        tool = Tool("GeneMark", "gmes_petap.pl", ["gmes_petap.pl"], ["annotate"],
                     install_hint="Requires a license from https://topaz.gatech.edu/GeneMark/license_download.cgi")
        result = check_tool(tool)
        if not result.found:
            output = format_results([], [result])
            assert "license" in output.lower()
            assert "hint:" in output
