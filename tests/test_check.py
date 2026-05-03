"""Tests for eukan.check — external tool and database pre-flight checks."""


from eukan.check import (
    PythonCheckResult,
    check_tool,
    format_results,
    generate_environment_yml,
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
    def test_all_pass(self):
        """Python dep checks should pass in a working environment."""
        from eukan.check import check_python_deps
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


class TestGenerateEnv:
    def test_generates_valid_yaml(self):
        """Output should be valid YAML with expected structure."""
        content = generate_environment_yml()
        assert "name: eukan" in content
        assert "bioconda" in content
        assert "augustus" in content
        assert "tools.toml" in content

    def test_deduplicates_packages(self):
        """spaln should appear only once despite spaln/makdbs sharing a package."""
        content = generate_environment_yml()
        assert content.count("- spaln") == 1
