"""Tests for scripts/generate-env.py, the environment.yml generator."""

import runpy
import tomllib
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
generate_environment_yml = runpy.run_path(
    str(_REPO / "scripts" / "generate-env.py"),
)["generate_environment_yml"]


def test_environment_yml_is_up_to_date():
    """The checked-in environment.yml is exactly what tools.toml generates."""
    assert (_REPO / "environment.yml").read_text() == generate_environment_yml()


def test_generates_valid_yaml():
    """Output should be valid YAML with expected structure."""
    content = generate_environment_yml()
    assert "name: eukan" in content
    assert "bioconda" in content
    assert "augustus" in content
    assert "tools.toml" in content


def test_no_hmmer():
    """hmmer is not a conda dependency: pyhmmer replaced it."""
    assert "hmmer" not in generate_environment_yml()


def test_lists_spaln_once():
    content = generate_environment_yml()
    assert content.count("- spaln") == 1


def test_includes_minimap2():
    """minimap2 is emitted from tools.toml with its bioconda version pin."""
    content = generate_environment_yml()
    assert "minimap2>=2.29" in content


def test_python_floor_matches_pyproject():
    """The conda python floor tracks pyproject's; a lower one fails at `pip: .`."""
    pyproject = tomllib.loads((_REPO / "pyproject.toml").read_text())
    floor = pyproject["tool"]["poetry"]["dependencies"]["python"].split(",")[0]
    assert f"  - python{floor}," in generate_environment_yml()
