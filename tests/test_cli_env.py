"""EUKAN_* environment variables through the CLI: an option left out falls
back to its variable, a given option wins, and with neither the config's
default -- the one --help shows -- applies."""

from __future__ import annotations

import os

import pytest
from click.testing import CliRunner

import eukan.annotation
import eukan.assembly
import eukan.functional
import eukan.repeats
import eukan.submission
from eukan.cli import cli

_DRIVERS = {
    "annotate": (eukan.annotation, "run_annotation_pipeline"),
    "assemble": (eukan.assembly, "run_assembly"),
    "func-annot": (eukan.functional, "run_functional_annotation"),
    "mask-repeats": (eukan.repeats, "run_repeats"),
    "prep-submission": (eukan.submission, "run_prep_submission"),
}
_ARGS = {
    "annotate": ["-g", "g.fa", "-p", "p.fa"],
    "assemble": ["-g", "g.fa", "-l", "l.fq", "-r", "r.fq"],
    "func-annot": ["-p", "p.fa"],
    "mask-repeats": ["-g", "g.fa"],
    "prep-submission": ["-t", "t.sbt", "--organism", "X y", "-g", "g.fa", "-i", "a.gff3"],
}


@pytest.fixture
def invoke(tmp_path, monkeypatch):
    for var in [v for v in os.environ if v.startswith("EUKAN_")]:
        monkeypatch.delenv(var)  # only what each test sets
    for name in ("g.fa", "p.fa", "l.fq", "r.fq", "t.sbt", "a.gff3"):
        (tmp_path / name).write_text(">x\nACGT\n")
    monkeypatch.chdir(tmp_path)
    seen = {}
    for command, (module, attr) in _DRIVERS.items():
        monkeypatch.setattr(module, attr, lambda c, _k=command, **kw: seen.setdefault(_k, c))

    def _invoke(command, *extra, env=None):
        result = CliRunner().invoke(cli, [command, *_ARGS[command], *extra], env=env)
        assert result.exit_code == 0, result.output
        return seen.pop(command)

    return _invoke


@pytest.mark.parametrize(("command", "var", "value", "field", "expected"), [
    ("annotate", "EUKAN_NUM_CPU", "3", "num_cpu", 3),
    ("annotate", "EUKAN_GENETIC_CODE", "6", "genetic_code", "6"),
    ("annotate", "EUKAN_WEIGHTS", "[5,5,5]", "weights", [5, 5, 5]),
    ("assemble", "EUKAN_ASSEMBLE_MAX_INTRON_LEN", "123", "max_intron_len", 123),
    ("assemble", "EUKAN_ASSEMBLE_DIAGNOSE_SOFTCLIPS", "false", "diagnose_softclips", False),
    ("func-annot", "EUKAN_FUNC_EVALUE", "1e-9", "evalue", 1e-9),
    ("mask-repeats", "EUKAN_REPEATS_ENGINE", "ncbi", "engine", "ncbi"),
    ("prep-submission", "EUKAN_SUBMIT_MODE", "b", "mode", "b"),
])
def test_env_var_fills_an_option_left_out(invoke, command, var, value, field, expected):
    assert getattr(invoke(command, env={var: value}), field) == expected


@pytest.mark.parametrize(("command", "flag", "env", "field", "expected"), [
    ("annotate", ["-n", "5"], {"EUKAN_NUM_CPU": "3"}, "num_cpu", 5),
    # given on the command line, even when equal to the default
    ("annotate", ["-c", "11"], {"EUKAN_GENETIC_CODE": "6"}, "genetic_code", "11"),
    ("assemble", ["-M", "777"], {"EUKAN_ASSEMBLE_MAX_INTRON_LEN": "123"}, "max_intron_len", 777),
])
def test_given_option_beats_env_var(invoke, command, flag, env, field, expected):
    assert getattr(invoke(command, *flag, env=env), field) == expected


@pytest.mark.parametrize("command", list(_DRIVERS))
def test_click_defaults_match_the_config_defaults(invoke, command):
    """--help shows the click default; with the option left out the config's
    default applies, so the two must agree."""
    config = invoke(command)
    explicit = []
    for p in cli.commands[command].params:
        if p.name in ("force", "print_command", "dry_run") or p.name.startswith("run_"):
            continue
        if p.default is False:  # a flag cannot be given as off; its field must be
            assert getattr(config, p.name) is False, p.name
        elif p.default is True:
            explicit.append(p.opts[0])
        elif isinstance(p.default, (int, float, str)) and p.default != "":
            explicit += [p.opts[0], str(p.default)]
        elif isinstance(p.default, tuple):
            for value in p.default:
                explicit += [p.opts[0], str(value)]

    def dump(config):
        return config.model_dump(exclude={"shortname"})

    assert dump(invoke(command, *explicit)) == dump(config)
