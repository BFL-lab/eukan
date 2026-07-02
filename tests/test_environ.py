"""Tests for eukan.infra.environ — subprocess environment construction."""

from __future__ import annotations

import os

from eukan.infra.environ import _prepend


class TestPrepend:
    def test_adds_component_when_absent(self):
        env = {"PATH": "/usr/bin:/bin"}
        _prepend(env, "PATH", "/opt/tool/bin")
        assert env["PATH"] == f"/opt/tool/bin{os.pathsep}/usr/bin:/bin"

    def test_skips_exact_duplicate(self):
        env = {"PATH": "/usr/bin:/bin"}
        _prepend(env, "PATH", "/usr/bin")
        assert env["PATH"] == "/usr/bin:/bin"

    def test_not_fooled_by_substring(self):
        """/foo must be added even though /foobar contains it as a substring."""
        env = {"PATH": f"/foobar{os.pathsep}/bin"}
        _prepend(env, "PATH", "/foo")
        components = env["PATH"].split(os.pathsep)
        assert components[0] == "/foo"
        assert "/foo" in components

    def test_prepend_into_empty_var(self):
        env: dict[str, str] = {}
        _prepend(env, "PATH", "/opt/bin")
        assert env["PATH"] == "/opt/bin"

    def test_prepend_into_empty_string_var(self):
        env = {"PATH": ""}
        _prepend(env, "PATH", "/opt/bin")
        assert env["PATH"] == "/opt/bin"
