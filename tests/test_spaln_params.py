"""Tests for eukan.annotation.spaln_params (the experimental --spsp builder)."""

import os
from types import SimpleNamespace

from eukan.annotation import spaln_params

# Stand-in for perl: logs each call's argv as one tab-separated line, fails
# like perl on a script path that does not exist, and writes the files that
# make_eij.pl and make_ssp.pl leave behind.
_FAKE_PERL = """#!/bin/sh
(IFS="$(printf '\\t')"; echo "$*") >> "$ARGV_LOG"
[ -f "$1" ] || { echo "Can't open perl script \\"$1\\"" >&2; exit 2; }
case "$1" in
  *make_eij.pl) printf '#header\\nc1 100 200\\n' > sp_c.eij ;;
  *make_ssp.pl) touch Splice3 Splice5 IntronPotTab ;;
esac
"""


def test_make_ssp_runs_as_an_argv_list(tmp_path, monkeypatch):
    """make_ssp.pl gets its script path as one argument even with a space in it.

    It used to go through a shell string, which split such a path in two, so
    perl was handed the wrong script and Splice3/5 were never written.
    """
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name, body in (("perl", _FAKE_PERL), ("spaln", "#!/bin/sh\nexit 0\n")):
        (bindir / name).write_text(body)
        (bindir / name).chmod(0o755)
    tables = tmp_path / "spaln tables"
    alndbs = tables / "alndbs"
    alndbs.mkdir(parents=True)
    for script in ("make_eij.pl", "make_ssp.pl"):
        (alndbs / script).write_text("")
    argv_log = tmp_path / "argv.log"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("ALN_TAB", str(tables / "table"))
    monkeypatch.setenv("ARGV_LOG", str(argv_log))
    monkeypatch.setattr(spaln_params.time, "sleep", lambda _s: None)

    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\nACGT\n")
    transcripts = tmp_path / "transcripts.fa"
    transcripts.write_text(">t1\nACGT\n")
    config = SimpleNamespace(
        shortname="SP", genome=genome, transcripts_fasta=transcripts, num_cpu=1,
    )

    assert spaln_params.build_ssp(config, tmp_path) == "sp"

    calls = [line.split("\t") for line in argv_log.read_text().splitlines()]
    assert calls[-1] == [str(alndbs / "make_ssp.pl"), "-dsp_g", "-S", "sp_c.eij"]
    assert (tables / "table" / "sp" / "Splice3").exists()
