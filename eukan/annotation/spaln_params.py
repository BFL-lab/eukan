"""Species-specific spaln parameter generation (experimental ``--spsp``).

Builds per-species splice and intron-length parameter files for spaln
from transcript evidence, an alternative to the default fitild-based
intron length distribution. Selected by :class:`PipelineConfig.spaln_ssp`
and invoked from :func:`eukan.annotation.alignment._run_spaln`.
"""

from __future__ import annotations

import contextlib
import os
import shutil
import time
from pathlib import Path

from eukan.exceptions import ExternalToolError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_shell
from eukan.infra.utils import symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def build_ssp(config: PipelineConfig, sdir: Path) -> str:
    """Build species-specific spaln parameters from transcript data.

    Creates a species parameter directory under $ALN_TAB and populates it
    using make_eij.pl + make_ssp.pl.  Returns the species name for spaln -T.

    The spaln makessp scripts expect genome databases named ``<species>_g``
    and transcript files named ``<species>_c.cf``, so we follow that
    convention using ``config.shortname`` as the species identifier.

    Caller is responsible for ensuring ``config.transcripts_fasta`` is set
    (this function is only invoked when ``use_ssp`` is true in _run_spaln).
    """
    assert config.transcripts_fasta is not None
    ssp_name = config.shortname.lower()
    aln_tab = os.environ.get("ALN_TAB", "")
    if not aln_tab:
        raise RuntimeError(
            "ALN_TAB environment variable is not set. "
            "Ensure spaln is installed and eukan's environment is configured."
        )

    ssp_dir = Path(aln_tab) / ssp_name
    if ssp_dir.exists():
        shutil.rmtree(ssp_dir)
    ssp_dir.mkdir(parents=True)

    # spaln convention: genome = <ssp>_g, transcripts = <ssp>_c.cf
    genome_name = f"{ssp_name}_g"
    transcript_name = f"{ssp_name}_c"

    # Copy transcripts into the species parameter directory
    shutil.copy2(config.transcripts_fasta, ssp_dir / f"{transcript_name}.cf")

    # Symlink genome FASTA into ssp_dir and format for spaln
    symlink(config.genome, ssp_dir / f"{genome_name}.fna")
    run_cmd(
        ["spaln", "-W", "-KD", f"-t{config.num_cpu}", f"{genome_name}.fna"],
        cwd=ssp_dir,
    )

    log.info(
        "Building species-specific spaln parameters (%s) from %s",
        ssp_name, config.transcripts_fasta.name,
    )

    # Locate spaln perl scripts — prefer bundled versions in ALN_TAB/../alndbs
    # (which include Util.pm and match the compiled utility versions) over the
    # system-installed copies which may lag behind.
    # We invoke via "perl <script>" rather than executing the script directly
    # to ensure the conda perl is used, since the scripts' shebangs may
    # reference /usr/bin/perl (system perl) which can differ in version.
    alndbs = Path(aln_tab).parent / "alndbs"
    _make_eij = str(alndbs / "make_eij.pl") if (alndbs / "make_eij.pl").exists() else "make_eij.pl"
    _make_ssp = str(alndbs / "make_ssp.pl") if (alndbs / "make_ssp.pl").exists() else "make_ssp.pl"

    # Step 1: Map transcripts to genome, producing exon-intron junctions
    run_cmd(
        [
            "perl", _make_eij, f"-d{genome_name}",
            f"-t{config.num_cpu}", f"{transcript_name}.cf",
        ],
        cwd=ssp_dir,
    )

    eij_file = ssp_dir / f"{transcript_name}.eij"
    if not eij_file.exists():
        raise RuntimeError(
            f"make_eij.pl did not produce {eij_file}. "
            "Check that transcripts align to the genome."
        )

    # Count introns (header line starts with #)
    with open(eij_file) as fh:
        n_introns = sum(1 for line in fh if not line.startswith("#"))
    if n_introns < 5000:
        log.warning(
            "Only %d unique introns in %s (recommend >= 5000 for good parameters). "
            "Species-specific parameters may be unreliable.",
            n_introns, eij_file.name,
        )

    # Ensure file modification times are distinct — make_ssp.pl uses
    # Perl's -M operator to compare input/output ages and may skip
    # generation when files are created within the same second.
    time.sleep(1)

    # Step 2: Generate species-specific parameter files.
    # When -S is set, make_ssp.pl runs levels 9, 13, and 16.  Level 16
    # requires the fitild binary; if fitild is absent the script exits 1
    # but the critical output (Splice3/5, IntronPotTab) from levels 9 and
    # 13 has already been written.  We therefore tolerate a non-zero exit
    # code as long as the required files were produced.
    # Run via shell to ensure Perl's system() calls and shell redirects
    # inside make_ssp.pl work correctly.
    with contextlib.suppress(ExternalToolError):
        # check for output files below — non-zero exit is tolerated when
        # the level 9/13 outputs (Splice3/5, IntronPotTab) were produced
        run_shell(
            f"perl {_make_ssp} -d{genome_name} -S {transcript_name}.eij",
            cwd=ssp_dir,
        )

    # Verify key output files were created (may be .dat, .dgz, or extensionless)
    for base in ("Splice3", "Splice5", "IntronPotTab"):
        if not any((ssp_dir / f"{base}{ext}").exists() for ext in ("", ".dat", ".dgz")):
            raise RuntimeError(
                f"make_ssp.pl did not produce {base} in {ssp_dir}. "
                "Species-specific parameter generation may have failed."
            )

    log.info("Species-specific spaln parameters written to %s", ssp_dir)
    return ssp_name
