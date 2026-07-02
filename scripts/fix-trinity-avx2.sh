#!/usr/bin/env bash
#
# Recompile Trinity's bundled C++ binaries for CPUs without AVX2.
#
# bioconda's Trinity (2.15.2) ships its Inchworm, Chrysalis and bamsifter
# binaries compiled with `-march=x86-64-v3`, which requires AVX2 + FMA + BMI.
# On pre-Haswell CPUs (e.g. the Ivy Bridge Xeon E5-26xx v2 family) every one
# dies with "Illegal instruction" — the first crash is Inchworm's
# `fastaToKmerCoverageStats` during read normalization.
#
# The conda package ships the full buildable source (CMakeLists/Makefiles), so
# this rebuilds those binaries in place at a safe `-march` (default x86-64-v2:
# SSE4.2/POPCNT, no AVX2 — runs on any x86-64-v2+ host). ParaFly and
# seqtk-trinity are already built to baseline and are left alone.
#
# Idempotent — safe to re-run (e.g. after `conda install/update trinity`).
#
# Requirements on PATH: g++, gcc, cmake, make (binutils' objdump for the
# optional verify). Run with the Trinity conda env active, or set
# TRINITY_LIBEXEC / CONDA_PREFIX explicitly.
#
# Env overrides:
#   MARCH            target arch (default: x86-64-v2)
#   CC / CXX         compilers (default: gcc / g++)
#   TRINITY_LIBEXEC  dir holding Inchworm/ Chrysalis/ trinity-plugins/
#                    (default: auto-detected from `command -v Trinity`)
#   CONDA_PREFIX     env prefix for htslib headers/libs (default: from env)
#
set -euo pipefail

MARCH="${MARCH:-x86-64-v2}"
CC="${CC:-gcc}"
CXX="${CXX:-g++}"

# --- locate Trinity's bundled source tree ----------------------------------
if [ -z "${TRINITY_LIBEXEC:-}" ]; then
    trin="$(command -v Trinity || true)"
    if [ -z "$trin" ]; then
        echo "ERROR: Trinity not on PATH. Activate the conda env or set TRINITY_LIBEXEC." >&2
        exit 1
    fi
    TRINITY_LIBEXEC="$(dirname "$(readlink -f "$trin")")"
fi

# bioconda has used two layouts: binaries directly in bin/, or under
# opt/trinity-*/. Fall back to a search if Inchworm/ isn't where we expect.
if [ ! -d "$TRINITY_LIBEXEC/Inchworm" ]; then
    search_root="${CONDA_PREFIX:-/opt/conda}"
    alt="$(find "$search_root" -maxdepth 7 -type f -name CMakeLists.txt \
             -path '*Inchworm*' 2>/dev/null | head -1 || true)"
    [ -n "$alt" ] && TRINITY_LIBEXEC="$(dirname "$(dirname "$alt")")"
fi
if [ ! -d "$TRINITY_LIBEXEC/Inchworm" ]; then
    echo "ERROR: cannot find Trinity's Inchworm/ source under '$TRINITY_LIBEXEC'." >&2
    exit 1
fi

PREFIX="${CONDA_PREFIX:-$(cd "$TRINITY_LIBEXEC/.." && pwd)}"

echo ">> Trinity bundled source : $TRINITY_LIBEXEC"
echo ">> conda prefix           : $PREFIX"
echo ">> target -march          : $MARCH"
echo ">> compiler               : $($CXX --version 2>/dev/null | head -1)"

FLAGS="-O3 -march=$MARCH -mtune=generic -fPIC"
export CC CXX
export CFLAGS="$FLAGS"
export CXXFLAGS="$FLAGS"
# cmake >= 4 rejects the old `cmake_minimum_required(VERSION 3.5)` floor
# unless this is set; harmless on cmake 3.x.
export CMAKE_POLICY_VERSION_MINIMUM="${CMAKE_POLICY_VERSION_MINIMUM:-3.5}"

# --- Inchworm + Chrysalis (cmake suites) -----------------------------------
# Each Makefile runs `cmake ... && make install`, picking up $CXXFLAGS. A clean
# build dir is required because the shipped CMakeCache points at the conda
# build host's (absent) compiler.
for suite in Inchworm Chrysalis; do
    echo ">> rebuilding $suite ..."
    ( cd "$TRINITY_LIBEXEC/$suite" && rm -rf build bin && make )
done

# --- bamsifter (single .cpp linking conda htslib) --------------------------
# Only the wrapped `_sift_bam_max_cov` binary carries v3 code; the htslib it
# loads (conda's libhts) uses runtime CPU dispatch and is already safe.
bs="$TRINITY_LIBEXEC/trinity-plugins/bamsifter"
if [ -f "$bs/sift_bam_max_cov.cpp" ]; then
    echo ">> rebuilding bamsifter (_sift_bam_max_cov) ..."
    ( cd "$bs" && "$CXX" -std=c++14 $CXXFLAGS -Wall \
        -o _sift_bam_max_cov sift_bam_max_cov.cpp \
        -I"$PREFIX/include" -L"$PREFIX/lib" -lhts -Wl,-rpath,"$PREFIX/lib" )
fi

# --- verify (best-effort) --------------------------------------------------
# Flag genuine AVX2/FMA/BMI mnemonics; SSE4.1 pextr[bwdq] are valid baseline.
if command -v objdump >/dev/null 2>&1; then
    echo ">> verifying rebuilt binaries are AVX2-free ..."
    v3re='^(vfmadd|vfmsub|vfnmadd|vfnmsub|vpbroadcast|vpgatherd|vperm2i128|vinserti128|vextracti128|vpblendd|vpsllv|vpsrlv|vpsrav|vpmaskmov|vpermd|vpermq|bzhi|mulx|pdep|sarx|shlx|shrx|rorx|tzcnt|lzcnt|blsr|blsi|blsmsk|bextr)'
    bad=0
    for f in "$TRINITY_LIBEXEC"/Inchworm/bin/* "$TRINITY_LIBEXEC"/Chrysalis/bin/* "$bs/_sift_bam_max_cov"; do
        [ -f "$f" ] || continue
        # `|| true`: a clean binary makes `grep -v` match nothing and exit 1,
        # which would abort the script under `set -e -o pipefail`.
        hits="$(objdump -d "$f" 2>/dev/null \
                  | awk -v p="$v3re" '{for(i=1;i<=NF;i++) if($i ~ p) print $i}' \
                  | { grep -vE '^pextr[bwdq]$' || true; } | sort -u | tr '\n' ' ')"
        if [ -n "$hits" ]; then echo "   !! $(basename "$f"): $hits"; bad=1; fi
    done
    if [ "$bad" -ne 0 ]; then
        echo ">> WARNING: residual v3 instructions detected (see above)." >&2
        exit 1
    fi
    echo ">> OK: all rebuilt Trinity binaries are AVX2-free (-march=$MARCH)."
else
    echo ">> (objdump not found — skipping verify)"
fi

echo ">> Done."
