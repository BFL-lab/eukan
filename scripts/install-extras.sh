#!/usr/bin/env bash
# Install tools that aren't available via conda: fitild and GeneMark.
#
# Run this after `conda activate eukan`. The script installs into
# $CONDA_PREFIX/opt/ so everything stays inside the conda environment.
#
# fitild is always built from source (GitHub).
#
# GeneMark requires a license — if gmes_linux_64_4.tar.gz and
# gm_key_64.gz are present in the project root, GeneMark is installed
# automatically. Otherwise it is skipped with a message.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: No conda environment is active. Run 'conda activate eukan' first." >&2
    exit 1
fi

OPT="${CONDA_PREFIX}/opt"
mkdir -p "$OPT"

# ---------------------------------------------------------------------------
# fitild — build from source
# ---------------------------------------------------------------------------
install_fitild() {
    echo "==> Installing fitild ..."

    if command -v fitild &>/dev/null; then
        echo "    fitild is already on PATH — skipping."
        return 0
    fi

    local dest="${OPT}/fitild"

    if [[ -d "$dest" ]]; then
        echo "    ${dest} already exists — rebuilding."
        rm -rf "$dest"
    fi

    env -u LD_LIBRARY_PATH git clone --depth 1 https://github.com/ogotoh/fitild "$dest"
    cd "$dest/src"

    # Fix C++ template ambiguity (max/fabs need explicit casts)
    sed -i 's/max(1\., fabs(b))/max((FTYPE)1., (FTYPE)fabs(b))/g' cmn.h

    CFLAGS="-O3 -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib" ./configure
    make -j"$(nproc)"

    # Merge ILD models into spaln's table directory so spaln can find them
    if [[ -d "${CONDA_PREFIX}/share/spaln/table" ]]; then
        cat "$dest"/table/IldModel*.txt > "${CONDA_PREFIX}/share/spaln/table/IldModel.txt" 2>/dev/null || true
        echo "    Merged IldModel tables into spaln table directory."
    fi

    echo "    fitild installed successfully."
}

# ---------------------------------------------------------------------------
# GeneMark — extract and configure (license required)
# ---------------------------------------------------------------------------
install_genemark() {
    if command -v gmes_petap.pl &>/dev/null; then
        echo "==> GeneMark is already on PATH — skipping."
        return 0
    fi

    local tar="${PROJECT_ROOT}/gmes_linux_64_4.tar.gz"
    local key="${PROJECT_ROOT}/gm_key_64.gz"

    if [[ ! -f "$tar" ]]; then
        echo "==> GeneMark: gmes_linux_64_4.tar.gz not found in project root — skipping."
        echo "    To install GeneMark, download the archive and license key from:"
        echo "      https://topaz.gatech.edu/GeneMark/license_download.cgi"
        echo "    Place gmes_linux_64_4.tar.gz and gm_key_64.gz in ${PROJECT_ROOT}/"
        echo "    Then re-run: ./scripts/install-extras.sh"
        return 0
    fi

    echo "==> Installing GeneMark ..."

    # Extract
    tar zxf "$tar" -C "$OPT"
    local gm_dir
    gm_dir=$(ls -d "${OPT}"/gmes_linux_64* 2>/dev/null | head -1)

    if [[ -z "$gm_dir" ]]; then
        echo "Error: extraction failed — no gmes_linux_64* directory in ${OPT}" >&2
        return 1
    fi

    # Symlink to a stable path
    ln -sfn "$gm_dir" "${OPT}/genemark"

    # Install the license key
    if [[ -f "$key" ]]; then
        gunzip -c "$key" > ~/.gm_key
        echo "    License key installed to ~/.gm_key"
    elif [[ ! -f ~/.gm_key ]]; then
        echo "Warning: gm_key_64.gz not found and ~/.gm_key does not exist." >&2
        echo "GeneMark will be installed but may not work without the license key." >&2
    fi

    # Fix shebangs to use conda's Perl (which has YAML, Hash::Merge, etc.)
    sed -i 's|#!/usr/bin/perl|#!/usr/bin/env perl|g' "${OPT}"/genemark/*.pl 2>/dev/null || true

    # Install MCE::Mutex (needed by GeneMark v4 parallel mode)
    if ! perl -MMCE::Mutex -e1 2>/dev/null; then
        echo "    Installing Perl module MCE::Mutex ..."
        cpanm --notest MCE::Mutex
    fi

    echo "    GeneMark installed successfully."
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
install_fitild
install_genemark

echo ""
echo "Done. Verify with: eukan check"
