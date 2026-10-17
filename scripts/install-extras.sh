#!/usr/bin/env bash
# Install tools that aren't available via conda: fitild, combinr, and GeneMark.
#
# Run this after `conda activate eukan`. The script installs into
# $CONDA_PREFIX/ so everything stays inside the conda environment.
#
# With no arguments it installs all three; name components to install only
# those, e.g. `install-extras.sh combinr` (which is how the Docker image
# installs combinr).
#
# fitild is always built from source (GitHub).
#
# combinr is downloaded as a pinned pre-built release binary and dropped onto
# PATH.
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

# Pinned combinr release (bump VERSION and the matching SHA-256 below together;
# after a new release builds, run scripts/refresh-combinr-shas.py <version> to
# fill the per-platform SHA-256s automatically).
COMBINR_VERSION="0.3.0"

# Verify $2's SHA-256 against $1 using whichever tool is present (sha256sum on
# Linux, shasum on macOS). Returns non-zero on mismatch.
verify_sha256() {
    local expected="$1" file="$2" actual
    if command -v sha256sum &>/dev/null; then
        actual="$(sha256sum "$file" | awk '{print $1}')"
    elif command -v shasum &>/dev/null; then
        actual="$(shasum -a 256 "$file" | awk '{print $1}')"
    else
        echo "Warning: no sha256 tool found — skipping checksum verification." >&2
        return 0
    fi
    [[ "$actual" == "$expected" ]]
}

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
# combinr — download the pinned pre-built release binary
# ---------------------------------------------------------------------------
install_combinr() {
    echo "==> Installing combinr ${COMBINR_VERSION} ..."

    if command -v combinr &>/dev/null; then
        echo "    combinr is already on PATH — skipping."
        return 0
    fi

    local os arch target sha
    os="$(uname -s)"
    arch="$(uname -m)"

    case "${os}:${arch}" in
        Linux:x86_64)
            target="x86_64-unknown-linux-musl"
            sha="f937ec44aa28d2ab50f98995a5d879c42a303d347f112c5fbdde479f113dde14" ;;
        Linux:aarch64|Linux:arm64)
            target="aarch64-unknown-linux-musl"
            sha="4ca31949a5ec430b6295bee7e2f69166a8822fff03fe5c13f06acac650d0bcca" ;;
        Darwin:x86_64)
            target="x86_64-apple-darwin"
            sha="74167aa9495ccd71f87d0905abef6ac4857cacc5362fc62e45d350a4e318d54a" ;;
        Darwin:arm64|Darwin:aarch64)
            target="aarch64-apple-darwin"
            sha="2966a502e5dd38deeab79347576f2ea2d316ba4a893cad28f83ac9ad98ef2b28" ;;
        *)
            echo "    No pre-built combinr binary for ${os}/${arch}." >&2
            echo "    Download a release from https://github.com/BFL-lab/combinr/releases" >&2
            echo "    and put 'combinr' on PATH, or set EUKAN_ASSEMBLE_COMBINR_PATH." >&2
            return 0 ;;
    esac

    # The release asset is .tar.xz, so unpacking needs xz on PATH. Check up
    # front for a clear message (the Docker build apt-installs xz-utils; conda
    # usually pulls xz in transitively, but it isn't a declared dependency).
    if ! command -v xz &>/dev/null; then
        echo "Error: 'xz' is required to unpack the combinr release tarball but was not found." >&2
        echo "    Install it (e.g. 'conda install -c conda-forge xz') and re-run." >&2
        return 1
    fi

    local tarball="combinr-${target}.tar.xz"
    local url="https://github.com/BFL-lab/combinr/releases/download/v${COMBINR_VERSION}/${tarball}"
    local tmp
    tmp="$(mktemp -d)"

    # Clean up the temp dir on every failure path (under `set -e` a bare
    # failure would otherwise abort before the trailing rm -rf).
    env -u LD_LIBRARY_PATH curl -fsSL -o "${tmp}/${tarball}" "$url" \
        || { echo "Error: failed to download combinr from ${url}." >&2; rm -rf "$tmp"; return 1; }
    if ! verify_sha256 "$sha" "${tmp}/${tarball}"; then
        echo "Error: combinr checksum mismatch for ${tarball}." >&2
        rm -rf "$tmp"
        return 1
    fi
    tar xJf "${tmp}/${tarball}" -C "$tmp" \
        || { echo "Error: failed to unpack ${tarball}." >&2; rm -rf "$tmp"; return 1; }
    install -m 0755 "${tmp}/combinr-${target}/combinr" "${CONDA_PREFIX}/bin/combinr" \
        || { echo "Error: failed to install combinr binary." >&2; rm -rf "$tmp"; return 1; }
    rm -rf "$tmp"

    echo "    combinr installed to ${CONDA_PREFIX}/bin/combinr."
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
if [[ $# -eq 0 ]]; then
    set -- fitild combinr genemark
fi
for component in "$@"; do
    case "$component" in
        fitild|combinr|genemark) "install_${component}" ;;
        *)
            echo "Error: unknown component '${component}' (choose from fitild, combinr, genemark)." >&2
            exit 2 ;;
    esac
done

echo ""
echo "Done. Verify with: eukan check"
