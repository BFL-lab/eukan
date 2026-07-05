# Dockerfile for eukan genome annotation pipeline.
#
# Multi-stage build:
#   Stage 1 (builder) — installs conda env, compiles fitild, installs GeneMark
#   Stage 2 (runtime) — copies only the conda env into a slim Debian base
#
# GeneMark is optional — if gmes_linux_64*.tar.gz and gm_key_64.gz are
# present in the build context, they'll be included. Otherwise the build
# succeeds and `eukan check` reports GeneMark as missing.

# ============================= BUILDER STAGE =================================
FROM condaforge/mambaforge:latest AS builder

# System packages needed for building fitild and spaln from source, plus
# curl/xz to fetch and unpack the pre-built combinr release tarball. jq must be
# on PATH before `mamba env create`: bioconductor data packages (pulled in
# transitively by trinity's R stack, e.g. bioconductor-genomeinfodbdata) run a
# post-link script that calls yq -> jq to parse their data-URL JSON, and the
# mambaforge base image does not ship jq — without it the env create aborts.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git autoconf curl ca-certificates xz-utils jq \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 1. Create the conda environment with all bioconda tools + eukan
# ---------------------------------------------------------------------------
WORKDIR /opt/eukan

# Copy only dependency specs first (layer caching)
COPY environment.yml /opt/eukan/
COPY pyproject.toml poetry.lock README.md /opt/eukan/
COPY eukan/ /opt/eukan/eukan/

RUN mamba env create -f environment.yml \
    && mamba clean -afy

# Activate the eukan environment for all subsequent RUN commands
SHELL ["conda", "run", "-n", "eukan", "/bin/bash", "-c"]

# Trinity's bundled C++ needs no special handling here: the conda env installs
# bioconda's current trinity build (compiled -march=x86-64-v3, AVX2/BMI2). The
# former x86-64-v2 pins (star/augustus/spaln/spades/stringtie/trinity) were
# dropped now that the deployment host is AVX2-capable, so this image requires
# an x86-64-v3 (Haswell+/Broadwell+) CPU at runtime.
# scripts/fix-trinity-avx2.sh remains only as a manual recovery tool for hosts
# without AVX2 (it rebuilds Trinity's binaries at a lower -march).
#
# TRINITY_HOME: bioconda's current trinity build (h077b44d_6) installs the whole
# Trinity tree directly under the env's bin/ (bin/Trinity, bin/util, bin/trinity-
# plugins, ...) — there is no longer an opt/trinity-<version> dir to symlink (the
# old x86-64-v2 _0 build had one). The package's own activate.d hook sets
# TRINITY_HOME=$PREFIX/bin, so the runtime ENV below points there; it is
# version-independent, so no per-version symlink is needed.

# ---------------------------------------------------------------------------
# 2. Build fitild from source (gsl + liblbfgs come from the conda env)
# ---------------------------------------------------------------------------
RUN git clone --depth 1 https://github.com/ogotoh/fitild /opt/conda/envs/eukan/opt/fitild \
    && cd /opt/conda/envs/eukan/opt/fitild/src \
    && sed -i 's/max(1\., fabs(b))/max((FTYPE)1., (FTYPE)fabs(b))/g' cmn.h \
    && ./configure \
    && make -j"$(nproc)" CFLAGS="-O3 -I/opt/conda/envs/eukan/include -L/opt/conda/envs/eukan/lib" \
    # Merge ILD models into spaln table directory
    && cat /opt/conda/envs/eukan/opt/fitild/table/IldModel*.txt \
       > /opt/conda/envs/eukan/share/spaln/table/IldModel.txt 2>/dev/null || true \
    # Clean build artifacts — keep only compiled binaries and table
    && cd /opt/conda/envs/eukan/opt/fitild \
    && rm -rf .git src/*.o

# ---------------------------------------------------------------------------
# 2b. Build spaln analysis utilities and install Perl scripts for makessp
#     The bioconda spaln package only ships the core binaries (spaln, sortgrcd,
#     makdbs, makmdm). Species-specific parameter generation (--spaln-ssp)
#     additionally needs: utn, dvn, exinpot, npssm, kmers, rdn, plus the
#     matching Perl scripts (make_eij.pl, make_ssp.pl, Util.pm, etc.).
# ---------------------------------------------------------------------------
RUN git clone --depth 1 https://github.com/ogotoh/spaln /tmp/spaln_src \
    && cd /tmp/spaln_src/src \
    && ./configure \
    && make -j"$(nproc)" all 2>/dev/null || true \
    && for bin in utn dvn exinpot npssm kmers rdn eij2ild; do \
        [ -f /tmp/spaln_src/src/$bin ] && cp /tmp/spaln_src/src/$bin /opt/conda/envs/eukan/bin/ || true; \
    done \
    && cp /tmp/spaln_src/perl/*.pl /tmp/spaln_src/perl/*.pm \
       /opt/conda/envs/eukan/share/spaln/alndbs/ \
    && chmod +x /opt/conda/envs/eukan/share/spaln/alndbs/*.pl \
    && rm -rf /tmp/spaln_src

# ---------------------------------------------------------------------------
# 2c. combinr — pre-built release binary (Rust tool, not on conda)
#
# combinr provides `combinr assemble` (transcript consolidation) and
# `combinr consensus` (weighted consensus gene models), both required by the
# assemble/annotate pipelines. Pull the pinned, statically-linked musl build
# from the GitHub release, verify its SHA-256, and install onto PATH (the conda
# env's bin/, already first on PATH in the runtime stage). Bump VERSION and
# SHA256 together to upgrade; they are validated against each other at build.
# ---------------------------------------------------------------------------
ARG COMBINR_VERSION=0.1.1
ARG COMBINR_SHA256=913c9698124cc3891c01d738cb775bcd5918f7af57ddd499108420d1aed7c222
RUN combinr_dir="combinr-x86_64-unknown-linux-musl" \
    && curl -fsSL -o "/tmp/${combinr_dir}.tar.xz" \
       "https://github.com/BFL-lab/combinr/releases/download/v${COMBINR_VERSION}/${combinr_dir}.tar.xz" \
    && echo "${COMBINR_SHA256}  /tmp/${combinr_dir}.tar.xz" | sha256sum -c - \
    && tar xJf "/tmp/${combinr_dir}.tar.xz" -C /tmp \
    && install -m 0755 "/tmp/${combinr_dir}/combinr" /opt/conda/envs/eukan/bin/combinr \
    && rm -rf "/tmp/${combinr_dir}.tar.xz" "/tmp/${combinr_dir}" \
    && combinr --version

# ---------------------------------------------------------------------------
# 3. GeneMark (optional — requires user-provided license files)
# ---------------------------------------------------------------------------
COPY gmes_linux_64*.tar.gz gm_key_64.gz* /tmp/_genemark/
RUN genemark_tar=$(ls /tmp/_genemark/gmes_linux_64*.tar.gz 2>/dev/null | head -1) \
    && if [ -n "$genemark_tar" ]; then \
        tar zxf "$genemark_tar" -C /opt/conda/envs/eukan/opt/ \
        && gm_dir=$(ls -d /opt/conda/envs/eukan/opt/gmes_linux_64* 2>/dev/null | head -1) \
        && ln -sfn "$gm_dir" /opt/conda/envs/eukan/opt/genemark \
        && gunzip -c /tmp/_genemark/gm_key_64.gz > /root/.gm_key 2>/dev/null || true \
        && sed -i 's|#!/usr/bin/perl|#!/usr/bin/env perl|g' /opt/conda/envs/eukan/opt/genemark/*.pl 2>/dev/null || true \
        && echo "GeneMark installed from $genemark_tar"; \
    else \
        echo "GeneMark files not found — skipping (install manually later)"; \
    fi && rm -rf /tmp/_genemark

# ---------------------------------------------------------------------------
# 4. Strip conda env of unnecessary files to shrink the runtime image
#
# Removed in addition to the obvious caches/test dirs:
#   - R interpreter + bioconductor (lib/R, bin/R, bin/Rscript) — never invoked
#   - pandoc (155 MB Haskell binary) — never invoked
#   - TEsorter HMM-database package — RepeatClassifier doesn't call it
#   - Cross-compile sysroot + gcc toolchain (x86_64-conda-linux-gnu,
#     libexec/gcc, lib/gcc) — build-time only; fitild/spaln are already built
#   - JDK source zip and jmods — JRE is what Trinity needs
#   - Berkeley DB docs, man pages, locale data, gir-1.0, fonts, cups/X11/alsa
#   - lib/ncbi-blast+ (391 MB BLAST+ dev SDK) — verified: no binary links to it
#   - Unused BLAST CLIs (tblastx, psiblast, deltablast, rps*blast, seedtop,
#     blast_formatter/_report, blastdbcheck) — not invoked by eukan or any
#     bundled tool (RepeatMasker, RepeatModeler, PASA, Trinity)
#   - matplotlib + fontTools + contourpy + PIL + RNA python packages — only
#     pulled in via lazy imports in scipy/biopython that we never trigger
#   - AUGUSTUS pre-shipped species (165 of 166) — eukan always trains from
#     scratch via new_species.pl using the `generic` template
#   - GeneMark GeneMark-E-tests + ProtHint — test data + EP+ mode (unused)
#   - bin/k8 — JS shell for minimap2's paftools.js only; eukan calls the
#     minimap2 binary directly (BAM via samtools/pysam), never paftools, so the
#     read/transcript mapping is unaffected by removing it
# ---------------------------------------------------------------------------
RUN find /opt/conda/envs/eukan -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda/envs/eukan -name "*.a" -delete 2>/dev/null || true \
    && find /opt/conda/envs/eukan -name "*.pyc" -delete 2>/dev/null || true \
    && find /opt/conda/envs/eukan -name "*.pyo" -delete 2>/dev/null || true \
    && find /opt/conda/envs/eukan -type d -name "tests" -path "*/site-packages/*" -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda/envs/eukan -type d -name "test" -path "*/site-packages/*" -exec rm -rf {} + 2>/dev/null || true \
    && rm -rf /opt/conda/envs/eukan/share/doc \
              /opt/conda/envs/eukan/share/man \
              /opt/conda/envs/eukan/share/info \
              /opt/conda/envs/eukan/share/gtk-doc \
              /opt/conda/envs/eukan/share/terminfo \
              /opt/conda/envs/eukan/share/locale \
              /opt/conda/envs/eukan/share/gir-1.0 \
              /opt/conda/envs/eukan/share/cups \
              /opt/conda/envs/eukan/share/X11 \
              /opt/conda/envs/eukan/share/alsa \
              /opt/conda/envs/eukan/share/hdf5_examples \
              /opt/conda/envs/eukan/include \
              /opt/conda/envs/eukan/conda-meta \
              /opt/conda/envs/eukan/man \
              /opt/conda/envs/eukan/docs \
              /opt/conda/envs/eukan/fonts \
              /opt/conda/envs/eukan/x86_64-conda-linux-gnu \
              /opt/conda/envs/eukan/libexec/gcc \
              /opt/conda/envs/eukan/lib/gcc \
              /opt/conda/envs/eukan/lib/ncbi-blast+ \
              /opt/conda/envs/eukan/lib/R \
              /opt/conda/envs/eukan/lib/jvm/jmods \
              /opt/conda/envs/eukan/lib/jvm/lib/src.zip \
              /opt/conda/envs/eukan/lib/python*/test \
              /opt/conda/envs/eukan/lib/python*/idlelib \
              /opt/conda/envs/eukan/lib/python*/tkinter \
              /opt/conda/envs/eukan/lib/python*/ensurepip \
              /opt/conda/envs/eukan/lib/python*/lib2to3 \
              /opt/conda/envs/eukan/lib/python*/site-packages/TEsorter \
              /opt/conda/envs/eukan/lib/python*/site-packages/matplotlib \
              /opt/conda/envs/eukan/lib/python*/site-packages/matplotlib-* \
              /opt/conda/envs/eukan/lib/python*/site-packages/mpl_toolkits \
              /opt/conda/envs/eukan/lib/python*/site-packages/fontTools \
              /opt/conda/envs/eukan/lib/python*/site-packages/contourpy \
              /opt/conda/envs/eukan/lib/python*/site-packages/contourpy-* \
              /opt/conda/envs/eukan/lib/python*/site-packages/PIL \
              /opt/conda/envs/eukan/lib/python*/site-packages/pillow-* \
              /opt/conda/envs/eukan/lib/python*/site-packages/RNA \
    && rm -f  /opt/conda/envs/eukan/bin/R \
              /opt/conda/envs/eukan/bin/Rscript \
              /opt/conda/envs/eukan/bin/pandoc \
              /opt/conda/envs/eukan/bin/TEsorter \
              /opt/conda/envs/eukan/bin/TEsorter-test \
              /opt/conda/envs/eukan/bin/concatenate_domains.py \
              /opt/conda/envs/eukan/bin/k8 \
              /opt/conda/envs/eukan/bin/tblastx \
              /opt/conda/envs/eukan/bin/psiblast \
              /opt/conda/envs/eukan/bin/deltablast \
              /opt/conda/envs/eukan/bin/rpsblast \
              /opt/conda/envs/eukan/bin/rpstblastn \
              /opt/conda/envs/eukan/bin/seedtop \
              /opt/conda/envs/eukan/bin/blast_formatter \
              /opt/conda/envs/eukan/bin/blast_report \
              /opt/conda/envs/eukan/bin/blastdbcheck \
    # AUGUSTUS species: keep only `generic` (used by new_species.pl as the
    # training template). Eukan always trains from scratch — pre-shipped
    # reference species (~140 MB) are never consulted.
    && find /opt/conda/envs/eukan/config/species -mindepth 1 -maxdepth 1 -type d \
              ! -name generic -exec rm -rf {} + 2>/dev/null || true \
    # GeneMark: drop test inputs (~11 MB) and ProtHint (~13 MB, EP+ mode unused).
    && rm -rf /opt/conda/envs/eukan/opt/gmes_linux_64_4/GeneMark-E-tests \
              /opt/conda/envs/eukan/opt/gmes_linux_64_4/ProtHint \
    && find /opt/conda/envs/eukan/lib -name "libpython*.a" -delete 2>/dev/null || true \
    && find /opt/conda/envs/eukan -name "*.egg-info" -type d -exec rm -rf {} + 2>/dev/null || true \
    # Pre-set writable permissions on dirs the runtime needs to mutate
    # (AUGUSTUS species training writes to config/, spaln --spsp writes to
    # share/spaln/). Doing this in the builder avoids a 191 MB chown layer
    # in the runtime stage since COPY --link preserves these as-is.
    && chmod -R a+rwX /opt/conda/envs/eukan/config /opt/conda/envs/eukan/share/spaln

# ============================= RUNTIME STAGE =================================
FROM debian:bookworm-slim AS runtime

ARG USER_ID=1000
ARG GROUP_ID=1000

# Runtime setup: minimal system deps + non-root user + GeneMark key (single layer)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/* \
    # Create non-root user
    && groupadd -g ${GROUP_ID} eukan 2>/dev/null || true \
    && useradd -u ${USER_ID} -g ${GROUP_ID} -m eukan 2>/dev/null || true

# Copy the fully-built conda environment from the builder.
# eukan itself (including tools.toml + configs/ under eukan/data/) lives in
# site-packages via `pip install .`, so no separate source copy is needed.
COPY --link --from=builder /opt/conda/envs/eukan /opt/conda/envs/eukan

# Copy GeneMark key — place in eukan home and /root, world-readable
COPY --from=builder /root/.gm_key /root/.gm_key
COPY --from=builder /root/.gm_key /home/eukan/.gm_key
RUN chmod 644 /root/.gm_key /home/eukan/.gm_key 2>/dev/null || true

# AUGUSTUS config + share/spaln are already chmod a+rwX in the builder so
# any UID can write during species training / spaln --spsp.

# ---------------------------------------------------------------------------
# Environment variables (mirrors conda-activate.sh)
# ---------------------------------------------------------------------------
ENV CONDA_DEFAULT_ENV=eukan \
    PATH="/opt/conda/envs/eukan/opt/genemark:/opt/conda/envs/eukan/opt/fitild/src:/opt/conda/envs/eukan/bin/EvmUtils:/opt/conda/envs/eukan/bin/EvmUtils/misc:/opt/conda/envs/eukan/opt/pasa-2.5.3/scripts:/opt/conda/envs/eukan/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/conda/envs/eukan/lib" \
    ZOE="/opt/conda/envs/eukan/share/snap" \
    AUGUSTUS_CONFIG_PATH="/opt/conda/envs/eukan/config" \
    QUARRY_PATH="/opt/conda/envs/eukan/opt/codingquarry-2.0/QuarryFiles" \
    EVM_HOME="/opt/conda/envs/eukan/opt/evidencemodeler-2.1.0" \
    TRINITY_HOME="/opt/conda/envs/eukan/bin" \
    PASAHOME="/opt/conda/envs/eukan/opt/pasa-2.5.3" \
    ALN_TAB="/opt/conda/envs/eukan/share/spaln/table" \
    ALN_DBS="/opt/conda/envs/eukan/share/spaln/seqdb" \
    PERL5LIB="/opt/conda/envs/eukan/share/spaln/alndbs" \
    HOME=/home/eukan \
    USER=eukan

WORKDIR /data
ENTRYPOINT ["eukan"]
