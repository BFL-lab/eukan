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

# System packages needed for building fitild
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git autoconf \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 1. Create the conda environment with all bioconda tools + eukan
# ---------------------------------------------------------------------------
WORKDIR /opt/eukan

# Copy only dependency specs first (layer caching)
COPY environment.yml /opt/eukan/
COPY pyproject.toml poetry.lock README.md /opt/eukan/
COPY eukan/ /opt/eukan/eukan/
COPY tools.toml /opt/eukan/
COPY configs/ /opt/eukan/configs/

RUN mamba env create -f environment.yml \
    && mamba clean -afy

# Activate the eukan environment for all subsequent RUN commands
SHELL ["conda", "run", "-n", "eukan", "/bin/bash", "-c"]

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
              /opt/conda/envs/eukan/include \
              /opt/conda/envs/eukan/conda-meta \
              /opt/conda/envs/eukan/lib/python*/test \
              /opt/conda/envs/eukan/lib/python*/unittest \
              /opt/conda/envs/eukan/lib/python*/idlelib \
              /opt/conda/envs/eukan/lib/python*/tkinter \
              /opt/conda/envs/eukan/lib/python*/ensurepip \
              /opt/conda/envs/eukan/lib/python*/lib2to3 \
    && find /opt/conda/envs/eukan/lib -name "libpython*.a" -delete 2>/dev/null || true \
    && find /opt/conda/envs/eukan -name "*.egg-info" -type d -exec rm -rf {} + 2>/dev/null || true

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

# Copy the fully-built conda environment from the builder
COPY --link --from=builder /opt/conda/envs/eukan /opt/conda/envs/eukan

# Copy eukan source (tools.toml, configs)
COPY --link --from=builder /opt/eukan /opt/eukan

# Copy GeneMark key — place in eukan home and /root, world-readable
COPY --from=builder /root/.gm_key /root/.gm_key
COPY --from=builder /root/.gm_key /home/eukan/.gm_key
RUN chmod 644 /root/.gm_key /home/eukan/.gm_key 2>/dev/null || true

# AUGUSTUS config + spaln must be writable for species training
RUN chown -R ${USER_ID}:${GROUP_ID} /opt/conda/envs/eukan/config 2>/dev/null || true \
    && chown -R ${USER_ID}:${GROUP_ID} /opt/conda/envs/eukan/share/spaln 2>/dev/null || true

# ---------------------------------------------------------------------------
# Environment variables (mirrors conda-activate.sh)
# ---------------------------------------------------------------------------
ENV CONDA_DEFAULT_ENV=eukan \
    PATH="/opt/conda/envs/eukan/opt/genemark:/opt/conda/envs/eukan/opt/fitild/src:/opt/conda/envs/eukan/bin/EvmUtils:/opt/conda/envs/eukan/bin/EvmUtils/misc:/opt/conda/envs/eukan/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/conda/envs/eukan/lib" \
    ZOE="/opt/conda/envs/eukan/share/snap" \
    AUGUSTUS_CONFIG_PATH="/opt/conda/envs/eukan/config" \
    QUARRY_PATH="/opt/conda/envs/eukan/opt/codingquarry-2.0/QuarryFiles" \
    EVM_HOME="/opt/conda/envs/eukan/opt/evidencemodeler-2.1.0" \
    TRINITY_HOME="/opt/conda/envs/eukan/opt/trinity-2.15.2" \
    PASAHOME="/opt/conda/envs/eukan/opt/pasa-2.5.3" \
    ALN_TAB="/opt/conda/envs/eukan/share/spaln/table" \
    ALN_DBS="/opt/conda/envs/eukan/share/spaln/seqdb" \
    PERL5LIB="/opt/conda/envs/eukan/share/spaln/alndbs" \
    EUKAN_ROOT=/opt/eukan \
    HOME=/home/eukan \
    USER=eukan

WORKDIR /data
ENTRYPOINT ["eukan"]
