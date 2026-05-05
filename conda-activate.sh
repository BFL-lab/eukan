#!/bin/bash
# Post-activation script for the eukan conda environment.
#
# NOTE: The `eukan` CLI sets these variables automatically at startup
# (see eukan/data/tools.toml and eukan/infra/environ.py).
# This script is only needed if you want to run the underlying tools
# (SNAP, spaln, EVM, etc.) directly outside of eukan.
#
# To install:
#   mkdir -p $CONDA_PREFIX/etc/conda/activate.d
#   cp conda-activate.sh $CONDA_PREFIX/etc/conda/activate.d/eukan.sh

# Ensure conda libraries are findable at runtime (needed by fitild).
# Note: this can break system tools (e.g. git) that link against system
# libcurl/libssh2. The eukan CLI sets this per-subprocess instead.
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# SNAP requires $ZOE to point to its data directory
export ZOE="${CONDA_PREFIX}/share/snap"

# EVM scripts are in a subdirectory — add to PATH
export PATH="${CONDA_PREFIX}/bin/EvmUtils:${CONDA_PREFIX}/bin/EvmUtils/misc:${PATH}"

# GeneMark (if installed into $CONDA_PREFIX/opt/genemark)
if [ -d "${CONDA_PREFIX}/opt/genemark" ]; then
    export PATH="${CONDA_PREFIX}/opt/genemark:${PATH}"
fi

# spaln table directory (needed by fitild for IldModel.txt)
if [ -d "${CONDA_PREFIX}/share/spaln/table" ]; then
    export ALN_TAB="${CONDA_PREFIX}/share/spaln/table"
fi
if [ -d "${CONDA_PREFIX}/share/spaln/seqdb" ]; then
    export ALN_DBS="${CONDA_PREFIX}/share/spaln/seqdb"
fi

# fitild (if installed into $CONDA_PREFIX/opt/fitild)
if [ -d "${CONDA_PREFIX}/opt/fitild/src" ]; then
    export PATH="${CONDA_PREFIX}/opt/fitild/src:${PATH}"
fi
