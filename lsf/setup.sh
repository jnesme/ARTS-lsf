#!/usr/bin/env bash
# ============================================================================
# ARTS environment setup
# Source this file to configure your shell for running artspipeline1.py.
#
# Usage:
#   source /work3/josne/github/arts/lsf/setup.sh
#
# Then run the pipeline directly, e.g.:
#   python "${ARTS_DIR}/artspipeline1.py" <input.gbk> "${ARTS_GAMMA_REFDIR}" \
#       -rd <resultdir> -cpu 4 -opt kres,duf,phyl -ast "${ASTRALJAR}"
# ============================================================================

# Resolve the ARTS repo root (this file lives in <repo>/lsf/)
ARTS_LSF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARTS_DIR="$(cd "${ARTS_LSF_DIR}/.." && pwd)"
export ARTS_DIR

# --- Conda ---
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate /work3/josne/miniconda3/envs/arts

# --- ARTS reference directories (phylum-specific, prebuilt) ---
export ARTS_GAMMA_REFDIR="${ARTS_DIR}/reference/gammaproteobacteria/"
export ARTS_ALPHA_REFDIR="${ARTS_DIR}/reference/alphaproteobacteria/"

# --- Astral jar (used by ARTS's -ast/--astral flag for species-tree building) ---
export ASTRALJAR="${ARTS_DIR}/astral/astral.5.7.7.jar"

# --- Known-resistance / DUF HMM databases (ARTS's own defaults, set explicitly
#     here for clarity/reproducibility in job logs) ---
export ARTS_KNOWNHMMS="${ARTS_DIR}/knownresistance.hmm"
export ARTS_DUFHMMS="${ARTS_DIR}/dufmodels.hmm"

echo "============================================"
echo " ARTS pipeline environment loaded"
echo "============================================"
echo " ARTS_DIR          : ${ARTS_DIR}"
echo " Gamma refdir      : ${ARTS_GAMMA_REFDIR}"
echo " Alpha refdir      : ${ARTS_ALPHA_REFDIR}"
echo " ASTRALJAR         : ${ASTRALJAR}"
echo " Known-resist HMMs : ${ARTS_KNOWNHMMS}"
echo " DUF HMMs          : ${ARTS_DUFHMMS}"
echo " Python            : $(python --version 2>&1)"
echo "============================================"
echo ""
