#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -J arts_smoketest_vibrio
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8400MB
#BSUB -W 3:00
#BSUB -o /work3/josne/github/arts/lsf/logs/arts_smoketest_vibrio_%J.out
#BSUB -e /work3/josne/github/arts/lsf/logs/arts_smoketest_vibrio_%J.err

set -euo pipefail

#==========================================================================
# ARTS smoke test — Vibrio sp. S0204 (solo run)
# Input is antiSMASH 8.0.1 output already produced by the funcscan workflow —
# NOT re-running antiSMASH here (ARTS reads the antiSMASH-annotated GenBank
# directly).
#==========================================================================
INPUT="/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/funcscan_results/bgc/antismash/S0204/S0204.gbk"
RESULTDIR="/work3/josne/Projects/Vibrio_Galathea3/arts_results/smoketest/vibrio_only"
CPU=4
#==========================================================================

ARTS_DIR="/work3/josne/github/arts"

if [ ! -f "${INPUT}" ]; then
    echo "ERROR: Input GenBank file not found: ${INPUT}"
    exit 1
fi

mkdir -p "${RESULTDIR}"
mkdir -p "${ARTS_DIR}/lsf/logs"

# Load environment (conda, ASTRALJAR, refdir vars)
source "${ARTS_DIR}/lsf/setup.sh"

echo "=========================================="
echo "ARTS smoke test - Vibrio sp. S0204 (solo)"
echo "Host       : $(hostname)"
echo "Job started: $(date)"
echo "LSF job ID : ${LSB_JOBID:-N/A}"
echo "Input      : ${INPUT}"
echo "Refdir     : ${ARTS_GAMMA_REFDIR}"
echo "Resultdir  : ${RESULTDIR}"
echo "CPUs       : ${CPU}"
echo "=========================================="

cd "${ARTS_DIR}"

set +e
python artspipeline1.py \
    "${INPUT}" \
    "${ARTS_GAMMA_REFDIR}" \
    -rd "${RESULTDIR}" \
    -cpu "${CPU}" \
    -opt kres,duf,phyl \
    -khmms "${ARTS_KNOWNHMMS}" \
    -duf "${ARTS_DUFHMMS}" \
    -ast "${ASTRALJAR}"
RUN_STATUS=$?
set -e

echo "=========================================="
echo "ARTS run finished with exit code: ${RUN_STATUS}"
echo "Host       : $(hostname)"
echo "Job ended  : $(date)"
echo "=========================================="

# --- Post-run sanity check ---
echo ""
echo "--- Sanity check: output tables ---"
CORETABLE="${RESULTDIR}/tables/coretable.tsv"
BGCTABLE="${RESULTDIR}/tables/bgctable.tsv"

check_table() {
    local f="$1"
    if [ ! -s "${f}" ]; then
        echo "FAIL: ${f} missing or empty"
        return 1
    fi
    local nlines
    nlines=$(wc -l < "${f}")
    if [ "${nlines}" -le 1 ]; then
        echo "FAIL: ${f} has only a header (or is empty) - ${nlines} line(s)"
        return 1
    fi
    echo "OK: ${f} - ${nlines} lines"
    return 0
}

SANITY_OK=1
check_table "${CORETABLE}" || SANITY_OK=0
check_table "${BGCTABLE}" || SANITY_OK=0

if [ "${SANITY_OK}" -eq 1 ]; then
    echo "Sanity check PASSED"
else
    echo "Sanity check FAILED - see above. This is the exact failure mode seen"
    echo "previously when ARTS was fed a non-antiSMASH-annotated GenBank file"
    echo "(bgctable.tsv came back header-only). Verify the input .gbk still"
    echo "carries the '##antiSMASH-Data-START##' COMMENT marker and"
    echo "region_number qualifiers."
fi

exit "${RUN_STATUS}"
