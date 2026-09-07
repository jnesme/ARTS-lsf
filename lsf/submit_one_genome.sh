#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8400MB
#BSUB -W 3:00

set -euo pipefail

#==========================================================================
# Generic single-genome ARTS submission script.
# Feeds an already-computed antiSMASH .gbk directly into artspipeline1.py
# (no -ras/antiSMASH re-run).
#
# PHYLUM IS MANDATORY, NOT DEFAULTED. ARTS's core-gene comparison is
# phylum-specific (marker HMMs + gene matrix built per phylum reference
# set) -- running a genome against the wrong phylum's refdir produces no
# error, just silently meaningless/garbage core-gene results. This matters
# concretely for this project: the pseudoalteromonas_seq collection is
# taxonomically MIXED (most strains are Gammaproteobacteria, but ~48/143
# are Alphaproteobacteria -- Tritonibacter/Paracoccus/Loktanella), so the
# correct refdir must be decided per-genome, never assumed from which
# source directory/collection a genome came from.
#
# Since this script takes positional CLI arguments, submit it with the
# #BSUB flags given explicitly on the bsub command line (bsub < script.sh
# only parses the header block and cannot pass arguments), e.g.:
#
#   bsub -q hpc -n 4 -R "span[hosts=1] rusage[mem=8GB]" -M 8400MB -W 3:00 \
#       -J arts_<STRAIN> \
#       -o /work3/josne/github/arts/lsf/logs/arts_<STRAIN>_%J.out \
#       -e /work3/josne/github/arts/lsf/logs/arts_<STRAIN>_%J.err \
#       ./submit_one_genome.sh <STRAIN> <GBK_PATH> <PHYLUM> <RESULTDIR>
#
# <PHYLUM> must be a directory name under reference/ containing a
# prebuilt coremodels.hmm/genematrix.txt/model_metadata.json (e.g.
# gammaproteobacteria, alphaproteobacteria) -- checked below.
#==========================================================================

ARTS_DIR="/work3/josne/github/arts"

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <STRAIN> <GBK_PATH> <PHYLUM> <RESULTDIR>"
    echo ""
    echo "Available phylum reference sets under ${ARTS_DIR}/reference/:"
    for d in "${ARTS_DIR}"/reference/*/; do
        p="$(basename "$d")"
        if [ -f "${d}coremodels.hmm" ]; then
            echo "  ${p}"
        fi
    done
    exit 1
fi

STRAIN="$1"
INPUT="$2"
PHYLUM="$3"
RESULTDIR="$4"
CPU=4

REFDIR="${ARTS_DIR}/reference/${PHYLUM}/"

if [ ! -f "${INPUT}" ]; then
    echo "ERROR: Input GenBank file not found: ${INPUT}"
    exit 1
fi

for required_file in coremodels.hmm genematrix.txt model_metadata.json; do
    if [ ! -f "${REFDIR}${required_file}" ]; then
        echo "ERROR: '${PHYLUM}' is not a valid/complete phylum reference set"
        echo "       (missing ${REFDIR}${required_file})."
        echo ""
        echo "Available phylum reference sets under ${ARTS_DIR}/reference/:"
        for d in "${ARTS_DIR}"/reference/*/; do
            p="$(basename "$d")"
            if [ -f "${d}coremodels.hmm" ]; then
                echo "  ${p}"
            fi
        done
        exit 1
    fi
done

mkdir -p "${RESULTDIR}"
mkdir -p "${ARTS_DIR}/lsf/logs"

source "${ARTS_DIR}/lsf/setup.sh"

echo "=========================================="
echo "ARTS run - ${STRAIN}"
echo "Host       : $(hostname)"
echo "Job started: $(date)"
echo "LSF job ID : ${LSB_JOBID:-N/A}"
echo "Input      : ${INPUT}"
echo "Phylum     : ${PHYLUM}"
echo "Refdir     : ${REFDIR}"
echo "Resultdir  : ${RESULTDIR}"
echo "CPUs       : ${CPU}"
echo "=========================================="

cd "${ARTS_DIR}"

set +e
python artspipeline1.py \
    "${INPUT}" \
    "${REFDIR}" \
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
    echo "Sanity check FAILED - see above."
fi

exit "${RUN_STATUS}"
