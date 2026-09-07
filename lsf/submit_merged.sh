#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -J arts_smoketest_merged
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8400MB
#BSUB -W 5:00
#BSUB -o /work3/josne/github/arts/lsf/logs/arts_smoketest_merged_%J.out
#BSUB -e /work3/josne/github/arts/lsf/logs/arts_smoketest_merged_%J.err

set -euo pipefail

#==========================================================================
# ARTS smoke test — merged run (S0204 + S1608), single shared refdir.
# ARTS processes comma-separated genomes strictly sequentially internally,
# then combine_results.py merges per-genome tables — walltime budget is
# roughly 2x the single-genome benchmark plus margin for the combine step.
# Both genomes are Gammaproteobacteria, so a single refdir is valid.
#==========================================================================
INPUT_VIBRIO="/work3/josne/Projects/Vibrio_Galathea3/vibrio_seq/funcscan_results/bgc/antismash/S0204/S0204.gbk"
INPUT_PSEUDOALT="/work3/josne/Projects/Vibrio_Galathea3/pseudoalteromonas_seq/funcscan_results/bgc/antismash/S1608/S1608.gbk"
RESULTDIR="/work3/josne/Projects/Vibrio_Galathea3/arts_results/smoketest/merged"
CPU=4
#==========================================================================

ARTS_DIR="/work3/josne/github/arts"

for f in "${INPUT_VIBRIO}" "${INPUT_PSEUDOALT}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Input GenBank file not found: ${f}"
        exit 1
    fi
done

# ARTS takes a single comma-separated (no spaces) input list for multi-genome
# combined mode; both genomes must share one refdir.
INPUT="${INPUT_VIBRIO},${INPUT_PSEUDOALT}"

mkdir -p "${RESULTDIR}"
mkdir -p "${ARTS_DIR}/lsf/logs"

source "${ARTS_DIR}/lsf/setup.sh"

echo "=========================================="
echo "ARTS smoke test - Merged (S0204 + S1608)"
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
# In multi-genome mode, ARTS writes each genome's own results under
# <RESULTDIR>/<basename(RESULTDIR)>_<N>/tables/ (N=0,1,... in input order),
# NOT directly under <RESULTDIR>/tables/. combine_results.py separately
# writes combined_core_table.*, combined_known_table.*, combined_dup_table.*
# under <RESULTDIR>/tables/ itself.
echo ""
echo "--- Sanity check: output tables ---"
RESULTDIR_BASENAME="$(basename "${RESULTDIR}")"
CORETABLE_0="${RESULTDIR}/${RESULTDIR_BASENAME}_0/tables/coretable.tsv"
CORETABLE_1="${RESULTDIR}/${RESULTDIR_BASENAME}_1/tables/coretable.tsv"
BGCTABLE_0="${RESULTDIR}/${RESULTDIR_BASENAME}_0/tables/bgctable.tsv"
BGCTABLE_1="${RESULTDIR}/${RESULTDIR_BASENAME}_1/tables/bgctable.tsv"
COMBINED_CORE="${RESULTDIR}/tables/combined_core_table.tsv"

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
check_table "${CORETABLE_0}" || SANITY_OK=0
check_table "${CORETABLE_1}" || SANITY_OK=0
check_table "${BGCTABLE_0}" || SANITY_OK=0
check_table "${BGCTABLE_1}" || SANITY_OK=0
check_table "${COMBINED_CORE}" || SANITY_OK=0

if [ "${SANITY_OK}" -eq 1 ]; then
    echo "Sanity check PASSED"
else
    echo "Sanity check FAILED - see above."
fi

exit "${RUN_STATUS}"
