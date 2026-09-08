#!/bin/bash
### General options
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8400MB
#BSUB -W 5:00
#BSUB -J "arts_batch[1-438]%20"
#BSUB -o /work3/josne/github/arts/lsf/logs/arts_batch_%J_%I.out
#BSUB -e /work3/josne/github/arts/lsf/logs/arts_batch_%J_%I.err

set -euo pipefail

#==========================================================================
# ARTS full-batch LSF job array driver.
# One array task per row of arts_phylum_samplesheet.csv (built by
# build_phylum_samplesheet.py from KmerFinder taxonomy). Each task looks up
# its own row via $LSB_JOBINDEX and delegates to submit_one_genome.sh,
# which does the actual phylum validation + ARTS invocation + sanity check.
#
# %20 throttles concurrency to be considerate of the shared hpc queue.
# The [1-438] array size must match the samplesheet's row count -- update
# both if the samplesheet is regenerated with a different count.
#
# Submit with: bsub < submit_batch_array.sh
#==========================================================================

CSV="/work3/josne/Projects/Vibrio_Galathea3/arts_phylum_samplesheet.csv"
RESULTS_BASE="/work3/josne/Projects/Vibrio_Galathea3/arts_results/batch"
ARTS_DIR="/work3/josne/github/arts"

if [ -z "${LSB_JOBINDEX:-}" ]; then
    echo "ERROR: this script must be run as an LSF job array task (LSB_JOBINDEX unset)."
    exit 1
fi

# CSV row 1 (after header) corresponds to array index 1.
ROW=$(tail -n +2 "${CSV}" | sed -n "${LSB_JOBINDEX}p")
if [ -z "${ROW}" ]; then
    echo "ERROR: no CSV row for job index ${LSB_JOBINDEX} (CSV has $(tail -n +2 "${CSV}" | wc -l) data rows)."
    exit 1
fi

COLLECTION=$(echo "${ROW}" | cut -d',' -f1)
STRAIN=$(echo "${ROW}" | cut -d',' -f2)
PHYLUM=$(echo "${ROW}" | cut -d',' -f3)
GBK=$(echo "${ROW}" | cut -d',' -f5)

RESULTDIR="${RESULTS_BASE}/${PHYLUM}/${STRAIN}"

echo "=========================================="
echo "Batch array task ${LSB_JOBINDEX}"
echo "Collection : ${COLLECTION}"
echo "Strain     : ${STRAIN}"
echo "Phylum     : ${PHYLUM}"
echo "GBK        : ${GBK}"
echo "Resultdir  : ${RESULTDIR}"
echo "=========================================="

exec "${ARTS_DIR}/lsf/submit_one_genome.sh" "${STRAIN}" "${GBK}" "${PHYLUM}" "${RESULTDIR}"
