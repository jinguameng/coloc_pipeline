#!/usr/bin/env bash
# =============================================================================
# check_lbf_index.sh
#
# Reports build progress and prints a --array= range for any failed tasks
# so they can be resubmitted without re-running successful ones.
#
# Sentinel file checked: {dataset_id}.lbf.genelist.done
# After all tasks complete, run: bash scripts/merge_lbf_geneindex.sh
#
# Usage:
#   bash scripts/check_lbf_index.sh [quant_method]
#
# Example resubmit output:
#   sbatch --array=3,17,45 scripts/build_lbf_index.slurm
# =============================================================================

set -euo pipefail

PIPELINE_DIR="/data/h_vmac/zhanm32/colocpipe"
QUANT="${1:-ge}"

TASKLIST="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_tasklist_${QUANT}.tsv"
INDEX_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_index"

[[ ! -f "$TASKLIST" ]] && { echo "[ERROR] Task list not found: $TASKLIST"; exit 1; }

TOTAL=$(tail -n +2 "$TASKLIST" | wc -l)
DONE=0
MISSING=()

while IFS=$'\t' read -r ds_id lbf_url study_label; do
    DONE_FILE="${INDEX_DIR}/${ds_id}.lbf.genelist.done"
    if [[ -f "$DONE_FILE" ]]; then
        (( DONE++ ))
    else
        MISSING+=("$ds_id")
    fi
done < <(tail -n +2 "$TASKLIST")

FAIL=$(( TOTAL - DONE ))

echo "============================================================"
echo " LBF index build status (quant=$QUANT)"
echo "============================================================"
echo " Total datasets : $TOTAL"
echo " Done           : $DONE"
echo " Remaining      : $FAIL"
echo ""

if [[ ${#MISSING[@]} -eq 0 ]]; then
    echo " All datasets successfully indexed."
    echo " fetch_eqtl_region.R will use local tabix queries for all datasets."
    exit 0
fi

echo " Missing datasets:"
for ds in "${MISSING[@]}"; do
    echo "   $ds"
done

echo ""
echo " Task IDs to resubmit:"
# Convert dataset names back to task IDs (line numbers in the task list)
TASK_IDS=()
while IFS=$'\t' read -r ds_id lbf_url study_label; do
    for missing in "${MISSING[@]}"; do
        if [[ "$ds_id" == "$missing" ]]; then
            # task ID = line number in file (header = line 1, first data = line 2 = task 1)
            TASK_IDS+=("$(grep -n "^${ds_id}"$'\t' "$TASKLIST" | cut -d: -f1 | awk '{print $1-1}')")
        fi
    done
done < <(tail -n +2 "$TASKLIST")

ARRAY_STR=$(IFS=,; echo "${TASK_IDS[*]}")
echo ""
echo "   sbatch --array=${ARRAY_STR}%8 ${PIPELINE_DIR}/scripts/build_lbf_index.slurm"
