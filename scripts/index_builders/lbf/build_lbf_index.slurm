#!/usr/bin/env bash
# =============================================================================
# build_lbf_index.slurm
#
# Thin SLURM array wrapper around process_lbf_dataset.sh.
# Each task downloads one LBF file, extracts a gene list, and deletes the
# full file — leaving only a tiny {dataset_id}.lbf.genelist.tsv (~KB).
#
# After ALL tasks complete, merge the gene lists:
#   bash scripts/merge_lbf_geneindex.sh
#
# Manual test first:
#   module load htslib
#   LINE=$(awk 'NR==2' data/eQTLcatalogue/lbf_tasklist_ge.tsv)
#   DS=$(echo "$LINE" | cut -f1)
#   URL=$(echo "$LINE" | cut -f2)
#   KEEP_RAW=1 \
#   WORK_DIR=/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/lbf_index/tmp_test \
#   bash scripts/process_lbf_dataset.sh "$DS" "$URL" \
#       /data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/lbf_index
#
# Submit:
#   N=$(tail -n +2 data/eQTLcatalogue/lbf_tasklist_ge.tsv | wc -l)
#   sbatch --array=1-${N}%10 scripts/build_lbf_index.slurm
# =============================================================================

#SBATCH --job-name=build_lbf_idx
#SBATCH --array=11-280%10
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/h_vmac/zhanm32/colocpipe/log/build_lbf_%A_%a.out
#SBATCH --error=/data/h_vmac/zhanm32/colocpipe/log/build_lbf_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengna.zhang@vumc.org

set -euo pipefail

PIPELINE_DIR="/data/h_vmac/zhanm32/colocpipe"
QUANT="ge"

TASKLIST="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_tasklist_${QUANT}.tsv"
INDEX_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_index"
SCRIPT="${PIPELINE_DIR}/scripts/lbf_index_build/process_lbf_dataset.sh"


# ── Read task ─────────────────────────────────────────────────────────────────
LINE=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR == n+1 { print; exit }' "$TASKLIST")
[[ -z "$LINE" ]] && { echo "[ERROR] No entry for task $SLURM_ARRAY_TASK_ID"; exit 1; }

DS_ID=$(echo   "$LINE" | cut -f1)
LBF_URL=$(echo "$LINE" | cut -f2)

echo "Job ${SLURM_ARRAY_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}: $DS_ID"

# ── Use panfs instead of /tmp ─────────────────────────────────────────────────
# /tmp on compute nodes is a small ramdisk — large LBF files overflow it
# silently (no stderr). Stay on panfs for the work directory.
export WORK_DIR="${INDEX_DIR}/tmp_${DS_ID}_${SLURM_ARRAY_JOB_ID}"
export N_CPUS="${SLURM_CPUS_PER_TASK:-2}"

# ── Clean up leftover manual test directory ───────────────────────────────────
TMP_TEST="${INDEX_DIR}/tmp_test"
if [[ -d "$TMP_TEST" ]]; then
    echo "[CLEANUP] Removing manual test directory: $TMP_TEST"
    rm -rf "$TMP_TEST"
fi

# ── Safety trap: remove WORK_DIR on any exit ──────────────────────────────────
trap 'echo "[SLURM CLEANUP] Removing: $WORK_DIR"; rm -rf "$WORK_DIR"' EXIT

# ── Delegate to process_lbf_dataset.sh ───────────────────────────────────────
bash "$SCRIPT" "$DS_ID" "$LBF_URL" "$INDEX_DIR"
