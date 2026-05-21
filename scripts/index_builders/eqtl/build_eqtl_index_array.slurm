#!/bin/bash
# =============================================================================
# build_eqtl_index_array.slurm
#
# SLURM array job: runs build_eqtl_index.R once per ge dataset.
# Processes 10 datasets concurrently to stay within EBI FTP limits.
# =============================================================================

#SBATCH --job-name=build_eqtl_idx
#SBATCH --array=59,76,77,79,80,83,90,91,95,101,110,112,120,142,146,160,164,165,168,172,182,183,185,189,200,202,211,215,216,236,237,263,264,265,271,272,273,280
#SBATCH --time=05:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/log/build_eqtl_idx_%A_%a.out
#SBATCH --error=/data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/log/build_eqtl_idx_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengna.zhang@vumc.org

# ── Adjust these paths to match your installation ─────────────────────────────
PIPELINE_DIR="/data/h_vmac/zhanm32/colocpipe"
SCRIPT_DIR="${PIPELINE_DIR}/scripts/eqtl_index_build"
METADATA_FILE="${PIPELINE_DIR}/data/eQTLcatalogue/dataset_metadata_r7.tsv"
TABIX_PATHS_FILE="${PIPELINE_DIR}/data/eQTLcatalogue/tabix_ftp_paths.tsv"
OUT_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/eqtl_index_build"
P_THRESH="1e-4"
# ─────────────────────────────────────────────────────────────────────────────

# Load R — same module used by the main pipeline
set +euo pipefail
module purge 2>/dev/null || true
module load r/4.5.0 2>/dev/null || echo "[WARN] r/4.5.0 module not found"
set -euo pipefail

echo "============================================================"
echo " Array task  : ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_MAX}"
echo " Job ID      : ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo " Node        : $(hostname)"
echo " Started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

# Small random jitter (0-30 s) to stagger simultaneous FTP connections
JITTER=$(( RANDOM % 30 ))
echo " Sleeping ${JITTER}s to stagger FTP connections..."
sleep "${JITTER}"

Rscript "${SCRIPT_DIR}/build_eqtl_index.R" \
    --dataset_idx  "${SLURM_ARRAY_TASK_ID}" \
    --tabix_paths  "${TABIX_PATHS_FILE}"    \
    --metadata     "${METADATA_FILE}"       \
    --out_dir      "${OUT_DIR}"             \
    --p_thresh     "${P_THRESH}"            \
    --chunk_size   500000                   \
    --timeout_secs 3600

EXIT_CODE=$?

echo ""
echo " Finished    : $(date '+%Y-%m-%d %H:%M:%S')"
echo " Exit code   : ${EXIT_CODE}"
exit "${EXIT_CODE}"
