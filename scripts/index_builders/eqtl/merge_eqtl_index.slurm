#!/bin/bash
# =============================================================================
# merge_eqtl_index.slurm
#
# Runs AFTER the build array job completes. Merges all per-dataset TSV parts
# into a single sorted, bgzipped, tabix-indexed file ready for fast lookups.
#
# Final output:
#   eqtl_index_ge_p1e4.tsv.gz       — bgzipped index
#   eqtl_index_ge_p1e4.tsv.gz.tbi   — tabix index
#
# Column order in final index (no chr prefix on chromosome values):
#   chromosome  position  dataset_id  study_label  tissue_label
#   gene_id  pvalue  beta  maf
# =============================================================================

#SBATCH --job-name=merge_eqtl_idx
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/log/merge_eqtl_idx_%j.out
#SBATCH --error=/data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/log/merge_eqtl_idx_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mengna.zhang@vumc.org

# ── Adjust these paths ────────────────────────────────────────────────────────
PIPELINE_DIR="/data/h_vmac/zhanm32/colocpipe"
OUT_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/eqtl_index_build"
FINAL_DIR="${PIPELINE_DIR}/data/eQTLcatalogue"
FINAL_NAME="eqtl_index_ge_p1e4"
# ─────────────────────────────────────────────────────────────────────────────

PARTS_DIR="${OUT_DIR}/parts"
STATUS_DIR="${OUT_DIR}/status"
MERGED_RAW="${OUT_DIR}/${FINAL_NAME}_unsorted.tsv"
MERGED_SORTED="${OUT_DIR}/${FINAL_NAME}_sorted.tsv"
FINAL_GZ="${FINAL_DIR}/${FINAL_NAME}.tsv.gz"

set +euo pipefail
module purge 2>/dev/null || true
module load r/4.5.0 2>/dev/null || echo "[WARN] r/4.5.0 not found"
module load htslib 2>/dev/null || \
echo "[WARN] htslib module not found — bgzip and tabix must be on PATH"
set -euo pipefail

echo "============================================================"
echo " merge_eqtl_index.slurm"
echo " Node     : $(hostname)"
echo " Started  : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

# ── 1. Report build outcomes ──────────────────────────────────────────────────
N_DONE=$(ls "${STATUS_DIR}"/*.done   2>/dev/null | wc -l)

N_PARTS=$(find "${PARTS_DIR}" -name "*.tsv" -size +0c 2>/dev/null | wc -l)

echo ""
echo "[1/5] Build array summary"
echo "      Tasks succeeded (have .done)  : ${N_DONE}"
echo "      Part files with data           : ${N_PARTS}"
echo ""

if [[ "${N_PARTS}" -eq 0 ]]; then
    echo "[ERROR] No part files found in ${PARTS_DIR} — aborting."
    exit 1
fi

# ── 2. Concatenate part files (no headers — parts were written without them) ──
echo "[2/5] Concatenating ${N_PARTS} part files..."

HEADER="#chromosome\tposition\tdataset_id\tstudy_label\ttissue_label\tgene_id\tpvalue\tbeta\tmaf"

# Write header then all non-empty part files
printf "%b\n" "${HEADER}" > "${MERGED_RAW}"

# Find all non-empty part files and cat them
find "${PARTS_DIR}" -name "*.tsv" -size +0c | sort | xargs cat >> "${MERGED_RAW}"

N_TOTAL=$(wc -l < "${MERGED_RAW}")
echo "      Total rows (incl. header): ${N_TOTAL}"

# ── 3. Sort by chromosome (numeric), then position ────────────────────────────
echo "[3/5] Sorting by chromosome and position..."

# Separate header from data, sort numerically, then reattach header.
# The awk step pads chromosome numbers for correct numeric ordering:
#   1 → 01, 2 → 02, ..., 9 → 09, 10 → 10, X → 23, Y → 24
# then the sort key is removed and original chromosome restored.
{
    head -1 "${MERGED_RAW}"          # keep header line first
    tail -n +2 "${MERGED_RAW}" |    # skip header
    awk 'BEGIN{OFS="\t"} {
        chr=$1
        if      (chr=="X") sortkey=23
        else if (chr=="Y") sortkey=24
        else              sortkey=int(chr)+0
        print sortkey, $0
    }' |
    sort --parallel=4 -k1,1n -k3,3n |   # sort on numeric chr key, then position
    cut -f2-                               # drop the sort key column
} > "${MERGED_SORTED}"

echo "      Sort complete"
rm -f "${MERGED_RAW}"

# ── 4. bgzip ──────────────────────────────────────────────────────────────────
echo "[4/5] Compressing with bgzip..."

bgzip --threads 4 --stdout "${MERGED_SORTED}" > "${FINAL_GZ}"
rm -f "${MERGED_SORTED}"

echo "      Written: ${FINAL_GZ}"
echo "      Size   : $(du -sh "${FINAL_GZ}" | cut -f1)"

# ── 5. Tabix index ────────────────────────────────────────────────────────────
# -s 1  : chromosome column is column 1
# -b 2  : begin position is column 2
# -e 2  : end position is column 2 (SNP-level, start == end)
# Lines starting with '#' are treated as header/comments automatically
echo "[5/5] Building tabix index..."

tabix -s 1 -b 2 -e 2 "${FINAL_GZ}"

echo "      Index : ${FINAL_GZ}.tbi"

# ── Summary ───────────────────────────────────────────────────────────────────
N_DATA=$(zcat "${FINAL_GZ}" | grep -v "^#" | wc -l)
N_GENES=$(zcat "${FINAL_GZ}" | grep -v "^#" | cut -f6 | sort -u | wc -l)
N_DS=$(zcat "${FINAL_GZ}" | grep -v "^#"   | cut -f3 | sort -u | wc -l)

echo ""
echo "============================================================"
echo " Index build complete"
echo "   File          : ${FINAL_GZ}"
echo "   Total rows    : ${N_DATA}"
echo "   Unique genes  : ${N_GENES}"
echo "   Datasets      : ${N_DS} / ${N_DONE} successfully processed"
echo "   Finished      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
echo "Next step — query the index:"
echo "  Rscript scripts/query_eqtl_index.R \\"
echo "      --index ${FINAL_GZ} \\"
echo "      --chr 19 --start 44658684 --end 45158684"
