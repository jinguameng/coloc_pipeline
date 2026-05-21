#!/usr/bin/env bash
# =============================================================================
# setup_lbf_tasklist.sh
#
# One-time setup: joins dataset_metadata_r7.tsv with tabix_ftp_paths.tsv,
# filters to datasets that have a non-empty ftp_lbf_path, and writes a
# task list used by build_lbf_index.slurm.
#
# Run ONCE from the colocpipe root directory before submitting the array job:
#   bash scripts/setup_lbf_tasklist.sh [quant_method]
#
# Arguments:
#   quant_method  — ge | exon | tx | txrev | all (default: ge)
#
# Outputs:
#   data/eQTLcatalogue/lbf_tasklist.tsv   — 3 columns: dataset_id, lbf_url, study_label
#   data/eQTLcatalogue/lbf_index/         — output directory (created, initially empty)
#
# After running this script, submit the array job:
#   N=$(tail -n +2 data/eQTLcatalogue/lbf_tasklist.tsv | wc -l)
#   sbatch --array=1-${N}%8 scripts/build_lbf_index.slurm
# =============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUANT="${1:-ge}"
QUANT="${QUANT,,}"

METADATA="${PIPELINE_DIR}/data/eQTLcatalogue/dataset_metadata_r7.tsv"
FTP_PATHS="${PIPELINE_DIR}/data/eQTLcatalogue/tabix_ftp_paths.tsv"
TASKLIST="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_tasklist${QUANT:+_${QUANT}}.tsv"
INDEX_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_index"

echo "============================================================"
echo " setup_lbf_tasklist.sh"
echo " Quant method : ${QUANT}"
echo " Metadata     : ${METADATA}"
echo " FTP paths    : ${FTP_PATHS}"
echo " Task list    : ${TASKLIST}"
echo " Index dir    : ${INDEX_DIR}"
echo "============================================================"

[[ ! -f "$METADATA"  ]] && { echo "[ERROR] Metadata file not found: $METADATA";  exit 1; }
[[ ! -f "$FTP_PATHS" ]] && { echo "[ERROR] FTP paths file not found: $FTP_PATHS"; exit 1; }

mkdir -p "$INDEX_DIR"

# ─── Join metadata + FTP paths, filter, emit task list ───────────────────────
# Uses awk to:
#   1. Read FTP paths file into a map: dataset_id → (lbf_url, study_label)
#   2. Read metadata, filter by quant_method, join on dataset_id
#   3. Skip rows where lbf_url is empty or "NA"
# Output: TSV with header: dataset_id | lbf_url | study_label

python3 - <<PYEOF
import csv, sys

quant = "${QUANT}"

# Load FTP paths: dataset_id → lbf_url
ftp = {}
with open("${FTP_PATHS}") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        ds  = row.get("dataset_id","").strip()
        url = row.get("ftp_lbf_path","").strip()
        if ds and url and url.lower() not in ("", "na"):
            ftp[ds] = url

# Load metadata, filter by quant_method, join
rows = []
with open("${METADATA}") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        ds = row.get("dataset_id","").strip()
        qm = row.get("quant_method","").strip().lower()
        sl = row.get("study_label","").strip()
        if quant != "all" and qm != quant:
            continue
        if ds not in ftp:
            continue
        rows.append((ds, ftp[ds], sl))

# Write
with open("${TASKLIST}", "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["dataset_id", "lbf_url", "study_label"])
    for r in rows:
        writer.writerow(r)

print(f"  Written {len(rows)} tasks to ${TASKLIST}")
PYEOF

N=$(tail -n +2 "$TASKLIST" | wc -l)
echo ""
echo " Task list written : $TASKLIST"
echo " Total tasks       : $N"
echo " Index directory   : $INDEX_DIR"
echo ""
echo " Next step — submit the array job:"
echo "   sbatch --array=1-${N}%8 ${PIPELINE_DIR}/scripts/build_lbf_index.slurm"
echo ""
echo " (The %8 limits concurrent tasks to 8 to respect EBI FTP rate limits."
echo "  Increase to %20 if EBI does not complain; decrease to %4 if you see"
echo "  frequent curl errors.)"
