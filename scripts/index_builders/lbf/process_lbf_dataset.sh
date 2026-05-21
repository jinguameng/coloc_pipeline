#!/usr/bin/env bash
# =============================================================================
# process_lbf_dataset.sh
#
# Downloads ONE eQTL Catalogue LBF file, extracts the list of fine-mapped
# gene IDs, then DELETES the full file. The output is a tiny gene list
# (~KB) rather than a full sorted+indexed LBF file (~GB).
#
# This lightweight index is used by fetch_eqtl_region.R to decide whether
# to stream a dataset from EBI FTP: if none of the requested genes appear
# in the gene list, the FTP download is skipped entirely.
#
# Manual test:
#   bash scripts/process_lbf_dataset.sh \
#       QTD000001 \
#       "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000001/QTD000001.lbf_variable.txt.gz" \
#       /data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/lbf_index
#
# Arguments:
#   $1  dataset_id
#   $2  lbf_url      FTP URL to .lbf_variable.txt.gz
#   $3  index_dir    destination for output gene list
#
# Optional env vars:
#   WORK_DIR   download location (default: index_dir/tmp_$dataset_id)
#   KEEP_RAW   set to "1" to keep raw download (for debugging)
#   N_CPUS     CPUs for sort (default: 2)
#
# Outputs:
#   {index_dir}/{dataset_id}.lbf.genelist.tsv   2 cols: dataset_id, gene_id
#   {index_dir}/{dataset_id}.lbf.genelist.done  empty sentinel on success
#
# After ALL tasks complete, run merge_lbf_geneindex.sh to produce the
# single lbf_gene_index.tsv consumed by fetch_eqtl_region.R.
# =============================================================================

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: bash process_lbf_dataset.sh <dataset_id> <lbf_url> <index_dir>"
    exit 1
fi

DS_ID="$1"
LBF_URL="$2"
INDEX_DIR="$3"
N_CPUS="${N_CPUS:-2}"
KEEP_RAW="${KEEP_RAW:-0}"

WORK_DIR="${WORK_DIR:-${INDEX_DIR}/tmp_${DS_ID}}"
mkdir -p "$WORK_DIR" "$INDEX_DIR"

RAW_GZ="${WORK_DIR}/${DS_ID}.lbf_variable.txt.gz"
GENE_LIST="${INDEX_DIR}/${DS_ID}.lbf.genelist.tsv"
DONE_FILE="${INDEX_DIR}/${DS_ID}.lbf.genelist.done"

cleanup() {
    if [[ "$KEEP_RAW" == "1" ]]; then
        echo "[CLEANUP] KEEP_RAW=1 — leaving: $WORK_DIR"
    else
        echo "[CLEANUP] Removing: $WORK_DIR"
        rm -rf "$WORK_DIR"
    fi
}
trap cleanup EXIT

echo "============================================================"
echo " process_lbf_dataset.sh"
echo " Dataset   : $DS_ID"
echo " URL       : $LBF_URL"
echo " Work dir  : $WORK_DIR"
echo " Gene list : $GENE_LIST"
echo " Keep raw  : $KEEP_RAW"
echo " Started   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

if [[ -f "$DONE_FILE" ]]; then
    echo "[SKIP] Already done ($DONE_FILE exists)"
    exit 0
fi

for tool in curl zcat awk sort; do
    command -v "$tool" &>/dev/null || { echo "[ERROR] Not on PATH: $tool"; exit 1; }
done
echo "[OK] Tools: curl zcat awk sort"
echo ""

# ─── STEP 1: Download ─────────────────────────────────────────────────────────
echo "[1/3] Downloading..."
echo "      URL  : $LBF_URL"
AVAIL_KB=$(df -k "$WORK_DIR" | awk 'NR==2{print $4}')
echo "      Available space : $(numfmt --to=iec --suffix=B $((AVAIL_KB * 1024)))"

curl -f --retry 5 --retry-delay 10 --retry-max-time 600 \
     --progress-bar -o "$RAW_GZ" "$LBF_URL"

RAW_BYTES=$(stat -c%s "$RAW_GZ")
echo "      Downloaded : $(numfmt --to=iec --suffix=B $RAW_BYTES)"

# ─── STEP 2: Inspect header ───────────────────────────────────────────────────
echo ""
echo "[2/3] Inspecting header..."

# set +o pipefail: awk 'NR==1{exit}' closes stdin early, sending SIGPIPE to
# zcat (exit 141). pipefail would treat that as fatal. Restore immediately.
set +o pipefail
HEADER=$(zcat "$RAW_GZ" | awk 'NR==1{print; exit}')
set -o pipefail

[[ -z "$HEADER" ]] && { echo "[ERROR] Empty header — file corrupt?"; exit 1; }

N_COLS=$(echo "$HEADER" | awk -F'\t' '{print NF}')
echo "      Columns : $N_COLS"
echo "      First 6 : $(echo "$HEADER" | cut -f1-6)"

GENE_COL=$(echo "$HEADER" | tr '\t' '\n' | awk '/^molecular_trait_id$/{print NR; exit}')
CHR_COL=$(echo  "$HEADER" | tr '\t' '\n' | awk '/^chromosome$/{print NR; exit}')
POS_COL=$(echo  "$HEADER" | tr '\t' '\n' | awk '/^position$/{print NR; exit}')

[[ -z "$GENE_COL" ]] && { echo "[ERROR] 'molecular_trait_id' column missing. Header: $HEADER"; exit 1; }
[[ -z "$CHR_COL"  ]] && { echo "[ERROR] 'chromosome' column missing";  exit 1; }
[[ -z "$POS_COL"  ]] && { echo "[ERROR] 'position' column missing";    exit 1; }

echo "      molecular_trait_id : col $GENE_COL"
echo "      chromosome         : col $CHR_COL"
echo "      position           : col $POS_COL"

# ─── STEP 3: Extract gene list, delete raw file ───────────────────────────────
echo ""
echo "[3/3] Extracting gene list..."

# awk reads to EOF (no early exit) — no SIGPIPE under pipefail.
# One output row per unique (dataset_id, gene_id) pair.
zcat "$RAW_GZ" \
    | awk -F'\t' -v ds="$DS_ID" -v gc="$GENE_COL" \
        'NR>1 && $gc != "" { print ds"\t"$gc }' \
    | sort -u \
    > "$GENE_LIST"

N_GENES=$(wc -l < "$GENE_LIST")
echo "      Unique genes : $N_GENES"

[[ "$N_GENES" -eq 0 ]] && echo "[WARN] No genes extracted — dataset may have no fine-mapped signals"

[[ "$KEEP_RAW" != "1" ]] && { rm -f "$RAW_GZ"; echo "      Raw download deleted"; }

touch "$DONE_FILE"

echo ""
echo "============================================================"
echo " Done!"
echo "   Gene list : $GENE_LIST  ($N_GENES genes)"
echo "   Finished  : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
