#!/usr/bin/env bash
# =============================================================================
# merge_lbf_geneindex.sh
#
# Run ONCE after the build_lbf_index.slurm array job completes.
# Concatenates all {dataset_id}.lbf.genelist.tsv files into a single
# lbf_gene_index.tsv that fetch_eqtl_region.R loads to decide whether
# to stream a dataset from EBI FTP.
#
# Usage:
#   bash scripts/merge_lbf_geneindex.sh [quant_method]
#
# Output:
#   data/eQTLcatalogue/lbf_gene_index.tsv   — header: dataset_id, gene_id
#
# fetch_eqtl_region.R logic:
#   Before streaming a dataset's LBF from FTP, it checks whether any of the
#   genes found in the nominal data appear in this index for that dataset.
#   If none match → skip FTP download entirely (gene was not fine-mapped).
#   If any match  → stream + filter via curl|zcat|awk as usual.
# =============================================================================

set -euo pipefail

PIPELINE_DIR="/data/h_vmac/zhanm32/colocpipe"
QUANT="${1:-ge}"

INDEX_DIR="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_index"
TASKLIST="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_tasklist_${QUANT}.tsv"
OUT_INDEX="${PIPELINE_DIR}/data/eQTLcatalogue/lbf_gene_index.tsv"

echo "============================================================"
echo " merge_lbf_geneindex.sh"
echo " Quant      : $QUANT"
echo " Index dir  : $INDEX_DIR"
echo " Output     : $OUT_INDEX"
echo " Started    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

[[ ! -f "$TASKLIST" ]] && { echo "[ERROR] Task list not found: $TASKLIST"; exit 1; }

TOTAL=$(tail -n +2 "$TASKLIST" | wc -l)
DONE=0
MISSING=()

# Check completion status
while IFS=$'\t' read -r ds_id lbf_url study_label; do
    DONE_FILE="${INDEX_DIR}/${ds_id}.lbf.genelist.done"
    if [[ -f "$DONE_FILE" ]]; then
        (( ++DONE ))
    else
        MISSING+=("$ds_id")
    fi
done < <(tail -n +2 "$TASKLIST")

echo " Build status: $DONE / $TOTAL complete"

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo ""
    echo " [WARN] ${#MISSING[@]} dataset(s) not yet done:"
    for ds in "${MISSING[@]}"; do echo "   $ds"; done
    echo ""
    echo " Merging only completed datasets. Re-run after remaining tasks finish."
fi
echo ""

# Write header then concatenate all available gene lists
echo -e "dataset_id\tgene_id" > "$OUT_INDEX"

N_FILES=0
N_ROWS=0
while IFS=$'\t' read -r ds_id lbf_url study_label; do
    GENE_LIST="${INDEX_DIR}/${ds_id}.lbf.genelist.tsv"
    if [[ -f "$GENE_LIST" && -s "$GENE_LIST" ]]; then
        cat "$GENE_LIST" >> "$OUT_INDEX"
        ROWS=$(wc -l < "$GENE_LIST")
        N_ROWS=$(( N_ROWS + ROWS ))
        (( ++N_FILES ))
    fi
done < <(tail -n +2 "$TASKLIST")

TOTAL_ROWS=$(( N_ROWS + 1 ))   # +1 for header

echo " Gene lists merged : $N_FILES files"
echo " Total gene-dataset pairs : $N_ROWS"
echo " Output : $OUT_INDEX  ($(numfmt --to=iec --suffix=B $(stat -c%s "$OUT_INDEX")))"
echo ""
echo " Next step — update colocpipe to use this index."
echo " fetch_eqtl_region.R will read it via --lbf_gene_index"
echo ""
echo " To clean up the existing large .lbf.tsv.gz files from the"
echo " previous indexing strategy (if any), run:"
echo "   rm -f ${INDEX_DIR}/*.lbf.tsv.gz ${INDEX_DIR}/*.lbf.tsv.gz.csi"
echo ""
echo " Finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
