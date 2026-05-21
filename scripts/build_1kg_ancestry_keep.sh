#!/usr/bin/env bash
# =============================================================================
# build_1kg_ancestry_keep.sh
#
# One-time setup script. Reads the 1000 Genomes race/ancestry file and
# creates PLINK-format keep files (FID IID) for each superpopulation.
#
# These keep files are then used by query_1kg_ld.sh every time a regional
# LD matrix needs to be computed — no large pre-computed matrices required.
#
# Expected race file format (no header, space/tab delimited):
#   Col 1 (V1) : FID
#   Col 2 (V2) : IID
#   Col 3 (V3) : Superpopulation code (AFR / AMR / EAS / EUR / SAS)
#
# Run once:
#   bash build_1kg_ancestry_keep.sh \
#       --race   /data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38_race.txt \
#       --out    /data/h_vmac/zhanm32/colocpipe/data/1KG
#
# Arguments:
#   --race   Path to race/ancestry file                  [required]
#   --out    Output directory for keep files             [required]
#            Files created: EUR.keep AMR.keep AFR.keep
#            (EAS.keep and SAS.keep are also created)
#   --help   Show this message
#
# Output keep file format (PLINK2 --keep compatible):
#   FID  IID
#   0    HG00096
#   ...
#
# Sample sizes in 1KG GRCh38 freeze (for reference):
#   EUR: 503    AMR: 347    AFR: 661    EAS: 504    SAS: 489
# =============================================================================

set -euo pipefail

# ─── Help ─────────────────────────────────────────────────────────────────────
if [[ $# -eq 0 ]] || echo "$*" | grep -qw -- '--help\|-h'; then
    sed -n '/^# =/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
    exit 0
fi

# ─── Parse arguments ─────────────────────────────────────────────────────────
RACE_FILE=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --race)   RACE_FILE="$2"; shift 2 ;;
        --out)    OUT_DIR="$2";   shift 2 ;;
        --help|-h)
            sed -n '/^# =/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# ─── Validate ─────────────────────────────────────────────────────────────────
if [[ -z "$RACE_FILE" || -z "$OUT_DIR" ]]; then
    echo "[ERROR] --race and --out are both required."
    exit 1
fi
if [[ ! -f "$RACE_FILE" ]]; then
    echo "[ERROR] Race file not found: $RACE_FILE"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "============================================================"
echo " build_1kg_ancestry_keep.sh"
echo "============================================================"
echo " Race file   : $RACE_FILE"
echo " Output dir  : $OUT_DIR"
echo "============================================================"
echo ""

# ─── Detect population column ─────────────────────────────────────────────────
# The race file has no header. Determine which column holds the population code
# by checking which column contains known superpopulation codes.
# Expected: col1=FID col2=IID col3=POP  (but verify in case format differs)

echo "[1/2] Detecting file format..."

POP_COL=$(awk '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "AFR" || $i == "EUR" || $i == "AMR" || $i == "EAS" || $i == "SAS") {
            print i
            exit
        }
    }
}
' "$RACE_FILE")

if [[ -z "$POP_COL" ]]; then
    echo "[ERROR] Could not detect population column in the first row."
    echo "        Expected one of: AFR AMR EAS EUR SAS"
    echo "        First row: $(head -n1 "$RACE_FILE")"
    exit 1
fi

# FID and IID are assumed to be the two columns before the population column
FID_COL=$(( POP_COL - 2 ))
IID_COL=$(( POP_COL - 1 ))

if [[ "$FID_COL" -lt 1 || "$IID_COL" -lt 1 ]]; then
    # Population is in col 1 or 2 — assume IID only, set FID=0
    FID_COL=0
    IID_COL=1
    echo "        Population column: $POP_COL"
    echo "        Format: IID-only (FID will be set to 0)"
else
    echo "        FID column        : $FID_COL"
    echo "        IID column        : $IID_COL"
    echo "        Population column : $POP_COL"
fi

# ─── Create keep files ────────────────────────────────────────────────────────

echo ""
echo "[2/2] Creating keep files..."

for POP in EUR AMR AFR EAS SAS; do
    OUT_FILE="${OUT_DIR}/${POP}.keep"

    if [[ "$FID_COL" -eq 0 ]]; then
        awk -v pop="$POP" -v pc="$POP_COL" -v ic="$IID_COL" \
            '$pc == pop { print "0\t" $ic }' "$RACE_FILE" > "$OUT_FILE"
    else
        awk -v pop="$POP" -v pc="$POP_COL" -v fc="$FID_COL" -v ic="$IID_COL" \
            '$pc == pop { print $fc "\t" $ic }' "$RACE_FILE" > "$OUT_FILE"
    fi

    N=$(wc -l < "$OUT_FILE")
    printf "        %-6s : %4d samples -> %s\n" "$POP" "$N" "$OUT_FILE"
done

echo ""
echo "============================================================"
echo " Done! Keep files created in: $OUT_DIR"
echo "============================================================"
echo ""
echo " Next step — compute a regional LD matrix:"
echo "   bash query_1kg_ld.sh \\"
echo "       --ancestry  EUR \\"
echo "       --bfile     /data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38 \\"
echo "       --keep_dir  ${OUT_DIR} \\"
echo "       --chr       19 \\"
echo "       --start     44783684 \\"
echo "       --end       45033684 \\"
echo "       --out       /path/to/output/EUR_chr19_region"
