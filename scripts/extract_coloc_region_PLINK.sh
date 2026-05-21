#!/bin/bash
# =============================================================================
# extract_coloc_region_PLINK.sh
#
# Purpose: Extract a genomic region around a lead SNP from PLINK 1.9
#          association results (.assoc.linear or .assoc.logistic), and
#          prepare the output for coloc (SNP, CHR, BP, A1, BETA, VARBETA, P).
#
# Handles:
#   --linear  output: BETA used directly;    VARBETA = SE^2
#   --logistic output: OR -> BETA = log(OR); VARBETA = SE^2
#              (PLINK SE for logistic = SE of log(OR), not SE of OR)
#
# Note on A2 (non-effect allele):
#   PLINK --linear / --logistic output only contains A1, so A2 is not
#   included here. Allele alignment for coloc is handled downstream in
#   prepare_coloc_datasets.R, which derives A2 from the eQTL Catalogue
#   ref/alt columns at the matching chr:pos position.
#
# Requirements:
#   - PLINK output must include SE column (run PLINK with --ci 0.95)
#   - Only ADD (additive) test rows are kept; covariate rows are dropped
#   - Rows with NA in BETA/OR, SE, or P are silently dropped
#
# Usage:
#   bash extract_coloc_region_PLINK.sh \
#       -f <assoc_file> \
#       -s <snp_id> \
#       [-w <window_kb>] \
#       [-o <output_file>]
#
# Arguments:
#   -f   Path to PLINK association results file  (required)
#   -s   rsID of the lead SNP                   (required)
#   -w   Total window size in KB, centred on lead SNP  (default: 250)
#   -o   Output file path  (default: <snp_id>_coloc_<window>kb.txt)
#   -h   Show this help message
#
# Output columns (tab-delimited):
#   SNP  CHR  BP  A1  BETA  VARBETA  P  N
#
#   N is populated from the NMISS column (number of non-missing observations),
#   which PLINK always emits for --linear / --logistic output. The downstream
#   prepare_coloc_datasets.R script takes median(N) across the region, matching
#   the same approach used for GWAMA input. No user-supplied N is required.
#   If NMISS is absent (unusual), N is written as 0 and the user-provided
#   fallback in inputs.txt is used instead.
#
# Example:
#   bash extract_coloc_region_PLINK.sh \
#       -f /data/.../ADNI_EUR.SPARE_AD.assoc.linear \
#       -s rs429358 \
#       -w 250 \
#       -o /data/.../rs429358_region250KB.txt
# =============================================================================

set -euo pipefail

# ─── Defaults ────────────────────────────────────────────────────────────────

WINDOW_KB=250
OUTPUT=""
ASSOC_FILE=""
TARGET_SNP=""

# ─── Parse arguments ─────────────────────────────────────────────────────────

while getopts "f:s:w:o:h" opt; do
    case $opt in
        f) ASSOC_FILE="$OPTARG" ;;
        s) TARGET_SNP="$OPTARG" ;;
        w) WINDOW_KB="$OPTARG"  ;;
        o) OUTPUT="$OPTARG"     ;;
        h)
            sed -n '/^# Usage/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "[ERROR] Unknown option. Use -h for help." >&2; exit 1 ;;
    esac
done

# ─── Validate required inputs ────────────────────────────────────────────────

if [[ -z "$ASSOC_FILE" || -z "$TARGET_SNP" ]]; then
    echo "[ERROR] -f (assoc file) and -s (target SNP) are both required." >&2
    exit 1
fi
if [[ ! -f "$ASSOC_FILE" ]]; then
    echo "[ERROR] Assoc file not found: $ASSOC_FILE" >&2
    exit 1
fi
if ! [[ "$WINDOW_KB" =~ ^[0-9]+$ ]] || [[ "$WINDOW_KB" -le 0 ]]; then
    echo "[ERROR] -w must be a positive integer (KB). Got: $WINDOW_KB" >&2
    exit 1
fi

# ─── Set default output name ─────────────────────────────────────────────────

if [[ -z "$OUTPUT" ]]; then
    OUTPUT="${TARGET_SNP}_coloc_${WINDOW_KB}kb.txt"
fi

# ─── Calculate half-window in base pairs ─────────────────────────────────────

HALF_WINDOW=$(( WINDOW_KB * 1000 / 2 ))

echo "============================================================"
echo " extract_coloc_region_PLINK.sh"
echo "============================================================"
echo " Assoc file  : $ASSOC_FILE"
echo " Target SNP  : $TARGET_SNP"
echo " Window      : ±${HALF_WINDOW} bp  (${WINDOW_KB} KB total)"
echo " Output file : $OUTPUT"
echo "============================================================"
echo ""

# ─── Step 1: Detect format and validate header ───────────────────────────────

echo "[1/3] Detecting file format from header..."

HEADER=$(head -n 1 "$ASSOC_FILE")

for col in CHR SNP BP A1 TEST SE P; do
    if ! echo "$HEADER" | grep -qw "$col"; then
        echo "[ERROR] Required column '$col' not found in header." >&2
        echo "        Header: $HEADER" >&2
        [[ "$col" == "SE" ]] && \
            echo "        HINT: SE requires running PLINK with --ci 0.95" >&2
        exit 1
    fi
done

# NMISS is expected but not mandatory — N column will be 0 if absent
if echo "$HEADER" | grep -qw "NMISS"; then
    echo "        NMISS column : found (will be written as N; median used by coloc)"
else
    echo "        NMISS column : NOT found — N will be 0 (provide N manually in inputs.txt)"
fi

if echo "$HEADER" | grep -qw "BETA"; then
    MODE="linear"
    echo "        Format detected: LINEAR (--linear)"
elif echo "$HEADER" | grep -qw "OR"; then
    MODE="logistic"
    echo "        Format detected: LOGISTIC (--logistic)"
else
    echo "[ERROR] Neither BETA nor OR column found in header." >&2
    exit 1
fi

# ─── Step 2: Find the target SNP's position ──────────────────────────────────

echo ""
echo "[2/3] Locating target SNP '$TARGET_SNP'..."

SNP_INFO=$(awk -v snp="$TARGET_SNP" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "SNP")  col_snp  = i
            if ($i == "CHR")  col_chr  = i
            if ($i == "BP")   col_bp   = i
            if ($i == "TEST") col_test = i
        }
        next
    }
    $col_snp == snp && $col_test == "ADD" {
        print $col_chr "\t" $col_bp
        exit
    }
' "$ASSOC_FILE")

if [[ -z "$SNP_INFO" ]]; then
    echo "[ERROR] SNP '$TARGET_SNP' not found in file (searched TEST=ADD rows)." >&2
    echo "        Check that the rsID matches exactly (case-sensitive)." >&2
    exit 1
fi

TARGET_CHR=$(echo "$SNP_INFO" | cut -f1)
TARGET_BP=$(echo  "$SNP_INFO" | cut -f2)
REGION_START=$(( TARGET_BP - HALF_WINDOW ))
REGION_END=$(( TARGET_BP + HALF_WINDOW ))
[[ "$REGION_START" -lt 1 ]] && REGION_START=1

echo "        Found   : CHR${TARGET_CHR}:${TARGET_BP}"
echo "        Region  : CHR${TARGET_CHR}:${REGION_START}-${REGION_END}"

# ─── Step 3: Extract region, compute BETA/VARBETA, write output ──────────────

echo ""
echo "[3/3] Extracting region and computing BETA/VARBETA..."

awk \
    -v mode="$MODE" \
    -v target_chr="$TARGET_CHR" \
    -v region_start="$REGION_START" \
    -v region_end="$REGION_END" \
    -v output_file="$OUTPUT" \
'
BEGIN {
    OFS = "\t"
    n_written      = 0
    n_skipped_na   = 0
    n_skipped_test = 0
}

# Parse header
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "CHR")                  col_chr    = i
        if ($i == "SNP")                  col_snp    = i
        if ($i == "BP")                   col_bp     = i
        if ($i == "A1")                   col_a1     = i
        if ($i == "TEST")                 col_test   = i
        if ($i == "NMISS")                col_nmiss  = i
        if ($i == "BETA" || $i == "OR")   col_effect = i
        if ($i == "SE")                   col_se     = i
        if ($i == "P")                    col_p      = i
    }
    print "SNP", "CHR", "BP", "A1", "BETA", "VARBETA", "P", "N" > output_file
    next
}

# Keep only ADD (additive) test rows; skip covariate rows
$col_test != "ADD" { n_skipped_test++; next }

# Filter to target chromosome and position window
$col_chr != target_chr || $col_bp < region_start || $col_bp > region_end { next }

# Skip rows with missing / sentinel values
$col_effect == "NA" || $col_effect == "." ||
$col_se     == "NA" || $col_se     == "." ||
$col_p      == "NA" || $col_p      == "." {
    n_skipped_na++
    next
}

{
    snp_val    = $col_snp
    chr_val    = $col_chr
    bp_val     = $col_bp
    a1_val     = $col_a1
    effect_val = $col_effect
    se_val     = $col_se
    p_val      = $col_p

    if (mode == "linear") {
        beta    = effect_val
        varbeta = se_val ^ 2
    } else {
        # Logistic: OR -> log(OR) = BETA (log-odds scale)
        # PLINK SE for logistic is already SE of log(OR)
        if (effect_val <= 0) { n_skipped_na++; next }
        beta    = log(effect_val)
        varbeta = se_val ^ 2
    }

    nmiss_val = (col_nmiss > 0 && $col_nmiss != "" && $col_nmiss != "NA") \
                ? ($col_nmiss + 0) : 0

    print snp_val, chr_val, bp_val, a1_val, beta, varbeta, p_val, nmiss_val > output_file
    n_written++
}

END {
    print "        SNPs written            : " n_written      > "/dev/stderr"
    print "        Skipped (NA/invalid)    : " n_skipped_na   > "/dev/stderr"
    print "        Skipped (non-ADD rows)  : " n_skipped_test > "/dev/stderr"

    if (n_written == 0) {
        print "[WARNING] No SNPs written. Check chromosome / region." \
            > "/dev/stderr"
    }
}
' "$ASSOC_FILE"

echo ""
echo "============================================================"
echo " Done!"
echo " Output : $OUTPUT"
echo "============================================================"
echo ""
echo " Columns : SNP  CHR  BP  A1  BETA  VARBETA  P  N"
echo " A1      : effect allele (PLINK minor/ALT allele after QC)"
echo " N       : NMISS (non-missing observations per SNP); median used by coloc"
echo " Note    : A2 is derived from eQTL ref/alt in prepare_coloc_datasets.R"
