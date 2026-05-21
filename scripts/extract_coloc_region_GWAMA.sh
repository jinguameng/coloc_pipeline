#!/usr/bin/env bash
# =============================================================================
# extract_coloc_region_GWAMA.sh
#
# Purpose: Extract a genomic region around a lead SNP from GWAMA meta-analysis
#          output (.out) and prepare for coloc:
#          SNP  CHR  BP  A1  BETA  VARBETA  P  N
#
# GWAMA output format (fixed column names, space- or tab-delimited):
#
#   Quantitative trait (--qt flag in GWAMA):
#     rs_number  reference_allele  other_allele  [eaf]
#     beta  [se]  beta_95L  beta_95U  z  p-value  ...  n_samples
#
#   Binary trait (default GWAMA):
#     rs_number  reference_allele  other_allele  [eaf]
#     OR  OR_95L  OR_95U  z  p-value  ...  n_samples
#
# Key design decisions:
#   1. CHR and BP are NEVER in GWAMA output (unless --map was used, which is
#      rare). This script always joins with a reference file by rsID.
#      Reference: /data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt
#      Format: CHR SNP BP  (space-delimited, one header row)
#      The reference is loaded into awk memory once — no temp files.
#
#   2. Column names are fixed by GWAMA documentation — no dynamic detection
#      loop needed. Only the delimiter and quant/binary mode are auto-detected.
#
#   3. SE derivation:
#        Quantitative — use 'se' column if present; otherwise
#          SE = (beta_95U - beta_95L) / 3.92
#        Binary — BETA = log(OR)
#          SE = (log(OR_95U) - log(OR_95L)) / 3.92
#      VARBETA = SE^2 in all cases.
#
#   4. n_samples from GWAMA is used as per-SNP N.
#      prepare_coloc_datasets.R takes median(N) across the region — no
#      user-supplied N is needed for GWAMA input.
#
# Usage:
#   bash extract_coloc_region_GWAMA.sh \
#       -f <gwama_out_file> \
#       -s <lead_snp_rsid> \
#       [-w <window_kb>]   \
#       [-o <output_file>] \
#       [-h]
#
# Arguments:
#   -f   Path to GWAMA .out file (.gz accepted)   [required]
#   -s   rsID of the lead SNP                     [required]
#   -w   Total window in KB, centred on lead SNP  [default: 250]
#   -o   Output file path                         [default: <snp>_coloc_<w>kb.txt]
#   -h   Show this help
#
# Output columns (tab-delimited):
#   SNP  CHR  BP  A1  BETA  VARBETA  P  N  MAF
#
#   MAF is derived from the GWAMA 'eaf' (effect allele frequency) column as:
#     MAF = eaf          if eaf <  0.5
#     MAF = 1 - eaf      if eaf >= 0.5
#   If 'eaf' is absent or NA, MAF is written as "" (empty).
#   prepare_coloc_datasets.R uses GWAS MAF for sdY estimation when sdY is not
#   provided, avoiding the less-precise fallback of using eQTL MAF.
# =============================================================================

set -euo pipefail

# ─── Fixed reference file path ───────────────────────────────────────────────
REF_FILE="/data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt"

# ─── Defaults ─────────────────────────────────────────────────────────────────
WINDOW_KB=250
OUTPUT=""
GWAMA_FILE=""
TARGET_SNP=""

# ─── Parse arguments ──────────────────────────────────────────────────────────
while getopts "f:s:w:o:h" opt; do
    case $opt in
        f) GWAMA_FILE="$OPTARG" ;;
        s) TARGET_SNP="$OPTARG" ;;
        w) WINDOW_KB="$OPTARG"  ;;
        o) OUTPUT="$OPTARG"     ;;
        h)
            sed -n '/^# Usage/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "[ERROR] Unknown option. Use -h for help." >&2; exit 1 ;;
    esac
done

# ─── Validate ─────────────────────────────────────────────────────────────────
if [[ -z "$GWAMA_FILE" || -z "$TARGET_SNP" ]]; then
    echo "[ERROR] -f (GWAMA file) and -s (target SNP) are required." >&2
    exit 1
fi
if [[ ! -f "$GWAMA_FILE" ]]; then
    echo "[ERROR] GWAMA file not found: $GWAMA_FILE" >&2
    exit 1
fi
if ! [[ "$WINDOW_KB" =~ ^[0-9]+$ ]] || [[ "$WINDOW_KB" -le 0 ]]; then
    echo "[ERROR] -w must be a positive integer (KB). Got: $WINDOW_KB" >&2
    exit 1
fi
if [[ ! -f "$REF_FILE" ]]; then
    echo "[ERROR] Reference file not found: $REF_FILE" >&2
    echo "        Expected format: CHR SNP BP (space-delimited, one header row)" >&2
    exit 1
fi

[[ -z "$OUTPUT" ]] && OUTPUT="${TARGET_SNP}_coloc_${WINDOW_KB}kb.txt"

HALF_WINDOW=$(( WINDOW_KB * 1000 / 2 ))

echo "============================================================"
echo " extract_coloc_region_GWAMA.sh"
echo "============================================================"
echo " GWAMA file  : $GWAMA_FILE"
echo " Target SNP  : $TARGET_SNP"
echo " Window      : ±${HALF_WINDOW} bp (${WINDOW_KB} KB total)"
echo " Reference   : $REF_FILE"
echo " Output      : $OUTPUT"
echo "============================================================"
echo ""

# ─── Handle .gz transparently ─────────────────────────────────────────────────
if file "$GWAMA_FILE" | grep -q 'gzip'; then
    CAT="zcat"
else
    CAT="cat"
fi

# ─── Step 1: Read header and detect format ────────────────────────────────────
echo "[1/3] Detecting GWAMA output format from header..."

# head -n1 causes SIGPIPE to cat/zcat under set -euo pipefail — disable it.
HEADER=$(set +o pipefail; $CAT "$GWAMA_FILE" | tr -d '\r' | head -n1)
echo "        Header: $HEADER"

# Detect delimiter: prefer tab, fall back to space.
# Use $'\t' (actual tab byte) and grep -F (fixed string) for portability.
TAB=$'\t'
if printf '%s' "$HEADER" | grep -qF "$TAB"; then
    DELIM="$TAB"
    echo "        Delimiter : TAB"
else
    DELIM=" "
    echo "        Delimiter : SPACE"
fi

# Detect quant vs binary and whether 'se' and 'n_samples' columns are present.
FORMAT=$(printf '%s' "$HEADER" | tr -d '\r' | awk -v delim="$DELIM" '
BEGIN { FS=delim }
{
    has_beta=0; has_or=0; has_se=0; has_n=0; has_eaf=0
    for (i=1; i<=NF; i++) {
        h = tolower($i)
        gsub(/[[:space:]]/, "", h)
        if (h == "beta")       has_beta = 1
        if (h == "or")         has_or   = 1
        if (h == "se")         has_se   = 1
        if (h == "n_samples")  has_n    = 1
        if (h == "eaf")        has_eaf  = 1
    }
    if (has_beta)    print "quant",  has_se, has_n, has_eaf
    else if (has_or) print "binary", 0,      has_n, has_eaf
    else             print "unknown", 0,     has_n, has_eaf
}')

TRAIT_MODE=$(echo "$FORMAT" | awk '{print $1}')
HAS_SE=$(    echo "$FORMAT" | awk '{print $2}')
HAS_N=$(     echo "$FORMAT" | awk '{print $3}')
HAS_EAF=$(   echo "$FORMAT" | awk '{print $4}')

if [[ "$TRAIT_MODE" == "unknown" ]]; then
    echo "[ERROR] Cannot detect trait type — header has neither 'beta' nor 'OR'." >&2
    echo "        Header: $HEADER" >&2
    exit 1
fi

echo "        Trait type : $TRAIT_MODE"
if [[ "$TRAIT_MODE" == "quant" ]]; then
    if [[ "$HAS_SE" == "1" ]]; then
        echo "        SE source  : 'se' column (direct)"
    else
        echo "        SE source  : derived from beta_95L / beta_95U  [SE=(U-L)/3.92]"
    fi
else
    echo "        BETA       : log(OR)"
    echo "        SE source  : derived from OR_95L / OR_95U  [SE=(log(U)-log(L))/3.92]"
fi
if [[ "$HAS_N" == "1" ]]; then
    echo "        N source   : n_samples column (median used by coloc)"
else
    echo "        N source   : n_samples not found — user-provided N used as fallback"
fi
if [[ "$HAS_EAF" == "1" ]]; then
    echo "        MAF source : eaf column  (MAF = eaf if eaf<0.5, else 1-eaf)"
else
    echo "        MAF source : eaf column NOT found — MAF column will be empty"
fi

# ─── Step 2: Locate lead SNP via reference file ───────────────────────────────
echo ""
echo "[2/3] Locating lead SNP '${TARGET_SNP}' in reference file..."

# GWAMA output never has CHR/BP. Look up the lead SNP in the reference.
# awk { exit } causes SIGPIPE — disable pipefail for this lookup.
SNP_INFO=$(set +o pipefail; awk -v snp="$TARGET_SNP" \
    'NR>1 && $2==snp { print $1 "\t" $3; exit }' "$REF_FILE")

if [[ -z "$SNP_INFO" ]]; then
    echo "[ERROR] SNP '${TARGET_SNP}' not found in reference file." >&2
    echo "        Reference : $REF_FILE" >&2
    exit 1
fi

TARGET_CHR=$(echo "$SNP_INFO" | cut -f1)
TARGET_BP=$(echo  "$SNP_INFO" | cut -f2)
REGION_START=$(( TARGET_BP - HALF_WINDOW ))
REGION_END=$(( TARGET_BP + HALF_WINDOW ))
[[ "$REGION_START" -lt 1 ]] && REGION_START=1

echo "        Found  : CHR${TARGET_CHR}:${TARGET_BP}"
echo "        Region : CHR${TARGET_CHR}:${REGION_START}-${REGION_END}"

# ─── Step 3: Extract region, join CHR/BP, compute BETA/VARBETA ────────────────
echo ""
echo "[3/3] Extracting region, joining with reference, computing BETA/VARBETA..."

$CAT "$GWAMA_FILE" | tr -d '\r' | awk \
    -v delim="$DELIM" \
    -v trait="$TRAIT_MODE" \
    -v has_se_col="$HAS_SE" \
    -v has_n_col="$HAS_N" \
    -v has_eaf_col="$HAS_EAF" \
    -v ref_file="$REF_FILE" \
    -v tchr="$TARGET_CHR" \
    -v rs="$REGION_START" \
    -v re="$REGION_END" \
    -v out="$OUTPUT" \
'
BEGIN {
    FS  = delim
    OFS = "\t"
    n_out = 0; n_skip = 0; n_noref = 0

    # Load reference file (CHR SNP BP, space-delimited) into memory.
    # Keyed by rsID so CHR/BP can be looked up for every GWAMA variant.
    while ((getline line < ref_file) > 0) {
        split(line, a, " ")
        if (a[1] == "CHR") continue
        ref_chr[a[2]] = a[1]
        ref_bp[a[2]]  = a[3] + 0
    }
    close(ref_file)

    print "SNP","CHR","BP","A1","BETA","VARBETA","P","N","MAF" > out
}

# ── Parse header: find column indices by known GWAMA column names ─────────────
# Column names are fixed by the GWAMA documentation:
#   rs_number  reference_allele  other_allele  [eaf]
#   beta / OR  [se]  beta_95L/beta_95U  OR_95L/OR_95U
#   z  p-value  -log10_p-value  q_statistic  q_p-value  i2
#   n_studies  n_samples  effects
NR == 1 {
    for (i = 1; i <= NF; i++) {
        h = tolower($i)
        gsub(/[[:space:]-]/, "", h)   # strip spaces and dashes (e.g. p-value → pvalue)
        if (h == "rsnumber"        || h == "rs_number")          col_snp   = i
        if (h == "referenceallele" || h == "reference_allele")   col_a1    = i
        if (h == "otherallele"     || h == "other_allele")       col_a2    = i
        if (h == "beta")                                          col_beta  = i
        if (h == "se")                                            col_se    = i
        if (h == "beta95l" || h == "beta_95l")                   col_b95l  = i
        if (h == "beta95u" || h == "beta_95u")                   col_b95u  = i
        if (h == "or")                                            col_or    = i
        if (h == "or95l"   || h == "or_95l")                     col_or95l = i
        if (h == "or95u"   || h == "or_95u")                     col_or95u = i
        if (h == "pvalue"  || h == "pvalues")                    col_p     = i
        if (h == "nsamples"|| h == "n_samples")                  col_n     = i
        if (h == "eaf")                                           col_eaf   = i
    }
    next
}

# ── Data rows ─────────────────────────────────────────────────────────────────
{
    snp_id = $col_snp

    # Look up CHR and BP from the reference
    chr_val = ref_chr[snp_id]
    bp_val  = ref_bp[snp_id] + 0
    if (chr_val == "" || bp_val == 0) { n_noref++; next }

    # Region filter
    if (chr_val != tchr || bp_val < rs || bp_val > re) next

    p_val = $col_p
    if (p_val == "" || p_val == "NA" || p_val == ".") { n_skip++; next }

    # ── BETA and SE by trait type ─────────────────────────────────────────────
    if (trait == "quant") {
        if ($col_beta == "" || $col_beta == "NA") { n_skip++; next }
        beta_val = $col_beta + 0

        if (has_se_col == "1" && $col_se != "" && $col_se != "NA") {
            # Direct SE column (present in some GWAMA versions)
            se_val = $col_se + 0
        } else {
            # Derive from 95% CI: SE = (beta_95U - beta_95L) / 3.92
            b95u = $col_b95u + 0
            b95l = $col_b95l + 0
            if ($col_b95u == "" || $col_b95l == "" || b95u == b95l) {
                n_skip++; next
            }
            se_val = (b95u - b95l) / 3.92
        }

    } else {
        # Binary: OR -> log(OR) = BETA; SE from log-scale CI
        or_val  = $col_or    + 0
        or95u   = $col_or95u + 0
        or95l   = $col_or95l + 0
        if (or_val <= 0 || or95u <= 0 || or95l <= 0 || or95u == or95l) {
            n_skip++; next
        }
        beta_val = log(or_val)
        se_val   = (log(or95u) - log(or95l)) / 3.92
    }

    if (se_val <= 0) { n_skip++; next }

    varbeta = se_val ^ 2

    n_val = (has_n_col == "1" && $col_n != "" && $col_n != "NA") \
            ? ($col_n + 0) : 0

    # MAF from EAF column (effect allele frequency)
    # MAF = eaf if eaf < 0.5, else 1 - eaf
    maf_val = ""
    if (has_eaf_col == "1" && col_eaf > 0 && $col_eaf != "" && $col_eaf != "NA") {
        eaf = $col_eaf + 0
        if (eaf > 0 && eaf < 1)
            maf_val = (eaf < 0.5) ? eaf : (1 - eaf)
    }

    print snp_id, chr_val, bp_val, $col_a1, beta_val, varbeta, p_val, n_val, maf_val > out
    n_out++
}

END {
    print "        SNPs written           : " n_out   > "/dev/stderr"
    print "        Skipped (NA/invalid)   : " n_skip  > "/dev/stderr"
    print "        Skipped (no ref match) : " n_noref > "/dev/stderr"
    if (n_out == 0)
        print "[WARNING] No SNPs written — check region and reference file" > "/dev/stderr"
}
'

# ─── Step 4: Calculate Median N ───────────────────────────────────────────────
# Extract column 8 (N), skip header, sort numerically, and use awk to find the median
MEDIAN_N=$(awk 'NR>1 {print $8}' "$OUTPUT" | sort -n | awk '
{ val[NR] = $1 }
END {
    if (NR == 0) {
        print "NA"
    } else if (NR % 2) {
        print val[(NR + 1) / 2]
    } else {
        print (val[NR / 2] + val[(NR / 2) + 1]) / 2.0
    }
}')

echo ""
echo "============================================================"
echo " Done! Output: $OUTPUT"
echo "============================================================"
echo " Columns : SNP  CHR  BP  A1  BETA  VARBETA  P  N  MAF"
if [[ "$TRAIT_MODE" == "quant" ]]; then
    echo " BETA    : 'beta' column"
    if [[ "$HAS_SE" == "1" ]]; then
        echo " VARBETA : se^2  (direct 'se' column)"
    else
        echo " VARBETA : ((beta_95U - beta_95L) / 3.92)^2"
    fi
else
    echo " BETA    : log(OR)"
    echo " VARBETA : ((log(OR_95U) - log(OR_95L)) / 3.92)^2"
fi
echo " N       : n_samples per SNP (median taken by prepare_coloc_datasets.R)"
echo " Median N: $MEDIAN_N"
echo " MAF     : derived from eaf column (MAF = eaf if eaf<0.5, else 1-eaf)"
echo " CHR/BP  : joined from reference file"