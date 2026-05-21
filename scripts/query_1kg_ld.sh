#!/usr/bin/env bash
# =============================================================================
# query_1kg_ld.sh
#
# Computes a signed LD matrix (Pearson r, REF-allele oriented) for a genomic
# region using 1000 Genomes Project GRCh38 data, subsetting to a specific
# ancestry superpopulation.
#
# Prerequisite (run once):
#   bash /data/h_vmac/zhanm32/colocpipe/scripts/build_1kg_ancestry_keep.sh \
#       --race /data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38_race.txt \
#       --out  /data/h_vmac/zhanm32/colocpipe/data/1KG
#
# ── Fixed internal paths (no user input required) ────────────────────────────
#   1KG bfile   : /data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38
#   Keep files  : /data/h_vmac/zhanm32/colocpipe/data/1KG/<POP>.keep
#   save_ld_rds : /data/h_vmac/zhanm32/colocpipe/scripts/save_ld_rds.R
#
# ── Usage ─────────────────────────────────────────────────────────────────────
#
#   bash query_1kg_ld.sh \
#       --ancestry  EUR \
#       --chr       19 \
#       --start     44783684 \
#       --end       45033684 \
#       --out       /path/to/output/EUR_chr19_region \
#       [--maf      0.01]
#
# ── Arguments ─────────────────────────────────────────────────────────────────
#
#   --ancestry   Superpopulation: EUR | AMR | AFR | EAS | SAS   [required]
#   --chr        Chromosome number (1-22)                        [required]
#   --start      Region start position (bp, GRCh38)              [required]
#   --end        Region end position (bp, GRCh38)                [required]
#   --out        Output file prefix (no extension)               [required]
#   --maf        Minimum MAF filter (default: 0.01)
#   --help / -h  Show this help message
#
# ── Output files ──────────────────────────────────────────────────────────────
#
#   <out>.ld       Signed LD matrix (tab-delimited, no header)
#   <out>.snplist  SNP IDs in matrix row/column order (one per line)
#   <out>.RDS      Named R matrix (SNP x SNP) — use this in coloc
#
# ── Sign convention ───────────────────────────────────────────────────────────
#
#   LD values are signed Pearson r relative to the REF allele
#   (--r-unphased ref-based), consistent with GWAS betas coded as effect
#   of ALT vs REF, and with eQTL Catalogue allele orientation.
#
# ── Example ───────────────────────────────────────────────────────────────────
#
#   bash /data/h_vmac/zhanm32/colocpipe/scripts/query_1kg_ld.sh \
#       --ancestry EUR \
#       --chr      19 \
#       --start    44783684 \
#       --end      45033684 \
#       --out      /data/h_vmac/zhanm32/colocpipe/tmp/1KG_EUR_chr19_APOE
#
#   In R:
#     LD <- readRDS("/data/h_vmac/zhanm32/colocpipe/tmp/1KG_EUR_chr19_APOE.RDS")
#     dim(LD)
#
# ── Note on sample sizes (for reference) ──────────────────────────────────────
#
#   EUR: 503    AMR: 347    AFR: 661    EAS: 504    SAS: 489
#
#   1KG has relatively small N per ancestry. For European analyses,
#   UKB-based LD (query_ukb_ld.sh) gives more stable estimates.
#   1KG is most valuable for AFR and AMR where UKB has no data.
#
# =============================================================================

set -euo pipefail

# ─── Fixed internal paths ─────────────────────────────────────────────────────
BFILE="/data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KEEP_DIR="$(dirname "$SCRIPT_DIR")/data/1KG"
SAVE_RDS_SCRIPT="${SCRIPT_DIR}/save_ld_rds.R"


# ─── Help ─────────────────────────────────────────────────────────────────────
if [[ $# -eq 0 ]] || echo "$*" | grep -qw -- '--help\|-h'; then
    sed -n '/^# =/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
    exit 0
fi

# ─── Defaults ─────────────────────────────────────────────────────────────────
ANCESTRY=""
CHR=""
START=""
END=""
OUT=""
MAF="0.01"

# ─── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ancestry) ANCESTRY="$2"; shift 2 ;;
        --chr)      CHR="$2";      shift 2 ;;
        --start)    START="$2";    shift 2 ;;
        --end)      END="$2";      shift 2 ;;
        --out)      OUT="$2";      shift 2 ;;
        --maf)      MAF="$2";      shift 2 ;;
        --help|-h)
            sed -n '/^# =/,/^# =/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            echo "        Run with --help for usage."
            exit 1 ;;
    esac
done

# ─── Validate user arguments ──────────────────────────────────────────────────
VALID_POPS="EUR AMR AFR EAS SAS"
for arg_name in ANCESTRY CHR START END OUT; do
    if [[ -z "${!arg_name}" ]]; then
        echo "[ERROR] Missing required argument: --${arg_name,,}"
        echo "        Run with --help for usage."
        exit 1
    fi
done

if ! echo "$VALID_POPS" | grep -qw "$ANCESTRY"; then
    echo "[ERROR] --ancestry must be one of: $VALID_POPS"
    exit 1
fi

# ─── Validate fixed paths ─────────────────────────────────────────────────────
KEEP_FILE="${KEEP_DIR}/${ANCESTRY}.keep"

if [[ ! -f "${BFILE}.bed" ]]; then
    echo "[ERROR] 1KG bfile not found: ${BFILE}.bed"
    echo "        Expected location: $BFILE"
    exit 1
fi
if [[ ! -f "$KEEP_FILE" ]]; then
    echo "[ERROR] Keep file not found: $KEEP_FILE"
    echo "        Run build_1kg_ancestry_keep.sh first:"
    echo "          bash ${SCRIPT_DIR}/build_1kg_ancestry_keep.sh \\"
    echo "              --race /data/h_vmac/GWAS_QC/1000G_data/1000G_final_b38_race.txt \\"
    echo "              --out  ${KEEP_DIR}"
    exit 1
fi
if [[ ! -f "$SAVE_RDS_SCRIPT" ]]; then
    echo "[ERROR] save_ld_rds.R not found: $SAVE_RDS_SCRIPT"
    exit 1
fi

mkdir -p "$(dirname "$OUT")"

# ─── Load modules (ACCRE HPC — safely skipped outside HPC) ───────────────────
if command -v module &>/dev/null; then
    module load plink/2.00-20251019-avx2 2>/dev/null || true
fi

N_SAMPLES=$(wc -l < "$KEEP_FILE")

echo "============================================================"
echo " query_1kg_ld.sh"
echo "============================================================"
echo " Ancestry    : $ANCESTRY  (N = $N_SAMPLES samples)"
echo " Region      : chr${CHR}:${START}-${END}"
echo " MAF filter  : >= $MAF"
echo " Output      : $OUT"
echo "============================================================"
echo ""

# ─── Step 1: Compute signed LD matrix with PLINK2 ────────────────────────────
echo "[1/3] Computing signed LD matrix with PLINK2..."

# The 1KG bfile uses rsIDs natively. We keep them as-is.
# prepare_coloc_datasets.R matches GWAS to LD by rsID (GWAS SNP column),
# and matches GWAS to eQTL by positional key (chr_pos). The LBF matrix
# column names are remapped from rsID to chr_pos before coloc.bf_bf.
plink2 \
    --bfile               "$BFILE"   \
    --keep                "$KEEP_FILE" \
    --chr                 "$CHR"     \
    --from-bp             "$START"   \
    --to-bp               "$END"     \
    --maf                 "$MAF"     \
    --snps-only           just-acgt  \
    --max-alleles         2          \
    --r-unphased          ref-based square \
    --out                 "$OUT"     \
    --silent

# Rename PLINK2 output to .ld
if [[ -f "${OUT}.unphased.vcor1" ]]; then
    mv "${OUT}.unphased.vcor1" "${OUT}.ld"
elif [[ -f "${OUT}.vcor1" ]]; then
    mv "${OUT}.vcor1" "${OUT}.ld"
else
    echo "[ERROR] PLINK2 LD output not found — check ${OUT}.log"
    exit 1
fi

N_SNPS=$(wc -l < "${OUT}.ld")
echo "      LD matrix: ${N_SNPS} x ${N_SNPS} SNPs"

# ─── Step 2: Extract SNP list in matrix order ─────────────────────────────────
echo ""
echo "[2/3] Extracting SNP list..."

plink2 \
    --bfile               "$BFILE"   \
    --keep                "$KEEP_FILE" \
    --chr                 "$CHR"     \
    --from-bp             "$START"   \
    --to-bp               "$END"     \
    --maf                 "$MAF"     \
    --snps-only           just-acgt  \
    --max-alleles         2          \
    --write-snplist                  \
    --out                 "$OUT"     \
    --silent

echo "      SNP list: ${N_SNPS} variants"

# ─── Step 3: Save as named RDS matrix ─────────────────────────────────────────
echo ""
echo "[3/3] Converting to named RDS matrix..."

if command -v module &>/dev/null; then
    module load r/4.5.0 2>/dev/null || true
fi

Rscript "$SAVE_RDS_SCRIPT" \
    --ld   "${OUT}.ld" \
    --snps "${OUT}.snplist" \
    --out  "${OUT}.RDS"

echo ""
echo "============================================================"
echo " Done!"
echo "============================================================"
echo " LD matrix  : ${OUT}.ld"
echo " SNP list   : ${OUT}.snplist"
echo " RDS matrix : ${OUT}.RDS"
echo ""
echo " Use in prepare_coloc_datasets.R:"
echo "   --ld ${OUT}.RDS"
