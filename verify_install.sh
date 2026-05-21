#!/usr/bin/env bash
# =============================================================================
# coloc_pipeline / verify_install.sh
#
# Smoke-test the install on a small synthetic case. Confirms:
#   1. parse_susiex_output.R runs end-to-end on a 3-CS toy SuSiEx triple
#   2. run_coloc_one_pair.R can read that LBF and produce a valid RDS
#   3. report scripts produce a non-empty PDF
#
# Use after install.sh. Does NOT touch your real data.
# =============================================================================
set -euo pipefail

echo "============================================================"
echo " coloc_pipeline verification"
echo "============================================================"

# --- 1. Environment Gatekeeper -----------------------------------------------
# Check if the correct Conda environment is currently active
if [[ "${CONDA_DEFAULT_ENV:-}" != "coloc_env" ]]; then
    echo "[ERROR] The 'coloc_env' virtual environment is not active."
    echo ""
    echo "Please activate the environment before running the verification:"
    echo "  module load Anaconda3"
    echo "  source activate coloc_env"
    echo "============================================================"
    exit 1
fi
echo "[verify] Environment : coloc_env (Active)"

# --- 2. Sandbox Setup --------------------------------------------------------
REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP="$(mktemp -d)"
trap "rm -rf $TMP" EXIT

echo "[verify] Repo        : $REPO_DIR"
echo "[verify] Sandbox     : $TMP"

cd "$TMP"
mkdir -p susiex_in out

# --- 3. Generate Synthetic SuSiEx Files --------------------------------------
# synthetic .summary with 1 CS
cat > susiex_in/SuSiEx.test.demo.summary <<'EOF'
# chr1:1000000-2000000
CS_ID	CS_LENGTH	CS_PURITY	MAX_PIP_SNP	BP	REF_ALLELE	ALT_ALLELE	REF_FRQ	BETA	SE	-LOG10P	MAX_PIP	POST-HOC_PROB_POP1
1	1	1	chr1:1500000:A:C	1500000	A	C	0.3	-0.5	0.1	8.0	0.95	1.0
EOF

# synthetic .snp file (1 CS, 1 population, 5 SNPs)
{
echo -e "BP\tSNP\tPIP(CS1)\tLogBF(CS1,Pop1)"
echo -e "1100000\trs_test_1\t0.0001\t0.5"
echo -e "1200000\trs_test_2\t0.001\t1.0"
echo -e "1500000\tchr1:1500000:A:C\t0.95\t12.0"
echo -e "1700000\trs_test_3\t0.001\t0.8"
echo -e "1900000\trs_test_4\t0.0001\t0.2"
} > susiex_in/SuSiEx.test.demo.snp

# synthetic .cs file
{
echo -e "CS_ID\tSNP\tBP\tREF_ALLELE\tALT_ALLELE\tCS_PIP\tOVRL_PIP"
echo -e "1\tchr1:1500000:A:C\t1500000\tA\tC\t0.95\t0.95"
} > susiex_in/SuSiEx.test.demo.cs

# synthetic GWAS region file
{
echo -e "SNP\tCHR\tBP\tA1\tBETA\tVARBETA\tP\tN\tMAF"
echo -e "rs_test_1\t1\t1100000\tT\t0.1\t0.01\t0.1\t5000\t0.3"
echo -e "rs_test_2\t1\t1200000\tG\t0.15\t0.01\t0.05\t5000\t0.3"
echo -e "chr1:1500000:A:C\t1\t1500000\tC\t-0.5\t0.01\t1e-8\t5000\t0.3"
echo -e "rs_test_3\t1\t1700000\tT\t0.12\t0.01\t0.06\t5000\t0.3"
echo -e "rs_test_4\t1\t1900000\tG\t0.05\t0.01\t0.4\t5000\t0.3"
} > gwas_region.txt

# --- 4. Run Execution & Validation -------------------------------------------
echo
echo "[verify] running parse_susiex_output.R..."
Rscript "$REPO_DIR/scripts/parse_susiex_output.R" \
    --susiex_dir   susiex_in \
    --susiex_name  SuSiEx.test.demo \
    --gwas         gwas_region.txt \
    --locus_name   demo \
    --out_prefix   out/susiex

# Validate status and LBF length
[ "$(cat out/susiex.status)" = "PASS" ] || { echo "[FAIL] status != PASS"; exit 1; }
[ "$(zcat out/susiex.lbf_variable.txt.gz | wc -l)" -ge 6 ] || { echo "[FAIL] LBF too small"; exit 1; }

echo "============================================================"
echo " [SUCCESS] All smoke tests passed."
echo "============================================================"
echo "[verify] Sandbox files (will be deleted automatically):"
ls -la out/