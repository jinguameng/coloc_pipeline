#!/usr/bin/env Rscript
# =============================================================================
# save_ld_rds.R
#
# Reads PLINK2 unphased LD output (.ld tab-delimited, no header) and the
# matching SNP list (.snplist, one ID per line), attaches row/column names,
# and saves the result as a named square matrix in RDS format.
#
# Called automatically by query_1kg_ld.sh — not intended for direct use.
#
# Usage (via query_1kg_ld.sh):
#   Rscript save_ld_rds.R --ld <file>.ld --snps <file>.snplist --out <file>.RDS
# =============================================================================

# Enforce shared ACCRE R library path
.libPaths(c("/data/h_vmac/waltes2/Coloc/Conda_Environment_New/rlib-4.5.0", .libPaths()))

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option("--ld",   type="character", help="Tab-delimited LD matrix (no header)"),
    make_option("--snps", type="character", help="SNP list file (one ID per line)"),
    make_option("--out",  type="character", help="Output RDS file path")
)
opts <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opts$ld), !is.null(opts$snps), !is.null(opts$out))

cat(sprintf("[save_ld_rds.R] Reading LD matrix: %s\n", opts$ld))
ld_mat <- as.matrix(read.table(opts$ld, header=FALSE, sep="\t",
                               check.names=FALSE, stringsAsFactors=FALSE))

snps <- trimws(readLines(opts$snps))
snps <- snps[nchar(snps) > 0]

if (nrow(ld_mat) != length(snps)) {
    stop(sprintf(
        "[save_ld_rds.R] Dimension mismatch: LD matrix is %d x %d but SNP list has %d entries.",
        nrow(ld_mat), ncol(ld_mat), length(snps)
    ))
}

# Report the ID format so prepare_coloc_datasets.R behaviour is predictable
example_id <- snps[1]
if (grepl("^rs[0-9]+$", example_id)) {
    cat(sprintf("[save_ld_rds.R] SNP ID format: rsID (e.g. %s)\n", example_id))
    cat("[save_ld_rds.R]   LD ↔ GWAS matching will use the GWAS SNP column (rsID)\n")
} else if (grepl("^chr[0-9XY]+_[0-9]+$", example_id)) {
    cat(sprintf("[save_ld_rds.R] SNP ID format: positional (e.g. %s)\n", example_id))
    cat("[save_ld_rds.R]   LD ↔ GWAS matching will use pos_key (chr_bp)\n")
} else {
    cat(sprintf("[save_ld_rds.R] SNP ID format: unrecognised (e.g. %s) — check manually\n",
                example_id))
}

rownames(ld_mat) <- snps
colnames(ld_mat) <- snps

# Ensure diagonal is exactly 1 (floating-point rounding from PLINK2)
diag(ld_mat) <- 1.0

# Symmetrise (PLINK2 upper-triangle output may have minor asymmetry)
ld_mat[lower.tri(ld_mat)] <- t(ld_mat)[lower.tri(ld_mat)]

cat(sprintf("[save_ld_rds.R] Saving %d x %d named LD matrix to: %s\n",
            nrow(ld_mat), ncol(ld_mat), opts$out))
saveRDS(ld_mat, opts$out)
cat("[save_ld_rds.R] Done.\n")
