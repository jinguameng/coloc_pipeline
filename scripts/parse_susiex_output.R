#!/usr/bin/env Rscript
# =============================================================================
# parse_susiex_output.R
#
# Convert SuSiEx output triple (.snp, .cs, .summary) into a single LBF table
# in the eQTL Catalogue `.lbf_variable.txt.gz` format. The downstream coloc
# step consumes this file with the same code path it uses for eQTL Catalogue
# SuSiE LBFs, so the conversion is the only SuSiEx-specific logic in the
# pipeline.
#
# ── What SuSiEx provides ──────────────────────────────────────────────────────
#
#   .summary  Header line `# chr19:44658684-45158684` + one row per credible
#             set (CS_ID, CS_LENGTH, CS_PURITY, MAX_PIP_SNP, BP, REF_ALLELE,
#             ALT_ALLELE, REF_FRQ, BETA, SE, -LOG10P, MAX_PIP,
#             POST-HOC_PROB_POP{i}). If SuSiEx did not converge the file
#             contains the literal string `FAIL`; if no CS survived filtering
#             it contains `NULL`. These two cases short-circuit this script.
#
#   .cs       One row per (CS, SNP). REF_ALLELE / ALT_ALLELE / BETA / SE /
#             -LOG10P / REF_FRQ are comma-separated per population. CS_PIP is
#             the SNP's PIP in the CS where it was assigned; OVRL_PIP is the
#             SNP's marginal PIP (i.e. the probability that the SNP is causal
#             in *any* signal in the locus).
#
#   .snp      One row per variant in the locus. Columns:
#               BP, SNP, then for each credible set k:
#                 PIP(CSk), LogBF(CSk, Pop1), LogBF(CSk, Pop2), ...,
#                 LogBF(CSk, PopP)
#             The number of CSes and populations is inferred from the header.
#             Missing-population entries are recorded as near-zero placeholders
#             (e.g. -9.25622e-09, -4.35127e-08).
#
# ── Cross-ancestry LBF combination ────────────────────────────────────────────
#
#   coloc.bf_bf() consumes a single L x N LBF matrix per side (L = signals,
#   N = variants). SuSiEx provides per-population LogBF instead. Under the
#   standard assumption of approximate independence between populations the
#   joint log Bayes factor for SNP j as the causal in signal k is the sum:
#
#       LBF_{k,j} = sum_{p=1..P} LogBF(CS_k, Pop_p, SNP_j)
#
#   Populations where the SNP is absent contribute the near-zero placeholders
#   above, so the summation behaves correctly without explicit NA handling.
#
#   No CS-filtering is needed: SuSiEx already only writes CSes that survived
#   its purity and marginal-p filters, so every CS column is a real signal.
#
# ── Allele resolution ─────────────────────────────────────────────────────────
#
#   The .snp file has no REF/ALT columns. Most SNP IDs are rsIDs, and only
#   the small subset assigned to credible sets carries allele info in .cs.
#   For coloc.bf_bf the LBF matrix must cover all variants (not just CS
#   members), so we need allele info for every variant. Resolution order:
#
#     1. PLINK2 .bim file (--bim, recommended). The 1KG reference is already
#        a pipeline dependency for LD plotting; using the SAME bim that the
#        SuSiEx reference panel was built from is ideal, but the 1KG bim is
#        a strong fallback that covers >99% of common variants.
#        plink2 .bim columns: chr id cm bp alt ref (ALT is the counted
#        allele; REF is the reference). Lookup is by SNP id first, then by
#        chromosome+position.
#     2. Parse from SNP id if it has chr:bp:ref:alt form
#        (e.g. "chr19:44909976:G:T") — SuSiEx writes some IDs this way.
#     3. .cs file. REF_ALLELE / ALT_ALLELE are comma-separated per
#        population; the first non-NA pair is used. Only CS members.
#     4. GWAS region file. Has A1 (effect allele) but not A2 — A2 stays NA
#        for these and the row gets dropped at the concordance check.
#     5. Otherwise ref/alt stay NA. The variant is kept in the output for
#        completeness but is dropped downstream when concordance is required.
#
# ── Output ────────────────────────────────────────────────────────────────────
#
#   Always writes:
#     <out_prefix>.lbf_variable.txt.gz   tabular LBF in eQTL Catalogue format
#     <out_prefix>.status                "PASS" | "FAIL" | "NULL"
#     <out_prefix>.cs_summary.tsv        per-CS metadata (purity, length,
#                                        max-PIP SNP, post-hoc per-pop probs)
#
#   On FAIL or NULL the LBF file is written with a header only and zero data
#   rows. Downstream code recognises this and falls back gracefully.
#
# ── Usage ─────────────────────────────────────────────────────────────────────
#
#   Rscript parse_susiex_output.R \
#       --susiex_dir   /path/to/susiex_output \
#       --susiex_name  SuSiEx.SPAREAD.apoe \
#       --gwas         /path/to/gwas_region.txt \
#       --bim          /path/to/1KG.bim \
#       --locus_name   apoe \
#       --out_prefix   /path/to/coloc_pipeline_out/loci/apoe/susiex
#
# Required R packages: optparse, dplyr, readr
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# ─── CLI arguments ────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--susiex_dir",  type="character",
                help="Directory containing SuSiEx outputs for this locus"),
    make_option("--susiex_name", type="character",
                help="Filename prefix without extension (e.g. SuSiEx.SPAREAD.apoe)"),
    make_option("--gwas",        type="character",
                help="GWAS region file (SNP CHR BP A1 BETA VARBETA P [N] [MAF])"),
    make_option("--bim",         type="character", default="",
                help="PLINK2 .bim file for allele lookup (recommended). Columns: chr id cm bp alt ref"),
    make_option("--locus_name",  type="character", default="GWAS",
                help="Locus name written to molecular_trait_id column [default: GWAS]"),
    make_option("--out_prefix",  type="character",
                help="Output file prefix (no extension)")
)
opts <- parse_args(OptionParser(option_list=option_list))

for (a in c("susiex_dir","susiex_name","gwas","out_prefix"))
    if (is.null(opts[[a]]) || nchar(trimws(opts[[a]])) == 0)
        stop(sprintf("[ERROR] --%s is required", a))

snp_path     <- file.path(opts$susiex_dir, paste0(opts$susiex_name, ".snp"))
cs_path      <- file.path(opts$susiex_dir, paste0(opts$susiex_name, ".cs"))
summary_path <- file.path(opts$susiex_dir, paste0(opts$susiex_name, ".summary"))

cat("============================================================\n")
cat(" parse_susiex_output.R\n")
cat(sprintf(" Locus     : %s\n", opts$locus_name))
cat(sprintf(" .snp      : %s\n", snp_path))
cat(sprintf(" .cs       : %s\n", cs_path))
cat(sprintf(" .summary  : %s\n", summary_path))
cat(sprintf(" GWAS file : %s\n", opts$gwas))
cat(sprintf(" Out prefix: %s\n", opts$out_prefix))
cat("============================================================\n\n")

dir.create(dirname(normalizePath(opts$out_prefix, mustWork=FALSE)),
           showWarnings=FALSE, recursive=TRUE)

status_path     <- paste0(opts$out_prefix, ".status")
lbf_path        <- paste0(opts$out_prefix, ".lbf_variable.txt.gz")
cs_summary_path <- paste0(opts$out_prefix, ".cs_summary.tsv")

# =============================================================================
# Helper: write_empty_lbf
# Writes a header-only LBF file plus a status marker. Used for FAIL / NULL.
# =============================================================================
write_empty_lbf <- function(status_str, n_cs_cols=1L) {
    writeLines(status_str, status_path)
    empty_lbf <- data.frame(
        molecular_trait_id = character(0),
        chromosome         = character(0),
        position           = integer(0),
        ref                = character(0),
        alt                = character(0),
        variant            = character(0),
        rsid               = character(0),
        cs_id              = character(0),
        pip                = numeric(0),
        stringsAsFactors   = FALSE
    )
    for (k in seq_len(n_cs_cols))
        empty_lbf[[paste0("lbf_variable", k)]] <- numeric(0)
    readr::write_tsv(empty_lbf, lbf_path)
    empty_cs <- data.frame(
        cs_id=character(0), n_snps=integer(0), purity=numeric(0),
        max_pip_snp=character(0), max_pip=numeric(0),
        stringsAsFactors=FALSE)
    readr::write_tsv(empty_cs, cs_summary_path)
}

# =============================================================================
# STEP 1 — Read .summary, detect status, get region chromosome
# =============================================================================
cat("[1/5] Reading .summary file...\n")
if (!file.exists(summary_path))
    stop(sprintf("[ERROR] .summary file not found: %s", summary_path))

summary_raw <- readLines(summary_path, warn=FALSE)
summary_raw <- summary_raw[nchar(trimws(summary_raw)) > 0]

if (length(summary_raw) == 0) {
    cat("      [WARN] .summary file is empty — treating as FAIL\n")
    write_empty_lbf("FAIL")
    cat(sprintf("\n      Wrote empty LBF + status='FAIL' to:\n        %s\n        %s\n",
                lbf_path, status_path))
    quit(status=0)
}

# Status line: a raw "FAIL" or "NULL" anywhere in the file is the signal
flat <- trimws(toupper(summary_raw))
if (any(flat == "FAIL")) {
    cat("      SuSiEx did not converge (FAIL) — writing empty LBF\n")
    write_empty_lbf("FAIL")
    cat(sprintf("\n      Wrote empty LBF + status='FAIL' to:\n        %s\n        %s\n",
                lbf_path, status_path))
    quit(status=0)
}
if (any(flat == "NULL")) {
    cat("      No credible set survived filtering (NULL) — writing empty LBF\n")
    write_empty_lbf("NULL")
    cat(sprintf("\n      Wrote empty LBF + status='NULL' to:\n        %s\n        %s\n",
                lbf_path, status_path))
    quit(status=0)
}

# Region header: first non-empty line should be `# chr19:44658684-45158684`
region_chr <- NA_character_
hdr_match  <- regmatches(summary_raw[1],
                         regexpr("chr[0-9XYM]+", summary_raw[1], ignore.case=TRUE))
if (length(hdr_match) > 0 && nchar(hdr_match) > 0) {
    region_chr <- sub("^chr", "", hdr_match, ignore.case=TRUE)
    cat(sprintf("      Region chromosome: %s\n", region_chr))
} else {
    cat("      [WARN] Could not parse chromosome from .summary header\n")
}

# Read summary body (skip lines starting with #)
body <- summary_raw[!grepl("^#", summary_raw)]
if (length(body) < 2) {
    cat("      [WARN] .summary has no credible sets after header — treating as NULL\n")
    write_empty_lbf("NULL")
    quit(status=0)
}
summary_df <- tryCatch(
    readr::read_tsv(I(paste(body, collapse="\n")), show_col_types=FALSE,
                    col_types=cols(.default="c")),
    error = function(e) NULL
)
if (is.null(summary_df) || nrow(summary_df) == 0) {
    cat("      [WARN] Could not parse .summary body — treating as NULL\n")
    write_empty_lbf("NULL")
    quit(status=0)
}
cat(sprintf("      %d credible set(s) in .summary\n", nrow(summary_df)))

# =============================================================================
# STEP 2 — Read .snp file, infer K credible sets and P populations
# =============================================================================
cat("\n[2/5] Reading .snp file...\n")
if (!file.exists(snp_path))
    stop(sprintf("[ERROR] .snp file not found: %s", snp_path))

snp_df <- readr::read_tsv(snp_path, show_col_types=FALSE,
                          col_types=cols(.default="c"))
cat(sprintf("      %d variants, %d columns\n", nrow(snp_df), ncol(snp_df)))

pip_cols <- grep("^PIP\\(CS",       colnames(snp_df), value=TRUE)
lbf_cols <- grep("^LogBF\\(CS",     colnames(snp_df), value=TRUE)
K <- length(pip_cols)
if (K == 0) {
    cat("      [WARN] No PIP(CS*) columns found in .snp — treating as NULL\n")
    write_empty_lbf("NULL")
    quit(status=0)
}
P <- length(lbf_cols) / K
if (P != round(P) || P < 1) {
    stop(sprintf("[ERROR] LogBF column count %d is not a multiple of K=%d CS",
                 length(lbf_cols), K))
}
P <- as.integer(P)
cat(sprintf("      Detected %d credible set(s) and %d population(s)\n", K, P))

# Re-derive CS ids from column names ("PIP(CS3)" -> 3) so we don't assume
# ordering. SuSiEx writes them sequentially but we are defensive.
cs_idx_in_file <- suppressWarnings(as.integer(
    sub(".*CS([0-9]+)\\).*", "\\1", pip_cols)))
cs_idx_in_file <- cs_idx_in_file[is.finite(cs_idx_in_file)]
if (length(cs_idx_in_file) != K)
    stop("[ERROR] Could not parse CS indices from .snp header")

# =============================================================================
# STEP 3 — Sum LogBF across populations per (SNP, CS)
# =============================================================================
cat("\n[3/5] Summing LogBF across populations...\n")

# Build a K x N matrix; rows = CSes, cols = variants.
# For each CS k, the population columns are the (k-1)*P + (1..P) entries of
# the LogBF columns *in the order they appear*. Be defensive and re-derive
# the column mapping from the LogBF header names.
lbf_meta <- data.frame(
    colname = lbf_cols,
    cs      = suppressWarnings(as.integer(sub(".*CS([0-9]+),.*", "\\1", lbf_cols))),
    pop     = suppressWarnings(as.integer(sub(".*Pop([0-9]+).*",  "\\1", lbf_cols))),
    stringsAsFactors = FALSE
)
if (anyNA(lbf_meta$cs) || anyNA(lbf_meta$pop))
    stop("[ERROR] Could not parse LogBF column names")

snp_id_vec <- snp_df$SNP
bp_vec     <- suppressWarnings(as.integer(snp_df$BP))

# K x N matrix in CS order (cs_idx_in_file)
lbf_mat <- matrix(0.0, nrow=K, ncol=nrow(snp_df))
for (i in seq_len(K)) {
    cs_k       <- cs_idx_in_file[i]
    cols_for_k <- lbf_meta$colname[lbf_meta$cs == cs_k]
    if (length(cols_for_k) == 0)
        stop(sprintf("[ERROR] No LogBF columns for CS%d", cs_k))
    # Sum across populations, coercing strings to numeric. Near-zero
    # placeholders contribute essentially nothing.
    mat_k <- vapply(cols_for_k, function(cn) {
        suppressWarnings(as.numeric(snp_df[[cn]]))
    }, numeric(nrow(snp_df)))
    if (is.null(dim(mat_k))) mat_k <- matrix(mat_k, ncol=1)
    s <- rowSums(mat_k, na.rm=TRUE)
    lbf_mat[i, ] <- s
}
rownames(lbf_mat) <- paste0("lbf_variable", seq_len(K))
colnames(lbf_mat) <- snp_id_vec
cat(sprintf("      LBF matrix : %d x %d (CSes x variants)\n", K, ncol(lbf_mat)))
cat(sprintf("      LBF range  : [%.3f, %.3f]\n",
            min(lbf_mat, na.rm=TRUE), max(lbf_mat, na.rm=TRUE)))

# Marginal PIP per SNP: sum PIP(CSi) across i (equivalent to OVRL_PIP for
# CS-member SNPs; safe approximation for non-CS SNPs).
pip_mat <- vapply(pip_cols, function(cn) {
    v <- suppressWarnings(as.numeric(snp_df[[cn]]))
    v[!is.finite(v)] <- 0
    v
}, numeric(nrow(snp_df)))
if (is.null(dim(pip_mat))) pip_mat <- matrix(pip_mat, ncol=1)
pip_marginal <- rowSums(pip_mat)
pip_marginal[pip_marginal > 1] <- 1   # numerical guard

# =============================================================================
# STEP 4 — Read .cs file, assign cs_id per SNP, gather allele info
# =============================================================================
cat("\n[4/5] Reading .cs file...\n")
cs_df <- if (file.exists(cs_path)) {
    tryCatch(
        readr::read_tsv(cs_path, show_col_types=FALSE,
                        col_types=cols(.default="c")),
        error = function(e) NULL
    )
} else NULL

cs_id_vec   <- rep(NA_character_, nrow(snp_df))
ref_from_cs <- rep(NA_character_, nrow(snp_df))
alt_from_cs <- rep(NA_character_, nrow(snp_df))

# Helper: take first non-NA element from a comma-separated allele string
first_non_na_allele <- function(x) {
    if (is.na(x) || nchar(trimws(x)) == 0) return(NA_character_)
    parts <- trimws(strsplit(x, ",", fixed=TRUE)[[1]])
    parts <- parts[parts != "" & toupper(parts) != "NA"]
    if (length(parts) == 0) NA_character_ else parts[1]
}

if (!is.null(cs_df) && nrow(cs_df) > 0 &&
    all(c("CS_ID","SNP") %in% colnames(cs_df))) {
    n_cs_rows <- nrow(cs_df)
    cat(sprintf("      %d (CS, SNP) rows\n", n_cs_rows))
    # First-CS-wins assignment, matching the convention used elsewhere.
    cs_df <- cs_df[order(suppressWarnings(as.integer(cs_df$CS_ID))), , drop=FALSE]
    m <- match(cs_df$SNP, snp_id_vec)
    keep <- !is.na(m) & is.na(cs_id_vec[m])
    cs_id_vec[m[keep]] <- paste0("L", cs_df$CS_ID[keep])
    if ("REF_ALLELE" %in% colnames(cs_df))
        ref_from_cs[m[keep]] <- vapply(cs_df$REF_ALLELE[keep],
                                       first_non_na_allele, character(1))
    if ("ALT_ALLELE" %in% colnames(cs_df))
        alt_from_cs[m[keep]] <- vapply(cs_df$ALT_ALLELE[keep],
                                       first_non_na_allele, character(1))
    n_assigned <- sum(!is.na(cs_id_vec))
    cat(sprintf("      Assigned cs_id to %d SNPs (unique)\n", n_assigned))
} else {
    cat("      [WARN] .cs file missing or unparseable — cs_id will be all NA\n")
}

# =============================================================================
# STEP 5 — Resolve alleles & chromosome per SNP, write output
# =============================================================================
cat("\n[5/5] Resolving alleles and chromosome, writing output...\n")

# (a) BIM lookup (primary). plink2 .bim columns: chr id cm bp alt ref
ref_from_bim <- rep(NA_character_, nrow(snp_df))
alt_from_bim <- rep(NA_character_, nrow(snp_df))
chr_from_bim <- rep(NA_character_, nrow(snp_df))
if (nchar(trimws(opts$bim)) > 0) {
    if (!file.exists(opts$bim)) {
        cat(sprintf("      [WARN] --bim file not found: %s — skipping BIM lookup\n",
                    opts$bim))
    } else {
        # Read only what we need; cheap even for 1KG-scale bim.
        # plink2 .bim has no header. Column types: chr(c) id(c) cm(d) bp(i) alt(c) ref(c)
        bim <- tryCatch(
            readr::read_tsv(opts$bim, col_names=c("chr","id","cm","bp","alt","ref"),
                            col_types="ccdicc", show_col_types=FALSE),
            error = function(e) NULL
        )
        # Fallback parse: column-oriented read with defaults that accept space or tab
        if (is.null(bim)) {
            bim <- tryCatch(
                read.table(opts$bim, header=FALSE,
                           col.names=c("chr","id","cm","bp","alt","ref"),
                           colClasses=c("character","character","numeric",
                                        "integer","character","character"),
                           comment.char="", stringsAsFactors=FALSE),
                error = function(e) NULL
            )
        }
        if (!is.null(bim) && nrow(bim) > 0) {
            cat(sprintf("      BIM loaded: %d variants\n", nrow(bim)))
            # Primary key: rsID
            m <- match(snp_id_vec, bim$id)
            hit <- !is.na(m)
            ref_from_bim[hit] <- bim$ref[m[hit]]
            alt_from_bim[hit] <- bim$alt[m[hit]]
            chr_from_bim[hit] <- bim$chr[m[hit]]
            # Fallback by chr:bp for remaining (within the locus's chromosome only)
            miss <- !hit
            if (any(miss) && !is.na(region_chr)) {
                bim_sub <- bim[bim$chr == region_chr, , drop=FALSE]
                key_bim <- paste0(bim_sub$chr, "_", bim_sub$bp)
                key_snp <- paste0(region_chr, "_", bp_vec[miss])
                m2 <- match(key_snp, key_bim)
                h2 <- !is.na(m2)
                ref_from_bim[which(miss)[h2]] <- bim_sub$ref[m2[h2]]
                alt_from_bim[which(miss)[h2]] <- bim_sub$alt[m2[h2]]
                chr_from_bim[which(miss)[h2]] <- bim_sub$chr[m2[h2]]
            }
            cat(sprintf("      BIM hits (id+pos): %d / %d\n",
                        sum(!is.na(ref_from_bim)), nrow(snp_df)))
        } else {
            cat("      [WARN] Could not read BIM file — skipping BIM lookup\n")
        }
    }
}

# (b) parse from SNP id if chr:bp:ref:alt format
ref_from_id <- rep(NA_character_, nrow(snp_df))
alt_from_id <- rep(NA_character_, nrow(snp_df))
chr_from_id <- rep(NA_character_, nrow(snp_df))
parts_lst   <- strsplit(snp_id_vec, ":", fixed=TRUE)
for (i in seq_along(parts_lst)) {
    p <- parts_lst[[i]]
    if (length(p) >= 4) {
        chr_from_id[i] <- sub("^chr", "", p[1], ignore.case=TRUE)
        ref_from_id[i] <- p[3]
        alt_from_id[i] <- p[4]
    }
}

# (c) GWAS region file lookup
gwas <- tryCatch(
    readr::read_tsv(opts$gwas, show_col_types=FALSE,
                    col_types=cols(.default="c")),
    error = function(e) NULL
)
ref_from_gwas <- rep(NA_character_, nrow(snp_df))
alt_from_gwas <- rep(NA_character_, nrow(snp_df))
chr_from_gwas <- rep(NA_character_, nrow(snp_df))
if (!is.null(gwas) && all(c("SNP","CHR","BP","A1") %in% colnames(gwas))) {
    # Primary key: rsID
    m <- match(snp_id_vec, gwas$SNP)
    hit <- !is.na(m)
    ref_from_gwas[hit] <- gwas$A1[m[hit]]   # A1 is GWAS effect allele; treat as REF placeholder
    chr_from_gwas[hit] <- gwas$CHR[m[hit]]
    # Fallback by BP — only one row per BP in the GWAS region file is assumed.
    miss <- !hit
    if (any(miss)) {
        m2 <- match(as.integer(bp_vec[miss]), suppressWarnings(as.integer(gwas$BP)))
        h2 <- !is.na(m2)
        ref_from_gwas[which(miss)[h2]] <- gwas$A1[m2[h2]]
        chr_from_gwas[which(miss)[h2]] <- gwas$CHR[m2[h2]]
    }
}

# Final allele resolution: BIM > id > cs > gwas; alt may stay NA after gwas
ref_final <- ref_from_bim
alt_final <- alt_from_bim
need <- is.na(ref_final) | is.na(alt_final)
ref_final[need & !is.na(ref_from_id)] <- ref_from_id[need & !is.na(ref_from_id)]
alt_final[need & !is.na(alt_from_id)] <- alt_from_id[need & !is.na(alt_from_id)]
need <- is.na(ref_final) | is.na(alt_final)
ref_final[need & !is.na(ref_from_cs)] <- ref_from_cs[need & !is.na(ref_from_cs)]
alt_final[need & !is.na(alt_from_cs)] <- alt_from_cs[need & !is.na(alt_from_cs)]
need <- is.na(ref_final)
ref_final[need & !is.na(ref_from_gwas)] <- ref_from_gwas[need & !is.na(ref_from_gwas)]
# alt remains NA if only GWAS hit — that variant gets dropped at concordance check

# Chromosome: prefer BIM, then SNP-id parse, then GWAS lookup, then region header
chr_final <- chr_from_bim
need <- is.na(chr_final) | nchar(chr_final) == 0
chr_final[need] <- chr_from_id[need]
need <- is.na(chr_final) | nchar(chr_final) == 0
chr_final[need] <- chr_from_gwas[need]
need <- is.na(chr_final) | nchar(chr_final) == 0
chr_final[need] <- region_chr

# Final variant id
variant_id <- ifelse(
    !is.na(chr_final) & !is.na(ref_final) & !is.na(alt_final),
    paste0("chr", chr_final, "_", bp_vec, "_", ref_final, "_", alt_final),
    NA_character_
)

# Quality counts for log
n_bim_alleles  <- sum(!is.na(ref_from_bim)  & !is.na(alt_from_bim))
n_id_alleles   <- sum(!is.na(ref_from_id)   & !is.na(alt_from_id))
n_cs_alleles   <- sum(!is.na(ref_from_cs)   & !is.na(alt_from_cs))
n_gwas_alleles <- sum(!is.na(ref_from_gwas))
n_full_alleles <- sum(!is.na(ref_final)     & !is.na(alt_final))
n_no_alleles   <- sum(is.na(ref_final)      | is.na(alt_final))
cat(sprintf("      Allele resolution:\n"))
cat(sprintf("        From BIM file       : %d\n", n_bim_alleles))
cat(sprintf("        Parsed from SNP id  : %d\n", n_id_alleles))
cat(sprintf("        From .cs file       : %d\n", n_cs_alleles))
cat(sprintf("        From GWAS file (A1) : %d\n", n_gwas_alleles))
cat(sprintf("        Full pair resolved  : %d / %d\n", n_full_alleles, nrow(snp_df)))
cat(sprintf("        Missing alleles     : %d  (will be dropped in coloc)\n",
            n_no_alleles))

# Build the eQTL Catalogue-format LBF data frame
out_df <- data.frame(
    molecular_trait_id = opts$locus_name,
    chromosome         = chr_final,
    position           = bp_vec,
    ref                = ref_final,
    alt                = alt_final,
    variant            = variant_id,
    rsid               = snp_id_vec,
    cs_id              = cs_id_vec,
    pip                = round(pip_marginal, 6),
    stringsAsFactors   = FALSE
)
for (i in seq_len(K))
    out_df[[paste0("lbf_variable", i)]] <- round(lbf_mat[i, ], 6)

readr::write_tsv(out_df, lbf_path)
writeLines("PASS", status_path)

# CS summary table for the report
cs_summary_out <- data.frame(
    cs_id          = paste0("L", summary_df$CS_ID),
    cs_length      = suppressWarnings(as.integer(summary_df$CS_LENGTH)),
    cs_purity      = suppressWarnings(as.numeric(summary_df$CS_PURITY)),
    max_pip_snp    = summary_df$MAX_PIP_SNP,
    max_pip_bp     = suppressWarnings(as.integer(summary_df$BP)),
    max_pip        = suppressWarnings(as.numeric(summary_df$MAX_PIP)),
    stringsAsFactors = FALSE
)
# Carry post-hoc population probabilities if present
posthoc_cols <- grep("^POST", colnames(summary_df), value=TRUE)
for (cn in posthoc_cols)
    cs_summary_out[[cn]] <- suppressWarnings(as.numeric(summary_df[[cn]]))
readr::write_tsv(cs_summary_out, cs_summary_path)

cat("\n============================================================\n")
cat(" Done!\n")
cat(sprintf(" Status     : PASS\n"))
cat(sprintf(" CSes       : %d\n", K))
cat(sprintf(" Variants   : %d\n", nrow(out_df)))
cat(sprintf(" LBF file   : %s\n", lbf_path))
cat(sprintf(" CS summary : %s\n", cs_summary_path))
cat(sprintf(" Status file: %s\n", status_path))
cat("============================================================\n")
