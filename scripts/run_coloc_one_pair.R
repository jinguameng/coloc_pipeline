#!/usr/bin/env Rscript
# =============================================================================
# run_coloc_one_pair.R
#
# Runs colocalization for ONE (locus, dataset, gene) triple. This is the
# Snakemake-friendly per-pair version of the legacy prepare_coloc_datasets.R
# (which looped through a manifest). Each Snakemake rule instance invokes
# this script once and produces one .coloc.RDS plus a single-row summary
# fragment.
#
# ── Inputs ────────────────────────────────────────────────────────────────────
#
#   --gwas              GWAS region file (SNP CHR BP A1 BETA VARBETA P [N] [MAF])
#                       produced by extract_coloc_region_PLINK.sh / _GWAMA.sh.
#   --eqtl_nominal      eQTL Catalogue nominal stats for the target gene
#                       (one gene, fetched by fetch_eqtl_region.R).
#   --eqtl_lbf          eQTL Catalogue LBF file for the target gene
#                       (.lbf.tsv, written by fetch_eqtl_region.R when
#                       method != ABF). Optional; if empty/missing the SuSiE
#                       arm is skipped for this pair.
#   --eqtl_cs           eQTL Catalogue credible-sets file for the target gene
#                       (.cs.tsv). Optional; if empty/missing the SuSiE arm
#                       is skipped.
#   --susiex_lbf        SuSiEx-derived GWAS LBF file in eQTL Catalogue format,
#                       produced by parse_susiex_output.R. Optional; if the
#                       paired .status file is FAIL/NULL or this file has 0
#                       data rows the SuSiE arm is skipped.
#   --susiex_status     SuSiEx status file (PASS / FAIL / NULL). Optional —
#                       inferred from --susiex_lbf if missing.
#
# ── Coloc params ──────────────────────────────────────────────────────────────
#
#   --method            ABF | SUSIE | ALL
#   --gwas_n            GWAS sample size N — if not provided, derived from
#                       the GWAS file's N column (median across the region).
#   --trait_type        quant | cc
#   --sdy               SD of phenotype (quant only; estimated from MAF if blank)
#   --s                 case proportion N_cases / N_total (cc only)
#
# ── Identity ──────────────────────────────────────────────────────────────────
#
#   --locus             Locus name (e.g. "apoe")
#   --dataset_id        eQTL Catalogue dataset id (e.g. "QTD000366")
#   --gene_id           Ensembl gene id (e.g. "ENSG00000130204")
#   --study_label       Study label (e.g. "GTEx") — for output
#   --tissue            Tissue label — for output
#   --cell_type         Sample group / condition — for output
#   --quant_method      ge | exon | tx | txrev — for output
#
# ── Output ────────────────────────────────────────────────────────────────────
#
#   --out               Output prefix; writes
#       <out>.coloc.RDS                full coloc result + aligned data + meta
#       <out>.summary_row.tsv          single-row summary fragment (one TSV
#                                      row aggregated later by aggregate_locus_results.R)
#
# ── Fallback logic ────────────────────────────────────────────────────────────
#
#   The user-chosen --method is respected where possible. If --method=SUSIE
#   but SuSiE inputs are unavailable (FAIL/NULL SuSiEx, missing eQTL CSes,
#   too few common variants, allele-concordance failure, etc.), ABF is run
#   automatically as a safety net so the pair still produces a result.
#   `susie_method` in the RDS / summary records the outcome.
#
# Required R packages: optparse, dplyr, readr, coloc, Matrix
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
    library(coloc)
    library(Matrix)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# ─── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--gwas",          type="character", help="GWAS region file"),
    make_option("--eqtl_nominal",  type="character", help="eQTL nominal stats file"),
    make_option("--eqtl_lbf",      type="character", default="", help="eQTL LBF file (optional)"),
    make_option("--eqtl_cs",       type="character", default="", help="eQTL credible-sets file (optional)"),
    make_option("--susiex_lbf",    type="character", default="", help="SuSiEx-parsed GWAS LBF file"),
    make_option("--susiex_status", type="character", default="", help="SuSiEx status file"),

    make_option("--method",        type="character", default="ABF",   help="ABF | SUSIE | ALL"),
    make_option("--gwas_n",        type="integer",   default=NA_integer_,
                help="GWAS sample size N (override; if absent, derived from N column)"),
    make_option("--trait_type",    type="character", default="quant", help="quant | cc"),
    make_option("--sdy",           type="character", default="",      help="SD of phenotype (quant only)"),
    make_option("--s",             type="character", default="",      help="Case proportion (cc only)"),

    make_option("--locus",         type="character", default="locus"),
    make_option("--dataset_id",    type="character", default=NA_character_),
    make_option("--gene_id",       type="character", default=NA_character_),
    make_option("--study_label",   type="character", default=NA_character_),
    make_option("--tissue",        type="character", default=NA_character_),
    make_option("--cell_type",     type="character", default=NA_character_),
    make_option("--quant_method",  type="character", default=NA_character_),

    make_option("--out",           type="character", help="Output file prefix")
)
opts <- parse_args(OptionParser(option_list=option_list))
opts$method     <- toupper(trimws(opts$method))
opts$trait_type <- tolower(trimws(opts$trait_type))

for (a in c("gwas","eqtl_nominal","out"))
    if (is.null(opts[[a]]) || nchar(trimws(opts[[a]])) == 0)
        stop(sprintf("[ERROR] --%s is required", a))

dir.create(dirname(normalizePath(opts$out, mustWork=FALSE)),
           showWarnings=FALSE, recursive=TRUE)

GWAS_SDY <- suppressWarnings(as.double(opts$sdy))
if (!is.finite(GWAS_SDY) || GWAS_SDY <= 0) GWAS_SDY <- NA_real_

GWAS_S <- suppressWarnings(as.double(opts$s))
if (!is.finite(GWAS_S) || GWAS_S <= 0 || GWAS_S >= 1) GWAS_S <- NA_real_

RUN_ABF   <- opts$method %in% c("ABF",   "ALL")
RUN_SUSIE <- opts$method %in% c("SUSIE", "ALL")

cat("============================================================\n")
cat(sprintf(" run_coloc_one_pair.R | %s.%s\n", opts$dataset_id, opts$gene_id))
cat(sprintf(" Locus       : %s\n", opts$locus))
cat(sprintf(" Method      : %s  (ABF=%s SUSIE=%s)\n",
            opts$method, ifelse(RUN_ABF,"y","n"), ifelse(RUN_SUSIE,"y","n")))
cat(sprintf(" Trait type  : %s\n", opts$trait_type))
cat("============================================================\n\n")

# =============================================================================
# Helpers (positional key, alignment, ABF, SuSiE arm)
# =============================================================================
pos_key <- function(chr, bp) paste0("chr", chr, "_", bp)

align_variants <- function(gwas_df, eqtl_df, min_snps=100L) {
    gwas_df$pos_key <- pos_key(gwas_df$CHR, gwas_df$BP)
    eqtl_df$pos_key <- pos_key(eqtl_df$chromosome, eqtl_df$position)

    common <- intersect(gwas_df$pos_key, eqtl_df$pos_key)
    if (length(common) < min_snps) return(NULL)

    g <- gwas_df[match(common, gwas_df$pos_key), ]
    e <- eqtl_df[match(common, eqtl_df$pos_key), ]

    dup <- duplicated(g$pos_key) | duplicated(e$pos_key)
    g <- g[!dup, ]; e <- e[!dup, ]

    ambig <- function(a, b) toupper(paste0(a, b)) %in% c("AT","TA","CG","GC")
    a1 <- toupper(g$A1); alt <- toupper(e$alt); ref <- toupper(e$ref)
    same_coding <- (a1 == alt) & !ambig(a1, ref)
    need_flip   <- (a1 == ref) & !ambig(a1, alt)
    same_coding[is.na(same_coding)] <- FALSE
    need_flip[is.na(need_flip)]     <- FALSE
    keep <- same_coding | need_flip

    g <- g[keep, ]; e <- e[keep, ]
    g$BETA[need_flip[keep]] <- -g$BETA[need_flip[keep]]
    list(gwas=g, eqtl=e)
}

# ── ABF arm (lifted from prepare_coloc_datasets.R, unchanged) ─────────────────
run_abf <- function(g, e, gwas_n, trait_type, sdy, s) {
    an_col <- if ("an" %in% colnames(e)) suppressWarnings(as.double(e$an)) else NULL
    n_eqtl <- if (!is.null(an_col)) round(median(an_col / 2, na.rm=TRUE)) else NA_real_
    if (!is.finite(n_eqtl) || n_eqtl < 10L) n_eqtl <- 200L

    if (trait_type == "cc") {
        d1 <- list(beta=g$BETA, varbeta=g$VARBETA, type="cc",
                   N=gwas_n, snp=g$pos_key)
        if (!is.na(s)) d1$s <- s
    } else {
        if (!is.na(sdy)) {
            d1 <- list(beta=g$BETA, varbeta=g$VARBETA, type="quant",
                       N=gwas_n, sdY=sdy, snp=g$pos_key)
        } else {
            if ("MAF" %in% colnames(g)) {
                gwas_maf <- suppressWarnings(as.double(g$MAF))
                maf_source <- "GWAS file MAF column"
            } else {
                gwas_maf <- suppressWarnings(as.double(e$maf))
                maf_source <- "eQTL MAF (GWAS MAF column not found)"
            }
            valid <- is.finite(gwas_maf) & gwas_maf > 0 & gwas_maf <= 0.5
            g <- g[valid, ]; e <- e[valid, ]; gwas_maf <- gwas_maf[valid]
            if (nrow(g) < 10L) {
                message("  [WARN] Too few SNPs with valid MAF — ABF skipped")
                return(list(result=NULL, warnings=character(0)))
            }
            message(sprintf("  [INFO] sdY estimated from MAF; source: %s (%d SNPs)",
                            maf_source, length(gwas_maf)))
            d1 <- list(beta=g$BETA, varbeta=g$VARBETA, type="quant",
                       N=gwas_n, MAF=gwas_maf, snp=g$pos_key)
        }
    }

    eqtl_maf <- suppressWarnings(as.double(e$maf))
    eqtl_maf[!is.finite(eqtl_maf) | eqtl_maf <= 0 | eqtl_maf > 0.5] <- NA_real_

    d2 <- list(beta=e$beta, varbeta=e$se^2, type="quant",
               N=n_eqtl, MAF=eqtl_maf, snp=e$pos_key)

    captured_warnings <- character(0)
    result <- withCallingHandlers(
        tryCatch(
            coloc::coloc.abf(d1, d2),
            error = function(err) {
                message("  [WARN] coloc.abf failed: ", err$message)
                NULL
            }
        ),
        warning = function(w) {
            captured_warnings <<- c(captured_warnings, conditionMessage(w))
            warning(w)
            invokeRestart("muffleWarning")
        }
    )
    list(result=result, warnings=captured_warnings)
}

# ── SuSiE arm (parsed-SuSiEx LBF × eQTL Catalogue LBF, via coloc.bf_bf) ───────
# The math is identical to the previous SuSiE arm: each side is an L x N LBF
# matrix; coloc.bf_bf does the heavy lifting. The only change is the *source*
# of the GWAS-side LBF — now SuSiEx (cross-ancestry, summed over populations
# in parse_susiex_output.R) instead of susie_rss against 1KG LD.

parse_cs_components <- function(cs_id_vec) {
    if (is.null(cs_id_vec)) return(integer(0))
    cs_id_vec <- cs_id_vec[!is.na(cs_id_vec) & nchar(trimws(cs_id_vec)) > 0]
    if (length(cs_id_vec) == 0) return(integer(0))
    m   <- regmatches(cs_id_vec, regexpr("[Ll](\\d+)$", cs_id_vec))
    idx <- suppressWarnings(as.integer(sub("^.*[Ll]", "", m)))
    idx <- idx[is.finite(idx) & idx > 0]
    sort(unique(idx))
}

read_eqtl_cs_file <- function(cs_path) {
    if (is.null(cs_path) || is.na(cs_path) || nchar(trimws(cs_path)) == 0) return(NULL)
    if (!file.exists(cs_path)) return(NULL)
    df <- tryCatch(
        readr::read_tsv(cs_path, col_types=cols(.default="c"), show_col_types=FALSE),
        error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
}

susiex_pass <- function(susiex_status_path, susiex_lbf_path) {
    # PASS only if status file says PASS (or absent and LBF file has data)
    if (nchar(trimws(susiex_status_path)) > 0 && file.exists(susiex_status_path)) {
        st <- trimws(readLines(susiex_status_path, warn=FALSE))
        if (length(st) == 0) return(FALSE)
        if (toupper(st[1]) != "PASS") return(FALSE)
    }
    if (nchar(trimws(susiex_lbf_path)) == 0 || !file.exists(susiex_lbf_path))
        return(FALSE)
    # Quick row count (header + at least 1 data line)
    con <- gzfile(susiex_lbf_path, "r"); on.exit(close(con))
    h <- readLines(con, n=2)
    length(h) >= 2
}

run_susie_coloc <- function(susiex_lbf_path, eqtl_lbf_path, eqtl_cs_path, gene_id) {
    if (!file.exists(susiex_lbf_path))
        return(list(result=NULL, method_used="ABF-fallback:no-gwas-lbf",
                    n_gwas_cs=NA_integer_, n_eqtl_cs=NA_integer_))
    if (!file.exists(eqtl_lbf_path))
        return(list(result=NULL, method_used="ABF-fallback:no-eqtl-lbf",
                    n_gwas_cs=NA_integer_, n_eqtl_cs=NA_integer_))

    gwas_lbf <- tryCatch(
        readr::read_tsv(susiex_lbf_path, col_types=cols(.default="c"), show_col_types=FALSE),
        error = function(e) NULL
    )
    eqtl_lbf <- tryCatch(
        readr::read_tsv(eqtl_lbf_path, col_types=cols(.default="c"), show_col_types=FALSE),
        error = function(e) NULL
    )
    if (is.null(gwas_lbf) || nrow(gwas_lbf) < 5L)
        return(list(result=NULL, method_used="ABF-fallback:gwas-lbf-empty",
                    n_gwas_cs=NA_integer_, n_eqtl_cs=NA_integer_))
    if (is.null(eqtl_lbf) || nrow(eqtl_lbf) < 5L)
        return(list(result=NULL, method_used="ABF-fallback:eqtl-lbf-empty",
                    n_gwas_cs=NA_integer_, n_eqtl_cs=NA_integer_))

    lbf_cols_g <- grep("^lbf_variable", colnames(gwas_lbf), value=TRUE)
    lbf_cols_e <- grep("^lbf_variable", colnames(eqtl_lbf), value=TRUE)
    if (length(lbf_cols_g) == 0 || length(lbf_cols_e) == 0)
        return(list(result=NULL, method_used="ABF-fallback:no-lbf-cols",
                    n_gwas_cs=NA_integer_, n_eqtl_cs=NA_integer_))

    # CS components on each side. The GWAS side's cs_id is written by
    # parse_susiex_output.R; every SuSiEx CS is "active" by construction
    # (SuSiEx already applied purity + p-value filters), so all GWAS LBF
    # columns are real. The eQTL side may have noise components.
    gwas_cs_idx <- if ("cs_id" %in% colnames(gwas_lbf))
        parse_cs_components(gwas_lbf$cs_id) else integer(0)
    # All GWAS LBF columns correspond to a SuSiEx-promoted CS — keep them all.
    # (If for some reason cs_id is empty, defensively keep all columns.)
    if (length(gwas_cs_idx) == 0) gwas_cs_idx <- seq_along(lbf_cols_g)

    eqtl_cs_idx <- integer(0)
    eqtl_cs_df  <- read_eqtl_cs_file(eqtl_cs_path)
    if (!is.null(eqtl_cs_df) && "cs_id" %in% colnames(eqtl_cs_df)) {
        if ("molecular_trait_id" %in% colnames(eqtl_cs_df))
            eqtl_cs_df <- dplyr::filter(eqtl_cs_df, molecular_trait_id == gene_id)
        eqtl_cs_idx <- parse_cs_components(eqtl_cs_df$cs_id)
    }

    n_gwas_cs <- length(gwas_cs_idx); n_eqtl_cs <- length(eqtl_cs_idx)
    cat(sprintf("    Credible sets: GWAS(SuSiEx) = %d | eQTL = %d (gene %s)\n",
                n_gwas_cs, n_eqtl_cs, gene_id))

    if (n_eqtl_cs == 0)
        return(list(result=NULL, method_used="ABF-only:no-credible-sets-eQTL",
                    n_gwas_cs=n_gwas_cs, n_eqtl_cs=0L))
    if (n_gwas_cs == 0)
        return(list(result=NULL, method_used="ABF-only:no-credible-sets-GWAS",
                    n_gwas_cs=0L, n_eqtl_cs=n_eqtl_cs))

    gwas_cs_idx <- gwas_cs_idx[gwas_cs_idx <= length(lbf_cols_g)]
    eqtl_cs_idx <- eqtl_cs_idx[eqtl_cs_idx <= length(lbf_cols_e)]
    if (length(gwas_cs_idx) == 0 || length(eqtl_cs_idx) == 0)
        return(list(result=NULL,
                    method_used="ABF-only:credible-set-index-out-of-range",
                    n_gwas_cs=n_gwas_cs, n_eqtl_cs=n_eqtl_cs))

    lbf_cols_g_active <- paste0("lbf_variable", gwas_cs_idx)
    lbf_cols_e_active <- paste0("lbf_variable", eqtl_cs_idx)

    # Filter eQTL LBF rows to target gene (one .lbf.tsv may contain many)
    if ("molecular_trait_id" %in% colnames(eqtl_lbf))
        eqtl_lbf <- dplyr::filter(eqtl_lbf, molecular_trait_id == gene_id)
    if (nrow(eqtl_lbf) == 0)
        return(list(result=NULL, method_used="ABF-only:no-credible-sets-eQTL",
                    n_gwas_cs=n_gwas_cs, n_eqtl_cs=0L))

    gwas_lbf$pos_key <- paste0("chr", gwas_lbf$chromosome, "_", gwas_lbf$position)
    eqtl_lbf$pos_key <- paste0("chr", eqtl_lbf$chromosome, "_", eqtl_lbf$position)

    common_pk <- intersect(gwas_lbf$pos_key, eqtl_lbf$pos_key)
    cat(sprintf("    LBF intersection: %d common pos_keys (GWAS: %d, eQTL: %d)\n",
                length(common_pk), nrow(gwas_lbf), nrow(eqtl_lbf)))
    if (length(common_pk) < 50L)
        return(list(result=NULL, method_used="ABF-fallback:too-few-common-snps",
                    n_gwas_cs=n_gwas_cs, n_eqtl_cs=n_eqtl_cs))

    g_meta <- gwas_lbf[match(common_pk, gwas_lbf$pos_key), ]
    e_meta <- eqtl_lbf[match(common_pk, eqtl_lbf$pos_key), ]

    # Allele concordance check (LBF is sign-invariant, so we just need the
    # same biallelic variant — ref/alt set must match, same or swapped)
    has_g_al <- all(c("ref","alt") %in% colnames(gwas_lbf)) &&
                !all(is.na(gwas_lbf$ref))
    has_e_al <- all(c("ref","alt") %in% colnames(eqtl_lbf)) &&
                !all(is.na(eqtl_lbf$ref))

    if (has_g_al && has_e_al) {
        gr <- toupper(g_meta$ref); ga <- toupper(g_meta$alt)
        er <- toupper(e_meta$ref); ea <- toupper(e_meta$alt)
        has_na <- is.na(gr) | is.na(ga) | is.na(er) | is.na(ea)
        ambig  <- (!has_na) & (paste0(gr, ga) %in% c("AT","TA","CG","GC") |
                               paste0(er, ea) %in% c("AT","TA","CG","GC"))
        same_orient <- (!has_na) & (gr == er) & (ga == ea)
        swapped     <- (!has_na) & (gr == ea) & (ga == er)
        concordant  <- (same_orient | swapped) & !ambig
        n_keep      <- sum(concordant, na.rm=TRUE)
        cat(sprintf("    Allele concordance: %d keep | %d ambig | %d discordant | %d NA\n",
                    n_keep, sum(ambig, na.rm=TRUE),
                    sum(!concordant & !ambig & !has_na, na.rm=TRUE),
                    sum(has_na, na.rm=TRUE)))
        if (n_keep < 50L)
            return(list(result=NULL, method_used="ABF-fallback:too-few-concordant-snps",
                        n_gwas_cs=n_gwas_cs, n_eqtl_cs=n_eqtl_cs))
        keep_pk <- common_pk[concordant]
    } else {
        cat("    Allele concordance: skipped (one side lacks ref/alt)\n")
        keep_pk <- common_pk
    }

    # Build LBF matrices for coloc.bf_bf (L x N: rows = signals, cols = variants)
    g_lbf_rows <- gwas_lbf[match(keep_pk, gwas_lbf$pos_key), lbf_cols_g_active, drop=FALSE]
    g_mat <- t(as.matrix(sapply(g_lbf_rows, as.numeric)))
    if (is.null(dim(g_mat))) g_mat <- matrix(g_mat, nrow=length(lbf_cols_g_active))
    storage.mode(g_mat) <- "double"
    colnames(g_mat) <- keep_pk; rownames(g_mat) <- lbf_cols_g_active

    e_lbf_rows <- eqtl_lbf[match(keep_pk, eqtl_lbf$pos_key), lbf_cols_e_active, drop=FALSE]
    e_mat <- t(as.matrix(sapply(e_lbf_rows, as.numeric)))
    if (is.null(dim(e_mat))) e_mat <- matrix(e_mat, nrow=length(lbf_cols_e_active))
    storage.mode(e_mat) <- "double"
    colnames(e_mat) <- keep_pk; rownames(e_mat) <- lbf_cols_e_active

    cat(sprintf("    coloc.bf_bf: GWAS %d x %d | eQTL %d x %d\n",
                nrow(g_mat), ncol(g_mat), nrow(e_mat), ncol(e_mat)))

    res <- tryCatch(
        coloc::coloc.bf_bf(g_mat, e_mat),
        error = function(err) {
            message("  [WARN] coloc.bf_bf failed: ", err$message)
            NULL
        }
    )
    if (is.null(res))
        return(list(result=NULL, method_used="ABF-fallback:bf-bf-error",
                    n_gwas_cs=n_gwas_cs, n_eqtl_cs=n_eqtl_cs))
    list(result=res, method_used="SuSiEx",
         n_gwas_cs=n_gwas_cs, n_eqtl_cs=n_eqtl_cs)
}

# =============================================================================
# STEP 1 — Load GWAS region (and derive median N if --gwas_n not given)
# =============================================================================
cat("[1/4] Loading GWAS region...\n")
gwas <- readr::read_tsv(opts$gwas, col_types=cols(.default="c"), show_col_types=FALSE) %>%
    dplyr::mutate(BP=as.integer(BP), BETA=as.double(BETA),
                  VARBETA=as.double(VARBETA), P=as.double(P)) %>%
    dplyr::filter(is.finite(BETA), is.finite(VARBETA), VARBETA > 0, is.finite(P))
if ("N" %in% colnames(gwas))   gwas$N   <- suppressWarnings(as.double(gwas$N))
if ("MAF" %in% colnames(gwas)) gwas$MAF <- suppressWarnings(as.double(gwas$MAF))
gwas$pos_key <- pos_key(gwas$CHR, gwas$BP)
N_GWAS_RAW <- nrow(gwas)
cat(sprintf("      %d GWAS SNPs loaded\n", N_GWAS_RAW))

# Derive effective N: user-provided > median of N column > error
EFFECTIVE_GWAS_N <- NA_integer_
if (!is.na(opts$gwas_n) && opts$gwas_n > 0L) {
    EFFECTIVE_GWAS_N <- opts$gwas_n
    cat(sprintf("      N source: user-provided (%d)\n", EFFECTIVE_GWAS_N))
} else if ("N" %in% colnames(gwas)) {
    vec <- gwas$N[is.finite(gwas$N) & gwas$N > 0]
    if (length(vec) > 0) {
        EFFECTIVE_GWAS_N <- as.integer(round(median(vec)))
        cat(sprintf("      N source: median of N column (%d; range %d–%d, %d SNPs)\n",
                    EFFECTIVE_GWAS_N, round(min(vec)), round(max(vec)), length(vec)))
    }
}
if (is.na(EFFECTIVE_GWAS_N) || EFFECTIVE_GWAS_N <= 0L)
    stop("[ERROR] Could not determine GWAS N (no --gwas_n and no usable N column)")

# =============================================================================
# STEP 2 — Load eQTL nominal and align
# =============================================================================
cat("\n[2/4] Loading eQTL nominal and aligning variants...\n")
eqtl <- readr::read_tsv(opts$eqtl_nominal, col_types=cols(.default="c"),
                        show_col_types=FALSE) %>%
    dplyr::mutate(
        position = as.integer(position),
        beta     = as.double(beta),
        se       = as.double(se),
        pvalue   = suppressWarnings(as.double(pvalue))
    ) %>%
    dplyr::filter(is.finite(beta), is.finite(se), se > 0)
N_EQTL_RAW <- nrow(eqtl)
cat(sprintf("      %d eQTL SNPs loaded\n", N_EQTL_RAW))

aligned <- align_variants(gwas, eqtl)
if (is.null(aligned))
    stop(sprintf("[ERROR] Too few common SNPs (< 100) — gwas=%d eqtl=%d",
                 N_GWAS_RAW, N_EQTL_RAW))
g <- aligned$gwas; e <- aligned$eqtl
cat(sprintf("      Aligned: %d common SNPs after harmonisation\n", nrow(g)))

# =============================================================================
# STEP 3 — Run coloc arms
# =============================================================================
cat("\n[3/4] Running coloc arms...\n")

abf_result   <- NULL; abf_warnings <- character(0)
pp_all_abf   <- rep(NA_real_, 5); names(pp_all_abf) <- paste0("PP.H", 0:4, ".abf")
pp4_abf      <- NA_real_

if (RUN_ABF) {
    cat("  → ABF arm\n")
    out <- run_abf(g, e, EFFECTIVE_GWAS_N, opts$trait_type, GWAS_SDY, GWAS_S)
    abf_result   <- out$result; abf_warnings <- out$warnings
    if (!is.null(abf_result)) {
        pp_all_abf <- as.numeric(abf_result$summary[paste0("PP.H", 0:4, ".abf")])
        names(pp_all_abf) <- paste0("PP.H", 0:4, ".abf")
        pp4_abf <- pp_all_abf["PP.H4.abf"]
        cat(sprintf("    ABF PP4: %.4f\n", pp4_abf))
    }
}

susie_result <- NULL; susie_method <- NA_character_
pp4_susie    <- NA_real_; n_gwas_cs <- NA_integer_; n_eqtl_cs <- NA_integer_

if (RUN_SUSIE) {
    cat("  → SuSiE arm (SuSiEx x eQTL Catalogue LBFs)\n")
    if (!susiex_pass(opts$susiex_status, opts$susiex_lbf)) {
        susie_method <- "ABF-fallback:susiex-fail-or-null"
        cat(sprintf("    SuSiE skipped: %s\n", susie_method))
    } else {
        sr <- run_susie_coloc(opts$susiex_lbf, opts$eqtl_lbf, opts$eqtl_cs,
                              opts$gene_id)
        susie_result <- sr$result
        susie_method <- sr$method_used
        n_gwas_cs    <- sr$n_gwas_cs %||% NA_integer_
        n_eqtl_cs    <- sr$n_eqtl_cs %||% NA_integer_
        if (!is.null(susie_result)) {
            pp4_susie <- max(susie_result$summary$PP.H4.abf, na.rm=TRUE)
            cat(sprintf("    SuSiE PP4: %.4f (best credible-set pair)\n", pp4_susie))
        } else {
            cat(sprintf("    SuSiE: %s\n", susie_method))
            # Safety net: if method=SUSIE only and SuSiE failed, run ABF.
            if (!RUN_ABF && is.null(abf_result)) {
                cat("    [INFO] SuSiE not available — running ABF as safety net\n")
                out <- run_abf(g, e, EFFECTIVE_GWAS_N, opts$trait_type, GWAS_SDY, GWAS_S)
                abf_result   <- out$result; abf_warnings <- out$warnings
                if (!is.null(abf_result)) {
                    pp_all_abf <- as.numeric(abf_result$summary[paste0("PP.H", 0:4, ".abf")])
                    names(pp_all_abf) <- paste0("PP.H", 0:4, ".abf")
                    pp4_abf <- pp_all_abf["PP.H4.abf"]
                    cat(sprintf("    Fallback ABF PP4: %.4f\n", pp4_abf))
                }
            }
        }
    }
}

if (is.null(abf_result) && is.null(susie_result))
    stop("[ERROR] All coloc arms failed — no result produced")

best_pp4 <- max(c(pp4_abf, pp4_susie), na.rm=TRUE)

# =============================================================================
# STEP 4 — Save RDS + single-row summary fragment
# =============================================================================
cat("\n[4/4] Saving outputs...\n")

result <- list(
    locus        = opts$locus,
    dataset_id   = opts$dataset_id,
    study_label  = opts$study_label,
    tissue       = opts$tissue,
    cell_type    = opts$cell_type,
    quant_method = opts$quant_method,
    gene_id      = opts$gene_id,
    method       = opts$method,
    susie_method = susie_method %||% NA_character_,
    n_gwas_cs    = n_gwas_cs,
    n_eqtl_cs    = n_eqtl_cs,
    trait_type   = opts$trait_type,
    gwas_n_used  = EFFECTIVE_GWAS_N,
    sdy_used     = GWAS_SDY,
    s_used       = GWAS_S,
    n_gwas_raw   = N_GWAS_RAW,
    n_eqtl_raw   = N_EQTL_RAW,
    n_snps       = nrow(g),
    pp_all_abf   = pp_all_abf,
    PP4_ABF      = pp4_abf,
    PP4_SuSiE    = pp4_susie,
    best_PP4     = best_pp4,
    abf_warnings = abf_warnings,
    abf_result   = abf_result,
    susie_result = susie_result,
    gwas_aligned = g,
    eqtl_aligned = e
)
saveRDS(result, paste0(opts$out, ".coloc.RDS"))


# --- NEW: Lookup EXACT tissue from manifest instead of using the broad user list ---
# Bulletproof manifest lookup based on the output directory path
manifest_path <- file.path(dirname(dirname(opts$out)), "eqtl_data", "manifest.tsv")
if (file.exists(manifest_path)) {
    manifest <- read.delim(manifest_path, sep="\t", stringsAsFactors=FALSE)
    exact_tissue <- manifest$tissue[manifest$dataset_id == opts$dataset_id][1]
    if (!is.na(exact_tissue) && exact_tissue != "") {
        opts$tissue <- exact_tissue
    }
}


summary_row <- data.frame(
    locus        = opts$locus,
    dataset_id   = opts$dataset_id %||% NA_character_,
    study_label  = opts$study_label %||% NA_character_,
    tissue       = opts$tissue %||% NA_character_,
    cell_type    = opts$cell_type %||% NA_character_,
    gene_id      = opts$gene_id %||% NA_character_,
    method       = opts$method,
    susie_detail = susie_method %||% NA_character_,
    n_gwas_cs    = n_gwas_cs,
    n_eqtl_cs    = n_eqtl_cs,
    n_snps       = nrow(g),
    PP4_ABF      = round(pp4_abf,   4),
    PP4_SuSiE    = round(pp4_susie, 4),
    best_PP4     = round(best_pp4,  4),
    result_file  = paste0(opts$out, ".coloc.RDS"),
    stringsAsFactors = FALSE
)
readr::write_tsv(summary_row, paste0(opts$out, ".summary_row.tsv"))

cat(sprintf("      Wrote: %s.coloc.RDS\n", opts$out))
cat(sprintf("      Wrote: %s.summary_row.tsv\n", opts$out))
cat("============================================================\n")
cat(sprintf(" Done | best_PP4 = %.4f | method = %s\n", best_pp4, opts$method))
cat("============================================================\n")
