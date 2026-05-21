#!/usr/bin/env Rscript
# =============================================================================
# query_eqtl_index.R
#
# Queries the pre-built eQTL index (eqtl_index_ge_p1e4.tsv.gz) for a genomic
# region of interest and reports candidate gene–tissue combinations to guide
# colocpipe inputs.txt configuration.
#
# Usage:
#   Rscript /data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/query_eqtl_index.R \
#       --index  /data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz \
#       --chr    19 \
#       --start  44658684 \
#       --end    45158684
#
# Or using a lead SNP + window (requires the pipeline's reference file):
#   Rscript /data/h_vmac/zhanm32/colocpipe/scripts/eqtl_index_build/query_eqtl_index.R \
#       --index   /data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz \
#       --snp     rs429358 \
#       --window  500 \
#       --ref     /data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt
#
# Output:
#   Console: summary table grouped by gene, ranked by min p-value.
#   Optionally: --out /path/to/results.tsv writes the full hit table.
#
# The summary table shows what to put in inputs.txt:
#   eQTL Catalogue tissue(s)  — tissue_label values to try
#   eQTL Catalogue gene ID(s) — gene_id values to test
#
# Required R packages: optparse, dplyr, readr, seqminer
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
    library(seqminer)
})

# ─── CLI arguments ────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--index",   type="character",
                default="/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz",
                help="Path to bgzipped + tabix-indexed eQTL index file"),

    # Region: provide either --chr/--start/--end OR --snp/--window/--ref
    make_option("--chr",     type="character", default=NULL,
                help="Chromosome (e.g. 19 — no 'chr' prefix)"),
    make_option("--start",   type="integer",   default=NULL,
                help="Region start bp (GRCh38)"),
    make_option("--end",     type="integer",   default=NULL,
                help="Region end bp (GRCh38)"),

    make_option("--snp",     type="character", default=NULL,
                help="Lead SNP rsID (alternative to --chr/--start/--end)"),
    make_option("--window",  type="integer",   default=500L,
                help="Window size in KB, centred on lead SNP [default: 500]"),
    make_option("--ref",     type="character",
                default="/data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt",
                help="SNP reference file (CHR SNP BP) — required if using --snp"),

    # Output options
    make_option("--out",     type="character", default=NULL,
                help="Optional: write full hit table to this TSV path"),
    make_option("--p_thresh",type="double",    default=1e-4,
                help="Max p-value to display [default: 1e-4 = index threshold]"),
    make_option("--min_hits",type="integer",   default=1L,
                help="Only show genes with at least this many index hits [default: 1]"),
    make_option("--top_n",   type="integer",   default=NULL,
                help="Show only top N gene-tissue combinations by min p-value")
)
opts <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opts$index))
    stop("Index file not found: ", opts$index)
if (!file.exists(paste0(opts$index, ".tbi")))
    stop("Tabix index not found: ", opts$index, ".tbi  -- run tabix on the index file first")

# ─── Resolve region ───────────────────────────────────────────────────────────
CHR   <- NULL
START <- NULL
END   <- NULL

if (!is.null(opts$snp)) {
    # Look up lead SNP position from reference file
    if (!file.exists(opts$ref))
        stop("Reference file not found: ", opts$ref,
             "\nProvide --ref or use --chr/--start/--end directly.")

    cat(sprintf("[lookup] Searching for %s in reference file...\n", opts$snp))

    # Reference file is space-delimited (not tab), with a header row.
    # read_table() handles one-or-more whitespace as delimiter.
    ref_df <- readr::read_table(opts$ref,
                                col_names    = c("CHR", "SNP", "BP"),
                                col_types    = "cci",
                                show_col_types = FALSE,
                                skip         = 1L)

    hit <- ref_df[ref_df$SNP == opts$snp, ]
    if (nrow(hit) == 0)
        stop(sprintf("SNP '%s' not found in reference file", opts$snp))

    CHR   <- as.character(hit$CHR[1])
    bp    <- hit$BP[1]
    half  <- (opts$window * 1000L) %/% 2L
    START <- max(1L, bp - half)
    END   <- bp + half
    cat(sprintf("[lookup] Found %s at CHR%s:%d\n", opts$snp, CHR, bp))
    cat(sprintf("[lookup] Window ±%d KB → %s:%d-%d\n", opts$window %/% 2L, CHR, START, END))

} else if (!is.null(opts$chr) && !is.null(opts$start) && !is.null(opts$end)) {
    CHR   <- as.character(opts$chr)
    CHR   <- sub("^chr", "", CHR, ignore.case = TRUE)  # strip chr prefix if present
    START <- opts$start
    END   <- opts$end
} else {
    stop("Provide either --chr/--start/--end OR --snp (+ optionally --window and --ref)")
}

# ─── Query the tabix-indexed index ───────────────────────────────────────────
region_str <- sprintf("%s:%d-%d", CHR, START, END)

cat(sprintf("\n[query] Region: %s\n", region_str))
cat(sprintf("[query] Index : %s\n", opts$index))
cat(sprintf("[query] p <= %.0e\n\n", opts$p_thresh))

# Read the column names from the '#' header line in the bgzipped file
read_index_header <- function(path) {
    con  <- gzfile(path, "rb")
    on.exit(try(close(con), silent = TRUE))
    line <- readLines(con, n = 1L, warn = FALSE)
    strsplit(sub("^#", "", trimws(line)), "\t")[[1]]
}

col_names <- tryCatch(
    read_index_header(opts$index),
    error = function(e) {
        message("[WARN] Could not read header: ", e$message,
                "\n       Using default column names.")
        c("chromosome","position","dataset_id","study_label",
          "tissue_label","gene_id","pvalue","beta","maf")
    }
)

# tabix query via seqminer
raw <- tryCatch(
    seqminer::tabix.read.table(
        tabixFile        = opts$index,
        tabixRange       = region_str,
        stringsAsFactors = FALSE
    ),
    error = function(e) {
        message("[ERROR] tabix query failed: ", e$message)
        NULL
    }
)

if (is.null(raw) || nrow(raw) == 0) {
    cat("  No eQTL associations found in this region at p < 1e-4.\n")
    cat("  Interpretation: either no gene has a cis-eQTL in this window in ge data,\n")
    cat("  or the signal is weaker than the index threshold (p > 1e-4).\n")
    cat("  Suggestion: widen your window or check a different quant method.\n")
    quit(status = 0, save = "no")
}

colnames(raw) <- col_names[seq_len(ncol(raw))]

# Coerce types
hits <- raw %>%
    dplyr::mutate(
        position = as.integer(position),
        pvalue   = suppressWarnings(as.double(pvalue)),
        beta     = suppressWarnings(as.double(beta)),
        maf      = suppressWarnings(as.double(maf))
    ) %>%
    dplyr::filter(is.finite(pvalue), pvalue <= opts$p_thresh)

if (nrow(hits) == 0) {
    cat(sprintf("  No hits at p <= %.0e in this region.\n", opts$p_thresh))
    quit(status = 0, save = "no")
}

# ─── Summarise: one row per gene × study × tissue combination ─────────────────
summary_df <- hits %>%
    dplyr::group_by(gene_id, study_label, tissue_label, dataset_id) %>%
    dplyr::summarise(
        n_hits    = dplyr::n(),
        min_pval  = min(pvalue,  na.rm = TRUE),
        lead_beta = beta[which.min(pvalue)],
        lead_maf  = maf[which.min(pvalue)],
        lead_pos  = position[which.min(pvalue)],
        .groups = "drop"
    ) %>%
    dplyr::filter(n_hits >= opts$min_hits) %>%
    dplyr::arrange(min_pval)

if (!is.null(opts$top_n) && nrow(summary_df) > opts$top_n) {
    summary_df <- head(summary_df, opts$top_n)
}

# ─── Print formatted output ───────────────────────────────────────────────────
cat(sprintf("Found %d gene × tissue combination(s) with p ≤ %.0e\n",
            nrow(summary_df), opts$p_thresh))
cat(sprintf("Region : chr%s:%d-%d\n\n", CHR, START, END))

# Format the table for console
display_df <- summary_df %>%
    dplyr::mutate(
        min_pval_fmt  = formatC(min_pval, format="e", digits=2),
        lead_beta_fmt = sprintf("%+.3f", lead_beta),
        lead_maf_fmt  = sprintf("%.3f",  lead_maf),
        lead_pos_fmt  = format(lead_pos, big.mark=",")
    ) %>%
    dplyr::select(
        Gene      = gene_id,
        Study     = study_label,
        Tissue    = tissue_label,
        N_hits    = n_hits,
        Min_p     = min_pval_fmt,
        Lead_beta = lead_beta_fmt,
        Lead_MAF  = lead_maf_fmt,
        Lead_pos  = lead_pos_fmt
    )

print(as.data.frame(display_df), row.names = FALSE)

# ─── Suggested inputs.txt configuration ──────────────────────────────────────
top_genes   <- unique(summary_df$gene_id)
top_tissues <- unique(summary_df$tissue_label)
top_studies <- unique(summary_df$study_label)

cat("\n────────────────────────────────────────────────────────────\n")
cat("Suggested inputs.txt configuration\n")
cat("────────────────────────────────────────────────────────────\n")

if (length(top_studies) == 1L) {
    cat(sprintf("eQTL Catalogue study:                    %s\n", top_studies))
} else {
    cat(sprintf("eQTL Catalogue study:                    (multiple: %s)\n",
                paste(top_studies, collapse=", ")))
    cat("NOTE: colocpipe accepts one study per run. Run separately per study.\n")
}

cat(sprintf("eQTL Catalogue tissue(s), comma-sep:     %s\n",
            paste(top_tissues, collapse=", ")))
cat(sprintf("eQTL Catalogue gene ID(s), comma-sep:    %s\n",
            paste(top_genes,   collapse=", ")))
cat("────────────────────────────────────────────────────────────\n")

# ─── Write full hit table if requested ───────────────────────────────────────
if (!is.null(opts$out)) {
    readr::write_tsv(hits, opts$out)
    cat(sprintf("\nFull hit table (%d rows) written to: %s\n", nrow(hits), opts$out))
}