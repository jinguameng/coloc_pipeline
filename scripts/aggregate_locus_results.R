#!/usr/bin/env Rscript
# =============================================================================
# aggregate_locus_results.R
#
# Concatenate all per-pair summary_row.tsv fragments produced by
# run_coloc_one_pair.R for a single locus into one coloc_summary.tsv,
# sorted by best_PP4 descending. Snakemake invokes this once per locus
# after all run_coloc_one_pair rules for that locus complete.
#
# Usage:
#   Rscript aggregate_locus_results.R --in_dir <dir> --out <coloc_summary.tsv>
# =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(dplyr); library(readr)
})

option_list <- list(
    make_option("--in_dir", type="character", help="Directory of *.summary_row.tsv files"),
    make_option("--out",    type="character", help="Output coloc_summary.tsv path")
)
opts <- parse_args(OptionParser(option_list=option_list))

frags <- list.files(opts$in_dir, pattern="\\.summary_row\\.tsv$", full.names=TRUE)
if (length(frags) == 0) stop(sprintf("[ERROR] No *.summary_row.tsv in %s", opts$in_dir))

dfs <- lapply(frags, function(f) readr::read_tsv(f, show_col_types=FALSE))
df  <- dplyr::bind_rows(dfs) %>%
    dplyr::arrange(dplyr::desc(best_PP4))

dir.create(dirname(normalizePath(opts$out, mustWork=FALSE)),
           showWarnings=FALSE, recursive=TRUE)
readr::write_tsv(df, opts$out)

cat(sprintf("[aggregate_locus_results.R] %d rows → %s\n", nrow(df), opts$out))
cat(sprintf("    best_PP4 ≥ 0.8 : %d\n", sum(df$best_PP4 >= 0.8, na.rm=TRUE)))
