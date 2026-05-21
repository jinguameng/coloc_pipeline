#!/usr/bin/env Rscript
# =============================================================================
# generate_summary_report.R
#
# Cross-locus PDF report. Reads the global all_coloc_summary.tsv produced
# by the Snakemake final aggregation step and lays out:
#
#   Page 1   Run parameters
#   Page 2   Cross-locus summary table (top pairs by best_PP4 across all loci)
#   Page 3   Per-locus best PP4 (one row per locus with its top pair)
#   Page 4   Histogram of best_PP4 across all (locus x dataset x gene) pairs
#
# Usage:
#   Rscript generate_summary_report.R \
#       --summary       <all_coloc_summary.tsv> \
#       --phenotype     SPAREAD \
#       --gwas_file     <path> \
#       --method        ABF|SUSIE|ALL \
#       --ancestry      EUR \
#       --study         GTEx \
#       --out           <summary_report.pdf>
# =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(dplyr); library(readr)
    library(ggplot2);  library(gridExtra); library(grid)
})

option_list <- list(
    make_option("--summary",   type="character"),
    make_option("--phenotype", type="character", default=""),
    make_option("--gwas_file", type="character", default=""),
    make_option("--method",    type="character", default=""),
    make_option("--ancestry",  type="character", default=""),
    make_option("--study",     type="character", default=""),
    make_option("--out",       type="character")
)
opts <- parse_args(OptionParser(option_list=option_list))

df <- readr::read_tsv(opts$summary, show_col_types=FALSE)

# --- NEW GRACEFUL EXIT FOR EMPTY RUNS ---
if (nrow(df) == 0) {
    cat("[WARN] Empty summary table. Generating placeholder report.\n")
    pdf(opts$out, width=8.5, height=11)
    grid::grid.newpage()
    grid::grid.text("No significant colocalization pairs found in this run.", 
                    gp=grid::gpar(fontsize=16, col="red"))
    dev.off()
    quit(save="no", status=0)
}
# ----------------------------------------

pdf(opts$out, width=11, height=8.5, onefile=TRUE)

# ─── Page 1: Run parameters ──────────────────────────────────────────────────
params <- rbind(
    c("Run date",         format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    c("Phenotype",        opts$phenotype),
    c("GWAS file",        opts$gwas_file),
    c("Method",           opts$method),
    c("Ancestry (plots)", opts$ancestry),
    c("Study",            opts$study),
    c("Loci tested",      sprintf("%d", length(unique(df$locus)))),
    c("Total pairs",      sprintf("%d", nrow(df))),
    c("Pairs PP4 >= 0.8", sprintf("%d", sum(df$best_PP4 >= 0.8, na.rm=TRUE))),
    c("Pairs PP4 >= 0.9", sprintf("%d", sum(df$best_PP4 >= 0.9, na.rm=TRUE)))
)
colnames(params) <- c("Parameter","Value")
title <- textGrob(sprintf("coloc_pipeline — Summary across %d loci",
                          length(unique(df$locus))),
                  gp=gpar(fontsize=13, fontface="bold"))
tbl <- tableGrob(params, rows=NULL,
                 theme=ttheme_minimal(base_size=10))
grid.arrange(title, tbl, ncol=1, heights=c(1, 14), newpage=FALSE)

# ─── Page 2: Top pairs across all loci ───────────────────────────────────────
top <- df %>% dplyr::arrange(dplyr::desc(best_PP4)) %>% head(40)
top_disp <- as.matrix(top %>% dplyr::transmute(
    locus, dataset_id, tissue, gene_id,
    PP4_ABF=sprintf("%.4f", PP4_ABF),
    PP4_SuSiE=sprintf("%.4f", PP4_SuSiE),
    best_PP4=sprintf("%.4f", best_PP4)
))
top_disp[is.na(top_disp)] <- ""
title2 <- textGrob(sprintf("Top %d pairs by best_PP4 (across all loci)",
                           nrow(top_disp)),
                   gp=gpar(fontsize=12, fontface="bold"))
tbl2 <- tableGrob(top_disp, rows=NULL,
                  theme=ttheme_minimal(base_size=8))
grid.arrange(title2, tbl2, ncol=1, heights=c(1, 14))

# ─── Page 3: Best pair per locus ─────────────────────────────────────────────
per_loc <- df %>% dplyr::group_by(locus) %>%
    dplyr::slice_max(best_PP4, n=1, with_ties=FALSE) %>% dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(best_PP4))
ploc_disp <- as.matrix(per_loc %>% dplyr::transmute(
    locus, dataset_id, tissue, gene_id,
    PP4_ABF=sprintf("%.4f", PP4_ABF),
    PP4_SuSiE=sprintf("%.4f", PP4_SuSiE),
    best_PP4=sprintf("%.4f", best_PP4),
    susie_detail
))
ploc_disp[is.na(ploc_disp)] <- ""
title3 <- textGrob("Best pair per locus",
                   gp=gpar(fontsize=12, fontface="bold"))
tbl3 <- tableGrob(ploc_disp, rows=NULL,
                  theme=ttheme_minimal(base_size=8))
grid.arrange(title3, tbl3, ncol=1, heights=c(1, 14))

# ─── Page 4: PP4 histogram ───────────────────────────────────────────────────
df_h <- df %>% dplyr::filter(is.finite(best_PP4))
if (nrow(df_h) > 0) {
    p <- ggplot(df_h, aes(best_PP4)) +
        geom_histogram(bins=40, fill="#2C7FB8", colour="white") +
        geom_vline(xintercept=c(0.5, 0.8, 0.9), linetype="dashed",
                   colour=c("grey50","darkorange","red")) +
        labs(x="best_PP4 (max of ABF and SuSiEx PP4)",
             y="Number of pairs",
             title=sprintf("Distribution of best_PP4 across %d pairs", nrow(df_h)),
             subtitle="Dashed lines: 0.5 / 0.8 / 0.9") +
        theme_bw(base_size=11)
    print(p)
}

invisible(dev.off())
cat(sprintf("[generate_summary_report.R] Wrote: %s\n", opts$out))
