#!/usr/bin/env Rscript
# =============================================================================
# generate_locus_report.R
#
# Per-locus PDF report. Reads every *.coloc.RDS file under --results_dir
# (each produced by run_coloc_one_pair.R for one dataset x gene), plus
# the parsed SuSiEx CS summary, and lays out:
#
#   Page 1   Run parameters (this locus's context)
#   Page 2   SuSiEx credible sets (from <prefix>.cs_summary.tsv)
#   Page 3   Cross-pair summary table (all dataset x gene tested at this locus)
#   Page 4+  One details page per (dataset, gene) — PP table + regional plot
#
# This is a minimal-viable port of the legacy generate_report.R. The legacy
# script's sensitivity plots, curl-progress log capture, and extensive ASCII
# normalisation are not reproduced here — those can be added back as needed
# without touching the parse/coloc core.
#
# Usage:
#   Rscript generate_locus_report.R \
#       --results_dir   <dir of *.coloc.RDS> \
#       --cs_summary    <path to susiex.cs_summary.tsv> \
#       --status        <path to susiex.status> \
#       --locus         apoe \
#       --lead_snp      rs429358 \
#       --window_kb     500 \
#       --phenotype     SPAREAD \
#       --gwas_file     <path> \
#       --gwas_type     PLINK|GWAMA \
#       --trait_type    quant|cc \
#       --method        ABF|SUSIE|ALL \
#       --ancestry      EUR \
#       --study         GTEx \
#       --out           <coloc_report.pdf>
# =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(dplyr); library(readr)
    library(ggplot2);  library(gridExtra); library(grid); library(scales)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

option_list <- list(
    make_option("--results_dir", type="character"),
    make_option("--cs_summary",  type="character", default=""),
    make_option("--status",      type="character", default=""),
    make_option("--ld",          type="character", default=""),
    make_option("--locus",       type="character"),
    make_option("--lead_snp",    type="character", default=""),
    make_option("--window_kb",   type="character", default=""),
    make_option("--phenotype",   type="character", default=""),
    make_option("--gwas_file",   type="character", default=""),
    make_option("--gwas_type",   type="character", default=""),
    make_option("--trait_type",  type="character", default=""),
    make_option("--method",      type="character", default=""),
    make_option("--ancestry",    type="character", default=""),
    make_option("--study",       type="character", default=""),
    make_option("--out",         type="character")
)
opts <- parse_args(OptionParser(option_list=option_list))


# ─── Load all RDS files for this locus ───────────────────────────────────────
rds_files <- list.files(opts$results_dir, pattern="\\.coloc\\.RDS$", full.names=TRUE)

# --- NEW GRACEFUL EXIT FOR EMPTY LOCI ---
if (length(rds_files) == 0) {
    cat("[WARN] No *.coloc.RDS found. Generating placeholder report.\n")
    suppressPackageStartupMessages(library(grid))
    pdf(opts$out, width=8.5, height=11)
    grid::grid.newpage()
    grid::grid.text(sprintf("Locus: %s\n\nNo significant colocalization pairs found.\n(eQTL pre-flight likely found no signals with p < 1e-4)", opts$locus), 
                    gp=grid::gpar(fontsize=14, col="red"))
    dev.off()
    quit(save="no", status=0)
}
# ----------------------------------------


# THIS IS THE CRUCIAL LINE THAT WAS LIKELY MISSING OR ALTERED:
results <- lapply(rds_files, readRDS)

# --- NEW FIX: Dynamically read exact tissue from manifest.tsv ---
manifest_path <- file.path(dirname(opts$results_dir), "eqtl_data", "manifest.tsv")
if (file.exists(manifest_path)) {
    manifest <- read.delim(manifest_path, sep="\t", stringsAsFactors=FALSE)
    results <- lapply(results, function(r) {
        # Find exact tissue for this specific dataset_id
        exact_tissue <- manifest$tissue[manifest$dataset_id == r$dataset_id][1]
        if (!is.na(exact_tissue) && exact_tissue != "") {
            r$tissue <- exact_tissue
        }
        return(r)
    })
}
# ----------------------------------------------------------------

cat(sprintf("[generate_locus_report.R] Loaded %d coloc results for locus %s\n",
            length(results), opts$locus))

# ─── Helpers ─────────────────────────────────────────────────────────────────
pp4_interpret <- function(p) {
    if (is.na(p))     return("?")
    if (p >= 0.9)     return("strong evidence for colocalization")
    if (p >= 0.8)     return("moderate evidence for colocalization")
    if (p >= 0.5)     return("weak evidence for colocalization")
    return("no evidence")
}

make_table_page <- function(rows, title_str, base_size=9) {
    title <- textGrob(title_str, gp=gpar(fontsize=13, fontface="bold"))
    tbl   <- tableGrob(
        rows, rows=NULL,
        theme=ttheme_minimal(
            base_size=base_size,
            core=list(fg_params=list(hjust=0, x=0.02)),
            colhead=list(fg_params=list(fontface="bold", hjust=0, x=0.02))))
    list(title=title, tbl=tbl)
}

regional_plot <- function(res, ld_mat, lead_snp) {
    g <- res$gwas_aligned; e <- res$eqtl_aligned
    if (is.null(g) || nrow(g) == 0) return(textGrob("No GWAS data"))
    
    # Calculate r2 against the lead SNP if LD is available
    r2_vec <- rep(NA_real_, nrow(g))
    if (!is.null(ld_mat) && lead_snp %in% rownames(ld_mat)) {
        ld_lead <- ld_mat[lead_snp, ]
        m <- match(g$SNP, names(ld_lead)) # Align LD to GWAS SNPs
        r2_vec <- (ld_lead[m])^2
    }

    df_g <- data.frame(BP=g$BP, neglogP=-log10(pmax(g$P, 1e-300)), src="GWAS", r2=r2_vec)
    
    if (!is.null(e) && nrow(e) > 0) {
        df_e <- data.frame(BP=e$position,
                           neglogP=-log10(pmax(e$pvalue, 1e-300)),
                           src=sprintf("eQTL (%s)", res$tissue %||% "?"),
                           r2=r2_vec) # Reuse aligned r2
        df <- rbind(df_g, df_e)
    } else df <- df_g

    ggplot(df, aes(BP, neglogP, color=r2)) +
        geom_point(alpha=0.8, size=1.5) +
        scale_color_gradient(low="blue", high="red", na.value="grey50", limits=c(0, 1)) +
        facet_wrap(~src, ncol=1, scales="free_y") +
        scale_x_continuous(labels=label_comma()) +
        labs(x=sprintf("Position on chr%s",
                       g$CHR[1] %||% e$chromosome[1] %||% "?"),
             y=expression(-log[10](P)),
             color=expression(r^2),
             title=sprintf("%s | %s", res$dataset_id, res$gene_id)) +
        theme_bw(base_size=10) +
        theme(strip.background=element_rect(fill="grey90"),
              legend.position="right")
}

# ─── Page 1: Run parameters ──────────────────────────────────────────────────
sus_status <- if (nchar(opts$status) > 0 && file.exists(opts$status))
    trimws(readLines(opts$status, warn=FALSE))[1] else "(not provided)"

params_rows <- rbind(
    c("Run date",      format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    c("Phenotype",     opts$phenotype),
    c("Locus",         opts$locus),
    c("Lead SNP",      opts$lead_snp),
    c("Window",        sprintf("%s KB", opts$window_kb)),
    c("", ""),
    c("--- GWAS ---", ""),
    c("GWAS file",     opts$gwas_file),
    c("GWAS format",   toupper(opts$gwas_type)),
    c("Trait type",    opts$trait_type),
    c("", ""),
    c("--- Coloc ---", ""),
    c("Method",        opts$method),
    c("SuSiEx status", sus_status),
    c("", ""),
    c("--- eQTL scope ---", ""),
    c("Study",         opts$study),
    c("Tissue(s)",     paste(unique(sapply(results, function(r) r$tissue %||% "?")), collapse=", ")),
    c("Genes tested",  paste(unique(sapply(results, function(r) r$gene_id %||% "?")), collapse=", ")),
    c("", ""),
    c("--- LD reference (regional plots only) ---", ""),
    c("Ancestry",      opts$ancestry),
    c("", ""),
    c("--- Results ---", ""),
    c("Pairs tested",       sprintf("%d", length(results))),
    c("Pairs PP4 >= 0.8",   sprintf("%d", sum(sapply(results, function(r) r$best_PP4 >= 0.8), na.rm=TRUE)))
)
colnames(params_rows) <- c("Parameter", "Value")

# --- NEW FIX: Wrap long text strings so they don't get cut off on Page 1 ---
params_rows[, "Value"] <- sapply(params_rows[, "Value"], function(x) {
    paste(strwrap(as.character(x), width=80), collapse="\n")
})
# ---------------------------------------------------------------------------

pdf(opts$out, width=11, height=8.5, onefile=TRUE)
p1 <- make_table_page(params_rows,
                      sprintf("coloc_pipeline — Run Parameters | Locus: %s", opts$locus))
grid.arrange(p1$title, p1$tbl, ncol=1, heights=c(1, 18), newpage=FALSE)

# ─── Page 2: SuSiEx CS metadata ──────────────────────────────────────────────
if (nchar(opts$cs_summary) > 0 && file.exists(opts$cs_summary)) {
    cs <- tryCatch(
        readr::read_tsv(opts$cs_summary, show_col_types=FALSE),
        error = function(e) NULL
    )
    if (!is.null(cs) && nrow(cs) > 0) {
        # Render as a wide table
        cs_disp <- as.matrix(cs)
        cs_disp[is.na(cs_disp)] <- ""
        # Round numeric-looking columns for display
        for (cn in c("cs_purity","max_pip")) {
            if (cn %in% colnames(cs))
                cs_disp[, cn] <- sprintf("%.4f", suppressWarnings(as.numeric(cs[[cn]])))
        }
        p2 <- make_table_page(cs_disp,
                              sprintf("SuSiEx credible sets at %s (status: %s)",
                                      opts$locus, sus_status),
                              base_size=8.5)
        grid.arrange(p2$title, p2$tbl, ncol=1, heights=c(1, 12))
    } else {
        grid.newpage()
        grid.text(sprintf("SuSiEx status: %s — no credible sets to display",
                          sus_status),
                  gp=gpar(fontsize=14, fontface="bold"))
    }
} else {
    grid.newpage()
    grid.text("SuSiEx CS summary not provided",
              gp=gpar(fontsize=14, fontface="bold"))
}

# ─── Page 3: Cross-pair summary table ────────────────────────────────────────
sum_rows <- do.call(rbind, lapply(results, function(r) c(
    r$dataset_id   %||% "",
    r$tissue       %||% "",
    r$gene_id      %||% "",
    sprintf("%.4f", r$PP4_ABF   %||% NA_real_),
    sprintf("%.4f", r$PP4_SuSiE %||% NA_real_),
    sprintf("%.4f", r$best_PP4  %||% NA_real_),
    r$susie_method %||% ""
)))
colnames(sum_rows) <- c("Dataset","Tissue","Gene","PP4_ABF","PP4_SuSiE","best_PP4","SuSiE detail")
sum_rows <- sum_rows[order(suppressWarnings(-as.numeric(sum_rows[,"best_PP4"]))), , drop=FALSE]
p3 <- make_table_page(sum_rows,
                      sprintf("Coloc summary across %d pair(s) at %s",
                              nrow(sum_rows), opts$locus))
grid.arrange(p3$title, p3$tbl, ncol=1, heights=c(1, 13))

# ─── Pages 4+: One details page per pair ─────────────────────────────────────
# Order by best_PP4 desc so the most interesting pairs come first
# ─── Pages 4+: One details page per pair ─────────────────────────────────────
# Order by best_PP4 desc so the most interesting pairs come first
ord <- order(-sapply(results, function(r) r$best_PP4 %||% 0))
for (i in ord) {
    r <- results[[i]]
    pp <- r$pp_all_abf
    detail_rows <- rbind(
        c("Dataset",      sprintf("%s | Study: %s | Tissue: %s",
                                  r$dataset_id %||% "?",
                                  r$study_label %||% "?",
                                  r$tissue %||% "?")),
        c("Gene",         r$gene_id %||% "?"),
        c("Quant method", r$quant_method %||% "?"),
        c("Coloc method", r$method %||% "?"),
        c("", ""),
        c("GWAS N used",  format(r$gwas_n_used %||% NA_integer_, big.mark=",")),
        c("SNPs in coloc",format(r$n_snps %||% NA_integer_, big.mark=",")),
        c("", ""),
        c("ABF PP(H0)",   if (!is.na(pp["PP.H0.abf"])) sprintf("%.4f", pp["PP.H0.abf"]) else "—"),
        c("ABF PP(H1)",   if (!is.na(pp["PP.H1.abf"])) sprintf("%.4f", pp["PP.H1.abf"]) else "—"),
        c("ABF PP(H2)",   if (!is.na(pp["PP.H2.abf"])) sprintf("%.4f", pp["PP.H2.abf"]) else "—"),
        c("ABF PP(H3)",   if (!is.na(pp["PP.H3.abf"])) sprintf("%.4f", pp["PP.H3.abf"]) else "—"),
        c("ABF PP4",
          if (!is.na(pp["PP.H4.abf"]))
              sprintf("%.4f  *** %s ***", pp["PP.H4.abf"], pp4_interpret(pp["PP.H4.abf"]))
          else "—"),
        c("SuSiEx CSes",  sprintf("GWAS: %s | eQTL: %s",
                                  ifelse(is.na(r$n_gwas_cs %||% NA), "?",
                                         as.character(r$n_gwas_cs)),
                                  ifelse(is.na(r$n_eqtl_cs %||% NA), "?",
                                         as.character(r$n_eqtl_cs)))),
        c("SuSiE best PP4",
          if (!is.na(r$PP4_SuSiE %||% NA_real_))
              sprintf("%.4f [%s] | detail: %s",
                      r$PP4_SuSiE, pp4_interpret(r$PP4_SuSiE),
                      r$susie_method %||% "?")
          else sprintf("— (%s)", r$susie_method %||% "not run")),
        c("Overall best PP4",
          sprintf("%.4f -> %s", r$best_PP4 %||% NA_real_, pp4_interpret(r$best_PP4 %||% NA_real_)))
    )
    colnames(detail_rows) <- c("Parameter","Value")

    title <- textGrob(sprintf("%s | %s | %s",
                              r$dataset_id %||% "?",
                              r$tissue     %||% "?",
                              r$gene_id    %||% "?"),
                      gp=gpar(fontsize=12, fontface="bold"))
    
    tbl   <- tableGrob(detail_rows, rows=NULL,
                       theme=ttheme_minimal(base_size=8.5))
    
    # Load LD Matrix once
    ld_mat <- NULL
    if (nchar(opts$ld) > 0 && file.exists(opts$ld)) {
        ld_mat <- tryCatch(readRDS(opts$ld), error=function(e) NULL)
    }
    
    plt   <- regional_plot(r, ld_mat, opts$lead_snp)
    
    # 1. DRAW TEXT TABLE TO PDF 
    # (Removed the plot from grid.arrange so it doesn't overlap)
    grid.arrange(title, tbl, ncol=1, heights=c(1, 12))
    
    # 2. DRAW PLOT TO PDF 
    # (Calling print() on a ggplot object forces it onto a fresh, clean page)
    print(plt)
    
    # 3. EXPORT HIGH-RES 600 DPI PNG
    png_file <- file.path(opts$results_dir, sprintf("%s.%s.regional_plot.png", r$dataset_id %||% "unknown", r$gene_id %||% "unknown"))
    png(png_file, width=8, height=6, units="in", res=600)
    print(plt)
    dev.off()
}

invisible(dev.off())
cat(sprintf("[generate_locus_report.R] Wrote: %s\n", opts$out))