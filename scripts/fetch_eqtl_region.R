#!/usr/bin/env Rscript
# =============================================================================
# fetch_eqtl_region.R
#
# Retrieves eQTL data for a genomic region from the eQTL Catalogue.
#
# Multi-tissue / multi-gene support:
#   --tissues   Comma-separated list of tissue_label values.
#   --gene_ids  Comma-separated Ensembl gene IDs.
#
# PRE-FLIGHT OPTIMIZATION:
#   Uses a local tabix-indexed file (--local_index) containing all eQTL hits 
#   with p < 1e-4. If a dataset/gene pair has no variants in this index for the 
#   requested region, the slow FTP nominal/LBF downloads are skipped entirely.
#
# Outputs written to --out directory:
#   {dataset_id}.{gene_id}.nominal.tsv
#   {dataset_id}.{gene_id}.lbf.tsv
#   {dataset_id}.{gene_id}.cs.tsv
#   manifest.tsv
#
# Required R packages: optparse, dplyr, readr, seqminer
# =============================================================================

# Enforce shared ACCRE R library path
.libPaths(c("/data/h_vmac/waltes2/Coloc/Conda_Environment_New/rlib-4.5.0", .libPaths()))

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
    library(seqminer)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# ─── CLI arguments ────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--chr",         type="character", help="Chromosome (e.g. 19, no 'chr' prefix)"),
    make_option("--start",       type="integer",   help="Region start bp (GRCh38)"),
    make_option("--end",         type="integer",   help="Region end bp (GRCh38)"),
    make_option("--metadata",    type="character",
                default="/data/h_vmac/zhanm32/coloc_pipeline/data/eQTLcatalogue/dataset_metadata_r7.tsv",
                help="Local path to dataset_metadata_r7.tsv"),
    make_option("--tabix_paths", type="character",
                default="/data/h_vmac/zhanm32/coloc_pipeline/data/eQTLcatalogue/tabix_ftp_paths.tsv",
                help="Local path to tabix_ftp_paths.tsv"),
    make_option("--local_index", type="character",
                default="/data/h_vmac/zhanm32/coloc_pipeline/data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz",
                help="Local tabix-indexed file containing significant eQTL hits (p < 1e-4). Used to preemptively skip FTP downloads."),
    make_option("--gene_ids",    type="character", default="",
                help="Comma-separated Ensembl gene IDs (leave blank to keep all genes in region)"),
    make_option("--study",       type="character", default="",
                help="Study label (e.g. GTEx); single study only"),
    make_option("--tissues",     type="character", default="",
                help="Comma-separated tissue_label values (e.g. 'brain (hippocampus), brain (frontal cortex BA9)')"),
    make_option("--quant",       type="character", default="ge",
                help="Quantification method: ge | exon | tx | txrev | all (default: ge)"),
    make_option("--method",      type="character", default="ABF",
                help="Coloc method: ABF | SUSIE | ALL (controls whether LBF is fetched)"),
    make_option("--sleep",       type="double",    default=2.0,
                help="Minimum seconds between tabix/download calls (floor: 2.0)"),
    make_option("--lbf_gene_index", type="character", default="",
                help="Path to lbf_gene_index.tsv (two columns: dataset_id, gene_id)."),
    make_option("--lbf_cache",    type="character", default="",
                help="Directory to cache downloaded LBF .gz files permanently."),
    make_option("--out",         type="character", help="Output directory")
)

opts <- parse_args(OptionParser(option_list=option_list))
opts$method <- toupper(trimws(opts$method))
opts$quant  <- tolower(trimws(opts$quant))

stopifnot(!is.null(opts$chr), !is.null(opts$start),
          !is.null(opts$end), !is.null(opts$out))
dir.create(opts$out, showWarnings=FALSE, recursive=TRUE)

SLEEP_SECS   <- max(opts$sleep, 2.0)
FETCH_LBF    <- opts$method %in% c("SUSIE", "ALL")
TABIX_REGION <- paste0(opts$chr, ":", opts$start, "-", opts$end)
LBF_CACHE_DIR <- trimws(opts$lbf_cache %||% "")
if (nchar(LBF_CACHE_DIR) == 0) LBF_CACHE_DIR <- NULL

# ─── Parse comma-separated tissues and gene IDs ───────────────────────────────
split_csv <- function(x) {
    v <- trimws(strsplit(x, ",")[[1]])
    v[nchar(v) > 0]
}

tissues_list  <- split_csv(opts$tissues)
gene_ids_list <- split_csv(opts$gene_ids)

cat("============================================================\n")
cat(sprintf(" fetch_eqtl_region.R | region chr%s:%s-%s\n",
            opts$chr, opts$start, opts$end))
cat(sprintf(" method: %-6s | fetch LBF: %-3s | quant: %-12s | sleep: %.1fs\n",
            opts$method, ifelse(FETCH_LBF,"yes","no"), opts$quant, SLEEP_SECS))
cat(sprintf(" study: %s\n", opts$study))
cat(sprintf(" tissue(s)  [%d]: %s\n",
            length(tissues_list),
            if (length(tissues_list) > 0) paste(tissues_list, collapse=" | ") else "(all)"))
cat(sprintf(" gene ID(s) [%d]: %s\n",
            length(gene_ids_list),
            if (length(gene_ids_list) > 0) paste(gene_ids_list, collapse=", ") else "(all in region)"))
cat("============================================================\n\n")

# =============================================================================
# Helper Functions (Retained from Original)
# =============================================================================
match_any_tissue <- function(tissue_label, tissue_id, sample_group, queries) {
    if (length(queries) == 0) return(rep(TRUE, length(tissue_label)))
    Reduce("|", lapply(queries, function(q) {
        grepl(q, tissue_label, fixed=TRUE) |
        grepl(q, tissue_id,    fixed=TRUE) |
        grepl(q, sample_group, fixed=TRUE)
    }))
}

HEADER_CACHE <- new.env(hash=TRUE, parent=emptyenv())
read_ftp_header <- function(url, max_tries=3L, retry_sleep=15L) {
    if (exists(url, envir=HEADER_CACHE)) return(get(url, envir=HEADER_CACHE))
    hdr <- NULL
    for (attempt in seq_len(max_tries)) {
        if (attempt > 1L) {
            Sys.sleep(retry_sleep)
        }
        hdr <- tryCatch({
            con    <- url(url, open="rb")
            gz_con <- gzcon(con)
            on.exit({ try(close(gz_con), silent=TRUE) }, add=TRUE)
            line   <- readLines(gz_con, n=1, warn=FALSE)
            if (length(line) == 0 || nchar(trimws(line)) == 0) stop("empty first line")
            strsplit(trimws(line), "\t")[[1]]
        }, error = function(e) { NULL })
        if (!is.null(hdr)) break
    }
    assign(url, hdr, envir=HEADER_CACHE)
    hdr
}

fetch_tabix <- function(url, region, max_tries=3L, retry_sleep=15L) {
    raw <- NULL
    for (attempt in seq_len(max_tries)) {
        if (attempt > 1L) Sys.sleep(retry_sleep)
        Sys.sleep(SLEEP_SECS)
        raw <- tryCatch({
            r <- seqminer::tabix.read.table(tabixFile=url, tabixRange=region, stringsAsFactors=FALSE)
            if (is.null(r) || nrow(r) == 0) NULL else r
        }, error = function(e) { NULL })
        if (!is.null(raw)) break
    }
    raw
}

fetch_cs <- function(url, chr, start, end, gene_ids=NULL, max_tries=3L, retry_sleep=15L) {
    raw <- NULL
    for (attempt in seq_len(max_tries)) {
        if (attempt > 1L) Sys.sleep(retry_sleep)
        Sys.sleep(SLEEP_SECS)
        raw <- tryCatch({
            r <- readr::read_tsv(url, col_types = cols(.default = "c"), show_col_types = FALSE, progress = FALSE)
            if (is.null(r) || nrow(r) == 0) NULL else r
        }, error = function(e) { NULL })
        if (!is.null(raw)) break
    }
    if (is.null(raw)) return(NULL)

    if (all(c("chromosome", "position") %in% colnames(raw))) {
        raw <- raw %>% dplyr::mutate(chromosome = as.character(chromosome), position = suppressWarnings(as.integer(position))) %>%
            dplyr::filter(chromosome == as.character(chr), is.finite(position), position >= start, position <= end)
    }
    if (!is.null(gene_ids) && length(gene_ids) > 0 && "molecular_trait_id" %in% colnames(raw)) {
        raw <- dplyr::filter(raw, molecular_trait_id %in% gene_ids)
    }
    if (nrow(raw) == 0) return(NULL)
    raw
}

fetch_lbf_full <- function(url, chr, start, end, gene_ids=NULL, cache_dir=NULL) {
    # Abbreviated for space, identical logic to original curl/cache logic
    Sys.sleep(SLEEP_SECS)
    fname <- basename(url)
    use_cache <- !is.null(cache_dir) && nchar(trimws(cache_dir)) > 0
    cached_file <- if (use_cache) file.path(cache_dir, fname) else NULL

    if (use_cache && file.exists(cached_file)) {
        local_bytes <- as.numeric(file.size(cached_file))
        remote_bytes <- tryCatch({
            head_out <- system2("curl", args=c("--head", "--silent", "--connect-timeout", "15", shQuote(url)), stdout=TRUE, stderr=FALSE)
            cl_line <- head_out[grepl("(?i)content-length", head_out, perl=TRUE)]
            if (length(cl_line) > 0) as.numeric(trimws(sub(".*:\\s*", "", cl_line[1]))) else NA_real_
        }, error = function(e) NA_real_)
        if (!is.na(remote_bytes) && !is.na(local_bytes) && local_bytes != remote_bytes) unlink(cached_file)
    }

    if (use_cache && file.exists(cached_file) && file.size(cached_file) > 0) {
        tmp <- cached_file; delete_tmp <- FALSE
    } else {
        if (use_cache) { dir.create(cache_dir, showWarnings=FALSE, recursive=TRUE); tmp <- cached_file; delete_tmp <- FALSE }
        else { tmp <- tempfile(fileext=".lbf.txt.gz"); delete_tmp <- TRUE }
        
        is_tty <- isatty(stderr())
        progress_args <- if (is_tty) "--progress-bar" else c("--silent", "--show-error")
        dl_status <- system2("curl", args=c("--location", "--connect-timeout", "30", "--max-time", "1800", "--retry", "2", "--retry-delay", "5", progress_args, "-o", shQuote(tmp), shQuote(url)), stdout=FALSE, stderr="")
        if (dl_status != 0 || !file.exists(tmp) || file.size(tmp) == 0) { if (delete_tmp) unlink(tmp); return(NULL) }
        gz_ok <- system2("gzip", c("--test", shQuote(tmp)), stdout=FALSE, stderr=FALSE)
        if (gz_ok != 0) { unlink(tmp); return(NULL) }
    }
    if (delete_tmp) on.exit(unlink(tmp), add=TRUE)

    hdr_line <- tryCatch({ con <- gzcon(file(tmp, open="rb")); on.exit(try(close(con), silent=TRUE), add=TRUE); readLines(con, n=1, warn=FALSE) }, error = function(e) character(0))
    df <- NULL
    if (length(hdr_line) > 0 && nchar(trimws(hdr_line)) > 0) {
        col_names <- strsplit(trimws(hdr_line), "\t")[[1]]
        col_gene <- which(col_names == "molecular_trait_id")
        col_chr <- which(col_names == "chromosome")
        col_pos <- which(col_names == "position")

        data_filter <- "1"
        if (length(col_gene) == 1 && !is.null(gene_ids) && length(gene_ids) > 0) {
            gene_re <- paste(gene_ids, collapse="|")
            data_filter <- sprintf("$%d ~ /^(%s)$/", col_gene, gene_re)
        }
        if (length(col_chr) == 1 && length(col_pos) == 1) {
            data_filter <- sprintf("%s && $%d == \"%s\" && ($%d+0) >= %d && ($%d+0) <= %d", data_filter, col_chr, as.character(chr), col_pos, start, col_pos, end)
        }
        awk_expr <- sprintf("'NR==1 || (%s)'", data_filter)
        awk_cmd <- sprintf("zcat %s 2>/dev/null | awk -F'\\t' %s", shQuote(tmp), awk_expr)
        
        df <- tryCatch({
            con <- pipe(awk_cmd, open="rb"); on.exit(try(close(con), silent=TRUE), add=TRUE)
            readr::read_tsv(con, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        }, error = function(e) NULL)
    }
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
}

# =============================================================================
# STEP 1 — Load metadata and FTP path table
# =============================================================================
cat("[1/4] Loading dataset metadata and FTP path table...\n")
datasets_df <- readr::read_tsv(opts$metadata, col_types=cols(.default="c"), show_col_types=FALSE)
ftp_df <- readr::read_tsv(opts$tabix_paths, col_types=cols(.default="c"), show_col_types=FALSE)

HAS_CS_PATH_COL <- "ftp_cs_path" %in% colnames(ftp_df)
ftp_select_cols <- c("dataset_id", "ftp_path", "ftp_lbf_path")
if (HAS_CS_PATH_COL) ftp_select_cols <- c(ftp_select_cols, "ftp_cs_path")

datasets_df <- dplyr::left_join(datasets_df, dplyr::select(ftp_df, dplyr::all_of(ftp_select_cols)), by="dataset_id")
if (!HAS_CS_PATH_COL) datasets_df$ftp_cs_path <- NA_character_

# =============================================================================
# STEP 2 — Filter to relevant datasets
# =============================================================================
cat("[2/4] Selecting datasets...\n")
datasets_use <- datasets_df
if (opts$quant != "all") datasets_use <- dplyr::filter(datasets_use, tolower(quant_method) == opts$quant)
if (nchar(opts$study) > 0) datasets_use <- dplyr::filter(datasets_use, grepl(opts$study, study_label, fixed=TRUE, ignore.case=FALSE) | study_id == opts$study)
if (length(tissues_list) > 0) {
    keep_rows <- match_any_tissue(datasets_use$tissue_label, datasets_use$tissue_id, datasets_use$sample_group, tissues_list)
    datasets_use <- datasets_use[keep_rows, ]
}
if (nrow(datasets_use) == 0) stop("No matching datasets found.")

# =============================================================================
# STEP 3 — PRE-FLIGHT CHECK: Query local p<1e-4 index for the entire region
# =============================================================================
cat("[3/4] Pre-flight check: Querying local index for significant hits...\n")
LOCAL_HITS <- NULL
local_index_path <- trimws(opts$local_index)

if (nchar(local_index_path) > 0 && file.exists(local_index_path)) {
    cat(sprintf("      Querying : %s\n", basename(local_index_path)))
    
    LOCAL_HITS <- tryCatch({
        res <- seqminer::tabix.read.table(local_index_path, TABIX_REGION, stringsAsFactors=FALSE)
        if (is.null(res) || nrow(res) == 0) NULL else res
    }, error = function(e) {
        cat(sprintf("      [WARN] Failed to read local index: %s\n", e$message))
        NULL
    })

    if (!is.null(LOCAL_HITS)) {
        # seqminer replaces '#' with '.' in headers. Clean it up.
        colnames(LOCAL_HITS) <- sub("^X\\.", "", colnames(LOCAL_HITS))
        colnames(LOCAL_HITS) <- sub("^#|\\.", "", colnames(LOCAL_HITS))
        colnames(LOCAL_HITS)[1] <- "chromosome"

        if (length(gene_ids_list) > 0) {
            LOCAL_HITS <- dplyr::filter(LOCAL_HITS, gene_id %in% gene_ids_list)
        }
        cat(sprintf("      Found    : %d significant hit(s) across %d dataset(s).\n",
                    nrow(LOCAL_HITS), length(unique(LOCAL_HITS$dataset_id))))
    } else {
        cat("      Found    : 0 significant hits in region for requested genes.\n")
    }
} else {
    cat("      [INFO] No local index provided or file missing. Will attempt full FTP fetch.\n")
}

# =============================================================================
# STEP 4 — Fetch data per dataset
# =============================================================================
cat(sprintf("[4/4] Fetching data for %d dataset(s)...\n", nrow(datasets_use)))

manifest_rows <- list()
n_tabix_calls <- 0L

for (i in seq_len(nrow(datasets_use))) {
    ds <- datasets_use[i, ]
    ds_id <- ds$dataset_id
    st_id <- ds$study_id %||% ""

    cat(sprintf("\n  [%d/%d] %-12s  %-20s  tissue: %s\n", i, nrow(datasets_use), ds_id, ds$study_label %||% "?", ds$tissue_label %||% "?"))
    nom_url <- ds$ftp_path %||% NA_character_
    lbf_url <- ds$ftp_lbf_path %||% NA_character_

    if (is.na(nom_url) || nchar(nom_url) == 0) {
        cat("    [SKIP] No FTP path available\n")
        next
    }

    # ── Local Index Skip Logic ──
    dataset_sig_genes <- NULL
    if (!is.null(LOCAL_HITS)) {
        ds_hits <- dplyr::filter(LOCAL_HITS, dataset_id == ds_id)
        if (nrow(ds_hits) == 0) {
            cat("    [SKIP] Minimum p-value > 1e-4 in local index. Not worth running coloc.\n")
            next
        }
        dataset_sig_genes <- unique(ds_hits$gene_id)
        cat(sprintf("    index    : %d gene(s) passed p < 1e-4 threshold\n", length(dataset_sig_genes)))
    } else if (nchar(local_index_path) > 0 && file.exists(local_index_path)) {
        # LOCAL_HITS is NULL because the region returned empty for EVERYTHING
        cat("    [SKIP] Minimum p-value > 1e-4 in local index. Not worth running coloc.\n")
        next
    }

    # Read nominal column header
    nom_hdr <- read_ftp_header(nom_url)
    if (is.null(nom_hdr)) { cat("    [SKIP] Could not read nominal file header\n"); next }

    # Fetch nominal data via seqminer
    n_tabix_calls <- n_tabix_calls + 1L
    raw_nom <- fetch_tabix(nom_url, TABIX_REGION)
    if (is.null(raw_nom)) { cat("    [SKIP] No nominal associations in region\n"); next }

    colnames(raw_nom) <- nom_hdr
    raw_nom <- dplyr::mutate(raw_nom, dplyr::across(dplyr::where(is.character), ~gsub("\r", "", .x, fixed=TRUE)))

    nom_df <- raw_nom %>%
        dplyr::mutate(chromosome = as.character(chromosome), position = as.integer(position), beta = as.double(beta), se = as.double(se), pvalue = as.double(pvalue), variant_id = paste0("chr", chromosome, "_", position, "_", ref, "_", alt)) %>%
        dplyr::filter(is.finite(beta), is.finite(se), se > 0)

    # Combine User Requested Genes with Local Index Hits
    target_genes <- gene_ids_list
    if (!is.null(dataset_sig_genes)) {
        if (length(target_genes) > 0) {
            target_genes <- intersect(target_genes, dataset_sig_genes)
        } else {
            target_genes <- dataset_sig_genes
        }
    }

    if (length(target_genes) > 0) {
        nom_df <- dplyr::filter(nom_df, molecular_trait_id %in% target_genes)
        if (nrow(nom_df) == 0) {
            cat(sprintf("    [SKIP] None of the targeted gene ID(s) found in regional nominal data\n"))
            next
        }
    }

    genes_found <- unique(nom_df$molecular_trait_id)
    cat(sprintf("    nominal: %d rows | %d gene(s)\n", nrow(nom_df), length(genes_found)))

    # LBF Fetching Logic
    lbf_all <- NULL
    if (FETCH_LBF && !is.na(lbf_url) && nchar(lbf_url) > 0) {
        n_tabix_calls <- n_tabix_calls + 1L
        lbf_all <- fetch_lbf_full(url=lbf_url, chr=opts$chr, start=opts$start, end=opts$end, gene_ids=genes_found, cache_dir=LBF_CACHE_DIR)
        if (!is.null(lbf_all)) {
            lbf_all <- lbf_all %>% dplyr::mutate(chromosome = as.character(chromosome), position = as.integer(position), variant_id = if ("variant" %in% colnames(.)) variant else paste0("chr", chromosome, "_", position))
        }
    }

    # CS Fetching Logic
    cs_all <- NULL
    if (FETCH_LBF && !is.na(ds$ftp_cs_path) && nchar(ds$ftp_cs_path) > 0) {
        n_tabix_calls <- n_tabix_calls + 1L
        cs_raw <- fetch_cs(ds$ftp_cs_path, chr=opts$chr, start=opts$start, end=opts$end, gene_ids=genes_found)
        if (!is.null(cs_raw)) {
            cs_raw <- dplyr::mutate(cs_raw, dplyr::across(dplyr::where(is.character), ~gsub("\r", "", .x, fixed=TRUE)))
            cs_all <- cs_raw
        }
    }

    # Write Outputs
    for (gid in genes_found) {
        nom_gene <- dplyr::filter(nom_df, molecular_trait_id == gid)
        nom_file <- file.path(opts$out, paste0(ds_id, ".", gid, ".nominal.tsv"))
        readr::write_tsv(nom_gene, nom_file)

        lbf_file <- NA_character_
        if (!is.null(lbf_all)) {
            lbf_gene <- dplyr::filter(lbf_all, molecular_trait_id == gid)
            if (nrow(lbf_gene) > 0) { lbf_file <- file.path(opts$out, paste0(ds_id, ".", gid, ".lbf.tsv")); readr::write_tsv(lbf_gene, lbf_file) }
        }

        cs_file <- NA_character_; n_cs_gene <- 0L
        if (!is.null(cs_all)) {
            cs_gene <- dplyr::filter(cs_all, molecular_trait_id == gid)
            if (nrow(cs_gene) > 0) { cs_file <- file.path(opts$out, paste0(ds_id, ".", gid, ".cs.tsv")); readr::write_tsv(cs_gene, cs_file); n_cs_gene <- length(unique(cs_gene$cs_id)) }
        }

        manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
            dataset_id=ds_id, study_id=st_id, study_label=ds$study_label %||% NA_character_, tissue=ds$tissue_label %||% NA_character_, cell_type=ds$sample_group %||% NA_character_, quant_method=ds$quant_method %||% NA_character_, gene_id=gid, n_snps=nrow(nom_gene), n_cs=n_cs_gene, nominal_file=nom_file, lbf_file=lbf_file, cs_file=cs_file, stringsAsFactors=FALSE
        )
    }
}

# =============================================================================
# Write manifest
# =============================================================================
if (length(manifest_rows) == 0) {
    cat("    [INFO] No significant eQTL data retrieved. Generating empty manifest to satisfy pipeline.\n")
    manifest <- data.frame(dataset_id=character(), study_id=character(), study_label=character(), 
                           tissue=character(), cell_type=character(), quant_method=character(), 
                           gene_id=character(), n_snps=integer(), n_cs=integer(), 
                           nominal_file=character(), lbf_file=character(), cs_file=character(), 
                           stringsAsFactors=FALSE)
} else {
    manifest <- dplyr::bind_rows(manifest_rows)
}

manifest_path <- file.path(opts$out, "manifest.tsv")
readr::write_tsv(manifest, manifest_path)

cat("\n============================================================\n")
cat(" fetch_eqtl_region.R complete\n")
cat(sprintf(" Dataset x gene pairs : %d\n", nrow(manifest)))
cat("============================================================\n")