#!/usr/bin/env Rscript
# =============================================================================
# build_eqtl_index.R
#
# Processes ONE ge dataset from the eQTL Catalogue and appends its
# sub-threshold associations (pvalue < p_thresh) to a per-dataset TSV file.
#
# Designed to be run as a SLURM array job — each task receives a unique
# --dataset_idx (1-based integer) that selects one dataset from the 280 ge
# entries in tabix_ftp_paths.tsv.
#
# Workflow per task:
#   1. Join tabix_ftp_paths.tsv + dataset_metadata_r7.tsv; filter to ge,
#      sort by dataset_id for reproducible indexing.
#   2. Select the dataset at --dataset_idx.
#   3. Download the full .all.tsv.gz to a temp file (DELETE AFTER USE).
#   4. Stream through in chunks of 500 k rows; filter pvalue < p_thresh.
#   5. Write filtered rows to {out_dir}/parts/{dataset_id}.tsv (no header).
#   6. Delete the temp download file.
#   7. Write {out_dir}/status/{dataset_id}.done on success, .failed on error.
#
# The per-dataset .tsv files are merged, sorted, bgzipped, and tabix-indexed
# by merge_eqtl_index.slurm AFTER the array completes.
#
# Required R packages: optparse, dplyr, readr
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
})

# ─── CLI arguments ────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--dataset_idx",  type="integer",   help="1-based index into sorted ge dataset list [required]"),
    make_option("--tabix_paths",  type="character",
                default="/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/tabix_ftp_paths.tsv",
                help="Path to tabix_ftp_paths.tsv"),
    make_option("--metadata",     type="character",
                default="/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/dataset_metadata_r7.tsv",
                help="Path to dataset_metadata_r7.tsv"),
    make_option("--out_dir",      type="character",
                default="/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/eqtl_index_build",
                help="Root output directory (parts/ and status/ subdirs created here)"),
    make_option("--p_thresh",     type="double",    default=1e-4,
                help="Nominal p-value threshold for inclusion [default: 1e-4]"),
    make_option("--chunk_size",   type="integer",   default=500000L,
                help="Rows per read chunk [default: 500000]"),
    make_option("--timeout_secs", type="integer",   default=3600L,
                help="Download timeout in seconds [default: 3600]")
)
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opts$dataset_idx))

# ─── Output directories ───────────────────────────────────────────────────────
parts_dir  <- file.path(opts$out_dir, "parts")
status_dir <- file.path(opts$out_dir, "status")
dir.create(parts_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(status_dir, showWarnings = FALSE, recursive = TRUE)

# ─── Helpers ─────────────────────────────────────────────────────────────────
write_status <- function(dataset_id, success, message = "") {
    ext  <- if (success) ".done" else ".failed"
    path <- file.path(status_dir, paste0(dataset_id, ext))
    writeLines(c(format(Sys.time()), message), path)
}

abort <- function(dataset_id, msg) {
    message("[ERROR] ", msg -
    write_status(dataset_id, success = FALSE, message = msg)
    quit(status = 1, save = "no")
}

# ─── STEP 1: Load and join metadata + FTP paths ───────────────────────────────
cat("============================================================\n")
cat(sprintf(" build_eqtl_index.R | task index: %d\n", opts$dataset_idx))
cat(sprintf(" p_thresh : %.0e  |  chunk_size: %d\n", opts$p_thresh, opts$chunk_size))
cat("============================================================\n\n")

cat("[1/6] Loading metadata and FTP path table...\n")

if (!file.exists(opts$metadata))
    stop("Metadata file not found: ", opts$metadata)
if (!file.exists(opts$tabix_paths))
    stop("tabix_ftp_paths not found: ", opts$tabix_paths)

meta_df <- readr::read_tsv(opts$metadata,    col_types = cols(.default = "c"), show_col_types = FALSE)
ftp_df  <- readr::read_tsv(opts$tabix_paths, col_types = cols(.default = "c"), show_col_types = FALSE)

# Join metadata + FTP paths on dataset_id
datasets_df <- dplyr::left_join(
    meta_df,
    dplyr::select(ftp_df, dataset_id, ftp_path),
    by = "dataset_id"
)

# Filter to ge datasets with a valid FTP path, sort for reproducible indexing
ge_datasets <- datasets_df %>%
    dplyr::filter(tolower(trimws(quant_method)) == "ge",
                  !is.na(ftp_path),
                  nchar(trimws(ftp_path)) > 0) %>%
    dplyr::arrange(dataset_id)

cat(sprintf("      %d ge datasets with FTP paths available\n", nrow(ge_datasets)))

if (opts$dataset_idx < 1 || opts$dataset_idx > nrow(ge_datasets)) {
    stop(sprintf("--dataset_idx %d is out of range (1-%d)",
                 opts$dataset_idx, nrow(ge_datasets)))
}

# ─── STEP 2: Select this task's dataset ──────────────────────────────────────
cat(sprintf("[2/6] Selecting dataset at index %d...\n", opts$dataset_idx))

ds         <- ge_datasets[opts$dataset_idx, ]
ds_id      <- ds$dataset_id
study_lbl  <- ds$study_label  %||% ""
tissue_lbl <- ds$tissue_label %||% ""
ftp_url    <- ds$ftp_path

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nchar(trimws(a)) > 0) a else b

cat(sprintf("      dataset_id   : %s\n", ds_id))
cat(sprintf("      study_label  : %s\n", study_lbl))
cat(sprintf("      tissue_label : %s\n", tissue_lbl))
cat(sprintf("      ftp_url      : %s\n", ftp_url))

# Skip if output already exists (allows resuming)
out_part   <- file.path(parts_dir,  paste0(ds_id, ".tsv"))
done_flag  <- file.path(status_dir, paste0(ds_id, ".done"))
if (file.exists(done_flag) && file.exists(out_part)) {
    cat(sprintf("      [SKIP] %s already processed (%s exists)\n", ds_id, done_flag))
    quit(status = 0, save = "no")
}

# ─── STEP 3: Download full nominal file to temp ───────────────────────────────
cat(sprintf("[3/6] Downloading %s to temp file...\n", basename(ftp_url)))

temp_file <- tempfile(pattern = paste0(ds_id, "_"), fileext = ".tsv.gz")
options(timeout = opts$timeout_secs)

dl_ok <- tryCatch({
    download.file(ftp_url, destfile = temp_file, method = "libcurl",
                  quiet = FALSE, mode = "wb")
    TRUE
}, error = function(e) {
    message("  [WARN] Download failed: ", e$message)
    FALSE
})

if (!dl_ok || !file.exists(temp_file) || file.size(temp_file) < 1000L) {
    if (file.exists(temp_file)) unlink(temp_file)
    abort(ds_id, sprintf("Download failed for %s (%s)", ds_id, ftp_url))
}

cat(sprintf("      Downloaded: %.1f MB\n", file.size(temp_file) / 1e6))

# ─── STEP 4: Stream through in chunks, filter pvalue < p_thresh ──────────────
cat(sprintf("[4/6] Filtering to pvalue < %.0e in chunks of %d rows...\n",
            opts$p_thresh, opts$chunk_size))

#  Columns we keep in the index (minimal footprint):
#    chromosome, position — for tabix positional lookup
#    gene_id              — Ensembl gene ID (same as molecular_trait_id for ge)
#    pvalue, beta, maf    — signal summary
#  dataset_id, study_label, tissue_label added from metadata

KEEP_COLS <- c("chromosome", "position", "gene_id", "pvalue", "beta", "maf")

collected  <- list()
n_rows_in  <- 0L
n_rows_out <- 0L

tryCatch({
    readr::read_tsv_chunked(
        file           = temp_file,
        callback       = readr::DataFrameCallback$new(function(chunk, pos) {
            n_rows_in <<- n_rows_in + nrow(chunk)

            # Coerce key columns; suppressWarnings because some rows may have NA
            chunk <- chunk %>%
                dplyr::mutate(
                    position = suppressWarnings(as.integer(position)),
                    pvalue   = suppressWarnings(as.double(pvalue)),
                    beta     = suppressWarnings(as.double(beta)),
                    maf      = suppressWarnings(as.double(maf))
                ) %>%
                dplyr::filter(is.finite(pvalue), pvalue < opts$p_thresh) %>%
                dplyr::mutate(
                    chromosome   = as.character(chromosome),
                    dataset_id   = ds_id,
                    study_label  = study_lbl,
                    tissue_label = tissue_lbl,
                    # Use gene_id column; fall back to molecular_trait_id for ge
                    gene_id = dplyr::coalesce(gene_id, molecular_trait_id)
                )

            # Reorder and select output columns
            out_cols <- c("chromosome", "position", "dataset_id",
                          "study_label", "tissue_label", "gene_id",
                          "pvalue", "beta", "maf")
            out_cols <- intersect(out_cols, colnames(chunk))
            chunk <- dplyr::select(chunk, dplyr::all_of(out_cols))

            if (nrow(chunk) > 0) {
                n_rows_out <<- n_rows_out + nrow(chunk)
                collected[[length(collected) + 1L]] <<- chunk
            }
        }),
        chunk_size     = opts$chunk_size,
        col_types      = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress       = FALSE
    )
}, error = function(e) {
    unlink(temp_file)
    abort(ds_id, sprintf("Chunked read failed for %s: %s", ds_id, e$message))
})

cat(sprintf("      Rows read   : %d\n", n_rows_in))
cat(sprintf("      Rows kept   : %d  (pvalue < %.0e)\n", n_rows_out, opts$p_thresh))

# ─── STEP 5: Write per-dataset output (no header — merge script adds it) ──────
cat(sprintf("[5/6] Writing %s...\n", out_part))

if (length(collected) > 0) {
    result_df <- dplyr::bind_rows(collected)
    # Remove any rows with NA in key positional columns
    result_df <- result_df[!is.na(result_df$chromosome) & !is.na(result_df$position), ]
    readr::write_tsv(result_df, out_part, col_names = FALSE)
    cat(sprintf("      Written: %d rows to %s\n", nrow(result_df), basename(out_part)))
} else {
    # Write empty file so merge script can detect it
    file.create(out_part)
    cat(sprintf("      No associations below threshold — empty file written\n"))
}

# ─── STEP 6: Delete temp download ─────────────────────────────────────────────
cat("[6/6] Deleting temp download file...\n")
unlink(temp_file)
cat(sprintf("      Deleted: %s\n", basename(temp_file)))

# ─── Done ─────────────────────────────────────────────────────────────────────
write_status(ds_id, success = TRUE,
             message = sprintf("%d rows written", n_rows_out))

cat(sprintf("\n[DONE] %s — %d associations at p < %.0e\n",
            ds_id, n_rows_out, opts$p_thresh))

