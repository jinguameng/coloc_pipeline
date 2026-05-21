#!/usr/bin/env Rscript
# =============================================================================
# find_failed_indices.R
#
# Reads *.failed files from the status directory, reconstructs the same
# sorted ge dataset list used by build_eqtl_index.R, and prints the
# 1-based array indices needed to resubmit just the failed tasks.
# Optionally filters for files modified AFTER a specific time.
#
# Usage (run from anywhere):
#   # To read ALL failed files:
#   Rscript find_failed_indices.R
#
#   # To read only files modified AFTER a specific time:
#   Rscript find_failed_indices.R "2026-04-06 09:28:00"
#
# Output:
#   A table of failed dataset IDs with their array indices, plus a ready-to-use
#   sbatch --array=... argument to copy-paste for resubmission.
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

TABIX_PATHS <- "/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/tabix_ftp_paths.tsv"
METADATA    <- "/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/dataset_metadata_r7.tsv"
STATUS_DIR  <- "/data/h_vmac/zhanm32/colocpipe/data/eQTLcatalogue/eqtl_index_build/status"
SLURM_SCRIPT <- "/data/h_vmac/zhanm32/colocpipe/scripts/build_eqtl_index_array.slurm"

# ── 0. Parse command line arguments for time cutoff
args <- commandArgs(trailingOnly = TRUE)
time_cutoff <- NULL

if (length(args) > 0) {
    # Combine arguments in case the user didn't quote the date string
    time_str <- paste(args, collapse = " ")
    time_cutoff <- as.POSIXct(time_str)
    
    if (is.na(time_cutoff)) {
        stop("Invalid time format provided. Please use standard format, e.g., 'YYYY-MM-DD HH:MM:SS'\nExample: Rscript find_failed_indices.R '2026-04-06 09:28:00'")
    }
    cat(sprintf("Filtering for .failed files modified AFTER: %s\n\n", time_cutoff))
} else {
    cat("No time cutoff provided. Reading ALL .failed files.\n\n")
}

# ── 1. Reconstruct the sorted ge dataset list (identical to build_eqtl_index.R)
meta_df <- readr::read_tsv(METADATA,    col_types = cols(.default = "c"), show_col_types = FALSE)
ftp_df  <- readr::read_tsv(TABIX_PATHS, col_types = cols(.default = "c"), show_col_types = FALSE)

ge_datasets <- dplyr::left_join(meta_df,
                                 dplyr::select(ftp_df, dataset_id, ftp_path),
                                 by = "dataset_id") %>%
    dplyr::filter(tolower(trimws(quant_method)) == "ge",
                  !is.na(ftp_path),
                  nchar(trimws(ftp_path)) > 0) %>%
    dplyr::arrange(dataset_id) %>%
    dplyr::mutate(array_index = dplyr::row_number())

cat(sprintf("Total ge datasets in sorted list: %d\n\n", nrow(ge_datasets)))

# ── 2. Read failed dataset IDs from the status directory
# We need full.names = TRUE to check file modification times
failed_files_full <- list.files(STATUS_DIR, pattern = "\\.failed$", full.names = TRUE)

if (length(failed_files_full) == 0) {
    cat("No .failed files found in:\n  ", STATUS_DIR, "\n")
    quit(status = 0, save = "no")
}

# Apply time filter if a valid time was provided
if (!is.null(time_cutoff)) {
    # Get file metadata, specifically modification time (mtime)
    f_info <- file.info(failed_files_full)
    
    # Filter paths where mtime is strictly greater than the cutoff
    failed_files_full <- failed_files_full[f_info$mtime >= time_cutoff]
    
    # Check for NA just in case file.info failed on some files
    failed_files_full <- failed_files_full[!is.na(failed_files_full)]
    
    if (length(failed_files_full) == 0) {
        cat(sprintf("No .failed files found that were modified after %s\n", time_cutoff))
        quit(status = 0, save = "no")
    }
}

# Extract just the filenames and remove the extension
failed_files <- basename(failed_files_full)
failed_ids <- sub("\\.failed$", "", failed_files)
cat(sprintf("Failed datasets to resubmit: %d\n\n", length(failed_ids)))

# ── 3. Look up their indices
result <- data.frame(dataset_id = failed_ids, stringsAsFactors = FALSE) %>%
    dplyr::left_join(dplyr::select(ge_datasets, dataset_id, array_index),
                     by = "dataset_id") %>%
    dplyr::arrange(array_index)

# Flag any IDs not found in the ge list (shouldn't happen, but worth checking)
not_found <- result[is.na(result$array_index), "dataset_id"]
if (length(not_found) > 0) {
    cat("[WARN] These dataset IDs were not found in the ge list:\n")
    cat(paste("  ", not_found, collapse = "\n"), "\n\n")
}

found <- result[!is.na(result$array_index), ]

# ── 4. Print the table
cat(sprintf("%-15s  %s\n", "dataset_id", "array_index"))
cat(strrep("-", 30), "\n")
for (i in seq_len(nrow(found))) {
    cat(sprintf("%-15s  %d\n", found$dataset_id[i], found$array_index[i]))
}

# ── 5. Print the ready-to-use sbatch command
indices_str <- paste(found$array_index, collapse = ",")
cat("\n────────────────────────────────────────────────────────────\n")
cat("Resubmit command (copy-paste):\n\n")
cat(sprintf("  sbatch --array=%s %s\n", indices_str, SLURM_SCRIPT))
cat("────────────────────────────────────────────────────────────\n")