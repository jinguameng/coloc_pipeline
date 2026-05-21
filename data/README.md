# Data Directory Overview

The `data/` directory contains all the static reference panels, metadata, and highly optimized pre-computed indexes required by the `coloc_pipeline`. These files allow the pipeline to perform cross-ancestry LD lookups and dynamically discover eQTLs without downloading massive datasets repeatedly from EBI FTP servers.

## Directory Structure & Key Files

### `1KG/` (1000 Genomes Reference)
Contains population-specific sample lists used to subset the GRCh38 1000 Genomes PLINK reference when computing regional LD matrices for the summary plots.
* **`*.keep` (AFR, AMR, EAS, EUR, SAS):** Text files containing the FID/IID of samples belonging to each superpopulation. This allows PLINK2 to quickly isolate the exact ancestry cohort needed for LD calculation without requiring pre-split `.bed` files.

### `eQTLcatalogue/` (eQTL Lookups & Indexes)
Contains the metadata and local indexes that power the pipeline's "Discovery Mode" and FTP download optimizations.

* **`dataset_metadata_r7.tsv` & `tabix_ftp_paths.tsv`**: Official eQTL Catalogue releases mapping `dataset_id`s to their corresponding studies, tissues (e.g., *brain (hippocampus)*), and remote FTP URLs.
* **`eqtl_index_ge_p1e4.tsv.gz` (and `.tbi`)**: A comprehensive, tabix-indexed database containing every nominal eQTL association across all gene expression (ge) datasets that passed a nominal *p* < 1e-4 threshold. The pipeline's `fetch_eqtl` checkpoint queries this index locally to instantly discover significant genes in a region before initiating any network downloads.
* **`lbf_tasklist_ge.tsv`**: An intermediate task list mapping datasets to their remote LBF (Log Bayes Factor) URLs. Used exclusively by the SLURM array jobs to coordinate the building of the LBF index.
* **`lbf_gene_index.tsv`**: A lightweight, merged lookup table identifying exactly which genes were successfully fine-mapped (have SuSiE LBF data) in which dataset. The pipeline uses this to skip massive LBF file downloads if the target gene is not present in the dataset.

## How these files were generated
To keep the main pipeline lightweight, the pre-computed indexes (`eqtl_index_ge_p1e4.tsv.gz` and `lbf_gene_index.tsv`) are provided as static files. 

If you need to rebuild these indexes from scratch (e.g., for a new eQTL Catalogue release), the backend generation scripts (SLURM arrays and R/Bash parsers) are preserved in the `scripts/index_builders/` directory.
