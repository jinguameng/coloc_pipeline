# eQTL Index Build

Pre-builds a local tabix-indexed eQTL summary statistics index from the eQTL
Catalogue FTP. This index is queried by `fetch_eqtl_region.R` during each
colocpipe run to retrieve nominal association statistics for the GWAS region
of interest — replacing live FTP downloads with near-instant local tabix
queries.

The index covers **nominal p-values** (`.all.tsv.gz` files) only. SuSiE LBF
data is handled separately by the `lbf_index_build/` pipeline.

---

## Scripts — in execution order

---

### 1. `build_eqtl_index_array.slurm`

**What it does**  
SLURM array job. Each task processes one dataset from the eQTL Catalogue:
downloads the nominal summary statistics file, filters to associations below
a p-value threshold (default p < 1×10⁻⁴), retains only the columns needed
for colocpipe, and writes a compressed part file to a staging directory.

One task = one dataset. With 280 ge-quantification datasets, this runs as a
280-task array. Tasks are parallelised (recommended: `%20` concurrent) and
each writes an independent part file, so failures can be retried without
re-running the whole array.

**Run once — before any colocpipe analysis.**

```bash
sbatch --array=1-280%20 scripts/eqtl_index_build/build_eqtl_index_array.slurm
```

**Outputs (per dataset, in staging directory)**
- `{dataset_id}_ge_p1e4.part.tsv.gz` — filtered, compressed part file
- `{dataset_id}_ge_p1e4.part.done` — empty sentinel written on success

---

### 2. `find_failed_indices.R`

**What it does**  
Scans the staging directory for datasets that did not produce a `.done`
sentinel file and prints the corresponding SLURM array task IDs. Run this
after `build_eqtl_index_array.slurm` completes to identify any tasks that
need to be resubmitted.

```bash
Rscript scripts/eqtl_index_build/find_failed_indices.R
```

**Output**  
Prints to console:
- List of failed dataset IDs
- Ready-to-paste `sbatch --array=...` command for resubmission

Resubmit failed tasks, then re-run `find_failed_indices.R` until all 280
show as complete before proceeding to step 3.

---

### 3. `merge_eqtl_index.slurm`

**What it does**  
Merges all 280 part files into a single genome-wide index file, sorts by
chromosome and position (required for tabix), compresses with bgzip, and
builds a tabix index (`.tbi`).

Run as a single SLURM job (not an array) after all part files are confirmed
complete by `find_failed_indices.R`.

```bash
sbatch scripts/eqtl_index_build/merge_eqtl_index.slurm
```

**Outputs**
- `data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz` — merged, sorted, bgzipped index
- `data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz.tbi` — tabix index

This file is the permanent index used by `fetch_eqtl_region.R` and
`query_eqtl_index.R` for all future colocpipe runs.

---

### 4. `build_eqtl_index.R`

**What it does**  
The R script called by each task inside `build_eqtl_index_array.slurm`.
Handles the per-dataset logic: reads the nominal FTP file via
`seqminer::tabix.read.table()`, filters rows, selects columns, and writes
the part file.

**Not called directly by users.** It is invoked automatically by the SLURM
array job. Documented here for reference.

---

### 5. `query_eqtl_index.R`

**What it does**  
Interactive query tool. Given a genomic region (or a lead SNP + window),
queries the merged index to identify which genes and tissues have eQTL
associations in that region. Outputs a formatted summary table and suggests
the exact values to paste into `inputs.txt` for a colocpipe run.

Use this before running colocpipe to decide which genes and tissues to test.

```bash
# Query by lead SNP and window
Rscript scripts/eqtl_index_build/query_eqtl_index.R \
    --index  data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz \
    --snp    rs429358 \
    --window 500 \
    --ref    /data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt \
    --out    results/rs429358_eqtl_hits.tsv

# Query by explicit coordinates
Rscript scripts/eqtl_index_build/query_eqtl_index.R \
    --index  data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz \
    --chr    19 \
    --start  44658684 \
    --end    45158684
```

**Key options**

| Option | Description |
|--------|-------------|
| `--index` | Path to the merged bgzipped index (from step 3) |
| `--snp` | Lead SNP rsID — script looks up CHR/BP from reference file |
| `--window` | Window size in KB centred on lead SNP (default: 500) |
| `--ref` | Space-delimited reference file with columns CHR SNP BP |
| `--chr / --start / --end` | Alternative to `--snp`: explicit region coordinates |
| `--out` | Optional TSV to write full hit table |
| `--p_thresh` | Max p-value to display (default: 1×10⁻⁴) |
| `--top_n` | Show only top N gene-tissue combinations by min p-value |

**Output**  
Console table with columns: Gene, Study, Tissue, N_hits, Min_p, Lead_beta,
Lead_MAF, Lead_pos — sorted by minimum p-value. Followed by a ready-to-paste
block for `inputs.txt`.

---

## One-time setup summary

```bash
# Step 1: submit array job (one task per dataset)
sbatch --array=1-280%20 scripts/eqtl_index_build/build_eqtl_index_array.slurm

# Step 2: check for failures and resubmit if needed
Rscript scripts/eqtl_index_build/find_failed_indices.R

# Step 3: merge all part files into the final index
sbatch scripts/eqtl_index_build/merge_eqtl_index.slurm

# Step 4 (ongoing): query the index before each colocpipe run
Rscript scripts/eqtl_index_build/query_eqtl_index.R \
    --index data/eQTLcatalogue/eqtl_index_ge_p1e4.tsv.gz \
    --snp   rs429358 --window 500 \
    --ref   /data/h_vmac/HelperScripts/GWAS_Reference_Variant_list.txt
```

---

## Disk usage

| File | Size (approx.) |
|------|---------------|
| Part files (280 × staging) | ~2–5 GB total, deleted after merge |
| `eqtl_index_ge_p1e4.tsv.gz` | ~2.1 GB |
| `eqtl_index_ge_p1e4.tsv.gz.tbi` | ~1 MB |
| **Total permanent** | **~2.1 GB** |
