# LBF Index Build

Pre-builds a lightweight local gene index from eQTL Catalogue SuSiE LBF
(Log Bayes Factor) files. The index records which genes were fine-mapped
(i.e. have an entry in the LBF file) for each dataset. It is used by
`fetch_eqtl_region.R` to decide whether to stream a dataset's LBF file from
EBI FTP: if none of the target genes appear in the index for a given dataset,
the download is skipped entirely.

**Why not store the full LBF files locally?**  
Each LBF file is 150 MB – 2.3 GB. Storing all 280 would require ~300 GB. The
gene index achieves the same skip-if-irrelevant benefit at ~26 MB total.

**What happens when a gene IS in the index?**  
`fetch_eqtl_region.R` streams the full LBF file from EBI FTP via a
`curl | zcat | awk` pipeline, filtering to the requested region and gene IDs
on the fly. Only matching rows enter R memory — the full file is never stored.

---

## Scripts — in execution order

---

### 1. `setup_lbf_tasklist.sh`

**What it does**  
One-time setup. Joins `dataset_metadata_r7.tsv` with `tabix_ftp_paths.tsv`,
filters to datasets that have a non-empty `ftp_lbf_path`, and writes a TSV
task list used by `build_lbf_index.slurm`. Also creates the `lbf_index/`
output directory.

Run once before submitting the array job. Accepts an optional quantification
method argument (default: `ge`).

```bash
# From the colocpipe root directory
bash scripts/lbf_index_build/setup_lbf_tasklist.sh ge
```

**Output**
- `data/eQTLcatalogue/lbf_tasklist_ge.tsv` — 3 columns: `dataset_id`, `lbf_url`, `study_label`
- `data/eQTLcatalogue/lbf_index/` — output directory (empty at this stage)

The script prints the total task count and the exact `sbatch` command for
step 2.

---

### 2. `build_lbf_index.slurm`

**What it does**  
Thin SLURM array wrapper. Reads one row from `lbf_tasklist_ge.tsv` per task
and delegates all processing to `process_lbf_dataset.sh`. Also cleans up any
leftover manual test directories (`tmp_test/`) and stale large `.lbf.tsv.gz`
files from previous pipeline versions.

**Note:** `htslib` is not required — bgzip and tabix are no longer used.

```bash
# N should match the number of rows in the task list
N=$(tail -n +2 data/eQTLcatalogue/lbf_tasklist_ge.tsv | wc -l)
sbatch --array=1-${N}%10 scripts/lbf_index_build/build_lbf_index.slurm
```

The `%10` limits concurrent tasks to 10 to respect EBI FTP rate limits.
Increase to `%20` if no curl errors occur; decrease to `%5` if downloads
fail frequently.

**Test a single dataset manually before submitting:**
```bash
LINE=$(awk 'NR==2' data/eQTLcatalogue/lbf_tasklist_ge.tsv)
DS=$(echo "$LINE" | cut -f1)
URL=$(echo "$LINE" | cut -f2)

KEEP_RAW=1 \
WORK_DIR=data/eQTLcatalogue/lbf_index/tmp_test \
bash scripts/lbf_index_build/process_lbf_dataset.sh \
    "$DS" "$URL" data/eQTLcatalogue/lbf_index
```

---

### 3. `process_lbf_dataset.sh`

**What it does**  
Core processing script for one dataset. Called by `build_lbf_index.slurm`
for each array task, but can also be run directly for testing. Steps:

1. **Download** — fetches `{dataset_id}.lbf_variable.txt.gz` from EBI FTP via
   `curl` with retry logic. Uses a panfs work directory (not `/tmp`) to avoid
   ramdisk overflow on compute nodes.
2. **Inspect header** — identifies column indices for `molecular_trait_id`,
   `chromosome`, and `position` dynamically from the header row (column order
   varies across eQTL Catalogue releases).
3. **Extract gene list** — pipes `zcat | awk | sort -u` to collect all unique
   gene IDs present in the file. Writes a two-column TSV:
   `dataset_id` + `gene_id`.
4. **Delete raw file** — removes the downloaded `.txt.gz` immediately after
   extraction. Only the tiny gene list is kept.
5. **Write sentinel** — touches `{dataset_id}.lbf.genelist.done` to mark
   success so resubmitted arrays skip completed tasks.

**Not normally called directly** — use `build_lbf_index.slurm` for batch
runs. Call directly only for manual testing.

**Optional environment variables:**

| Variable | Default | Description |
|----------|---------|-------------|
| `WORK_DIR` | `lbf_index/tmp_{dataset_id}` | Where to download the raw file |
| `KEEP_RAW` | `0` | Set to `1` to keep raw download for debugging |
| `N_CPUS` | `2` | CPUs for `sort` |

**Outputs (per dataset)**
- `{index_dir}/{dataset_id}.lbf.genelist.tsv` — unique gene IDs (~KB)
- `{index_dir}/{dataset_id}.lbf.genelist.done` — empty sentinel

---

### 4. `check_lbf_index.sh`

**What it does**  
Progress check and resubmission helper. Scans `lbf_index/` for `.done`
sentinel files and reports how many datasets are complete. For any that are
missing, it prints the corresponding SLURM array task IDs and the exact
`sbatch --array=...` command to resubmit only the failed tasks.

Run any time after submitting the array job to monitor progress, and again
after resubmission to confirm all tasks completed.

```bash
bash scripts/lbf_index_build/check_lbf_index.sh ge
```

**Output (example)**
```
============================================================
 LBF index build status (quant=ge)
============================================================
 Total datasets : 280
 Done           : 267
 Remaining      : 13

 Missing datasets:
   QTD000031
   QTD000056
   ...

   sbatch --array=6,11,...%10 scripts/lbf_index_build/build_lbf_index.slurm
```

---

### 5. `merge_lbf_geneindex.sh`

**What it does**  
Run once after all 280 tasks are confirmed complete by `check_lbf_index.sh`.
Concatenates all per-dataset gene lists into a single master index file used
by `fetch_eqtl_region.R` at pipeline runtime.

```bash
bash scripts/lbf_index_build/merge_lbf_geneindex.sh ge
```

**Output**
- `data/eQTLcatalogue/lbf_gene_index.tsv` — header: `dataset_id`, `gene_id`
  (~16 MB, one row per dataset × fine-mapped gene pair)

After this file exists, `colocpipe.sh` picks it up automatically via the
`--lbf_gene_index` argument passed to `fetch_eqtl_region.R`. No further
configuration is needed.

The script also prints a reminder command to delete any stale large
`.lbf.tsv.gz` files if they exist from a previous indexing strategy:
```bash
rm -f data/eQTLcatalogue/lbf_index/*.lbf.tsv.gz \
       data/eQTLcatalogue/lbf_index/*.lbf.tsv.gz.csi
```

---

## Full workflow summary

```bash
# 1. Generate task list (once)
bash scripts/lbf_index_build/setup_lbf_tasklist.sh ge

# 2. Test one dataset manually before committing to all 280
LINE=$(awk 'NR==2' data/eQTLcatalogue/lbf_tasklist_ge.tsv)
DS=$(echo "$LINE" | cut -f1); URL=$(echo "$LINE" | cut -f2)
KEEP_RAW=1 WORK_DIR=data/eQTLcatalogue/lbf_index/tmp_test \
bash scripts/lbf_index_build/process_lbf_dataset.sh \
    "$DS" "$URL" data/eQTLcatalogue/lbf_index

# 3. Submit full array
N=$(tail -n +2 data/eQTLcatalogue/lbf_tasklist_ge.tsv | wc -l)
sbatch --array=1-${N}%10 scripts/lbf_index_build/build_lbf_index.slurm

# 4. Monitor and resubmit failures as needed
bash scripts/lbf_index_build/check_lbf_index.sh ge

# 5. Merge into master index (once all tasks done)
bash scripts/lbf_index_build/merge_lbf_geneindex.sh ge
```

---

## Disk usage

| Files | Count | Size each | Total |
|-------|-------|-----------|-------|
| `*.lbf.genelist.tsv` | 280 | ~10–50 KB | ~10 MB |
| `*.lbf.genelist.done` | 280 | 0 bytes | ~0 |
| `lbf_gene_index.tsv` | 1 | ~16 MB | ~16 MB |
| Raw downloads (temp, deleted) | 1 at a time | 150 MB – 2.3 GB | 0 MB permanent |
| **Total permanent** | | | **~26 MB** |

---

## How `fetch_eqtl_region.R` uses the gene index

At pipeline runtime, for each dataset selected by study/tissue/quant filters:

1. Nominal associations are fetched via tabix from the local eQTL index
   (`eqtl_index_ge_p1e4.tsv.gz`) — this is fast.
2. Genes found in the nominal data are looked up in `lbf_gene_index.tsv`
   for that dataset.
3. **If no match** → LBF fetch skipped entirely. Log message:
   `lbf: [SKIP] gene not fine-mapped in this dataset — FTP download skipped`
4. **If match** → LBF file streamed from EBI FTP via `curl | zcat | awk`,
   filtered to the GWAS region and matching genes on the fly. Only matching
   rows enter R memory.
