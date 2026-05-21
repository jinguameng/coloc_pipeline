# coloc_pipeline

GWAS × eQTL colocalization pipeline that consumes **externally-run SuSiEx**
cross-ancestry fine-mapping outputs and runs `coloc::coloc.bf_bf()` against
pre-computed eQTL Catalogue SuSiE LBFs. Replaces the legacy `colocpipe` whose
SuSiE-on-GWAMA step was unreliable due to LD reference mismatch.

## Pipeline process

```
   GWAMA / PLINK  ──►  SuSiEx pipeline  ──►  coloc_pipeline (this repo)
   (sumstats)         (cross-ancestry             (coloc.bf_bf,
                       fine-mapping)               ABF, reports)
```

Run the SuSiEx pipeline first (see https://github.com/jinguameng/susiex_pipeline).
Then point `coloc_pipeline` at its output directory. If SuSiEx returns
`FAIL` or `NULL` for a locus, the SuSiE arm is skipped for that locus and
ABF still runs (unless `--method=SUSIE` is enforced).

## Key methodological note

`coloc.bf_bf()` requires an L×N log-Bayes-factor matrix per side. SuSiEx
gives one LBF per (CS, population, SNP). For each SNP and CS we **asssum LBF
across populations**:

```
    LBF[k, j] = Σ_p  LogBF(CS_k, Pop_p, SNP_j)
```

Populations where a SNP is absent contribute SuSiEx's near-zero placeholders
(`~1e-8`), so the summation behaves correctly without explicit NA handling.
This combined LBF is treated as the GWAS-side input to `coloc.bf_bf`; the
eQTL side stays single-ancestry (eQTL Catalogue SuSiE LBFs) — there is no
SuSiEx step on the eQTL side.

## Layout

```
coloc_pipeline/
├── bin/colocpipe              # launcher
├── Snakefile                  # main workflow
├── config/                    # admin defaults (you edit per-analysis copy)
├── templates/                 # what `colocpipe init` scaffolds
├── snakemake_slurm_profile/   # SLURM config
├── scripts/
│   ├── parse_susiex_output.R       # SuSiEx → LBF (eQTL Catalogue format)
│   ├── run_coloc_one_pair.R        # one (locus, dataset, gene) per call
│   ├── aggregate_locus_results.R   # combine per-pair rows
│   ├── generate_locus_report.R     # per-locus PDF
│   ├── generate_summary_report.R   # cross-locus PDF
│   ├── extract_coloc_region_PLINK.sh
│   ├── extract_coloc_region_GWAMA.sh
│   ├── query_1kg_ld.sh             # LD for plots only
│   ├── save_ld_rds.R
│   ├── fetch_eqtl_region.R         # eQTL Catalogue puller
│   └── build_1kg_ancestry_keep.sh
└── data/                      # 1KG keeps + eQTL Catalogue metadata
```

## 🛠️ Admin Installation (One-Time Setup)

Clone this repository onto a shared filesystem where your group has read and execute permissions.

```bash
git clone <your-github-repo-url> coloc_pipeline
cd coloc_pipeline

## install
bash install.sh
export PATH="$PWD/bin:$PATH"

## Smoke-test
bash verify_install.sh

```

## 🚀 User Workflow

Users do not need to copy the entire repository. They simply use the pipeline launcher to scaffold an analysis in their own workspace.

**1. Initialize an Analysis Directory**
```bash
mkdir ~/my_coloc_analysis
cd ~/my_coloc_analysis

## initialize the pipeline
/path/to/coloc_pipeline/bin/colocpipe init .
```

**2. Configure Your Analysis**
Edit the generated files to match your dataset:
* `config/pipeline.yaml`: point susiex_dir at your SuSiEx outputs
* `config/loci.tsv`: one row per (locus, lead_snp, tissues, genes)

**3. Dry-run * Submit**
```bash
## dry-run
/path/to/coloc_pipeline/bin/colocpipe dry-run

## local run for small job
/path/to/coloc_pipeline/bin/colocpipe run -j 4 ## request 4 CPUs

## submit to SLURM cluster

sbatch submit          # SLURM (uses snakemake_slurm_profile/)
```

## Expected SuSiEx layout

For locus `apoe` and phenotype `SPAREAD`, the pipeline looks for:

```
{susiex_dir}/apoe/SuSiEx.SPAREAD.apoe.snp
{susiex_dir}/apoe/SuSiEx.SPAREAD.apoe.cs
{susiex_dir}/apoe/SuSiEx.SPAREAD.apoe.summary
```

This matches the default layout of `jinguameng/susiex_pipeline`. If yours
differs, override `params.sx_dir` / `params.sx_name` in `Snakefile`.

## 📊 Output tree

```
{output_dir}/{phenotype}/
├── loci/
│   └── {locus_name}/
│       ├── gwas_region.txt
│       ├── ld_{ancestry}.RDS
│       ├── susiex.lbf_variable.txt.gz
│       ├── susiex.cs_summary.tsv
│       ├── susiex.status                  # PASS | FAIL | NULL
│       ├── eqtl_data/{dataset}.{gene}.{nominal,lbf,cs}.tsv
│       ├── coloc_results/{dataset}.{gene}.coloc.RDS
│       ├── coloc_summary.tsv
│       └── coloc_report.pdf               # per-locus PDF
├── all_coloc_summary.tsv                  # cross-locus table
└── summary_report.pdf                     # cross-locus PDF
```

## 📚 Citation
If you use this pipeline, please cite the original SuSiEx and Coloc methodologies:
> Yuan K, Longchamps RJ, Pardiñas AF, et al. Fine-mapping across diverse ancestries drives the discovery of putative causal variants underlying human complex traits and diseases. *Nat Genet* 56, 1841–1850 (2024). doi:10.1038/s41588-024-01870-z

> Gibson, J. et al. (2019). A meta-analysis of genome-wide association studies of epigenetic age acceleration. PLoS Genetics, 15(11), e1008104. https://doi.org/10.1371/journal.pgen.1008104
Cited by: 129

> Wallace, C. (2021). A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genetics, 17(9), e1009440. https://doi.org/10.1371/journal.pgen.1009440
Cited by: 576

