# =============================================================================
# coloc_pipeline / Snakefile
#
# Orchestrates the post-SuSiEx colocalization pipeline. The SuSiEx fine-
# mapping step is OUT-OF-SCOPE for this Snakefile; users run the SuSiEx
# pipeline (https://github.com/jinguameng/susiex_pipeline or equivalent)
# first, then point this pipeline at the resulting directory.
#
# ─── Pipeline Modes ──────────────────────────────────────────────────────────
#   1. Targeted Mode : Tests specific user-requested genes listed in loci.tsv.
#   2. Discovery Mode: If gene_ids is left blank or "ALL", the pipeline queries
#                      a local eQTL index to discover genes with p < 1e-4,
#                      dynamically spawning colocalization jobs only for hits.
#
# ─── Rule Graph ──────────────────────────────────────────────────────────────
#
#   extract_gwas_region  (per locus)
#       └─→ gwas_region.txt
#
#   compute_ld           (per locus, single ancestry)
#       └─→ ld_{ancestry}.RDS                        (for regional plots)
#
#   parse_susiex         (per locus)
#       └─→ susiex.lbf_variable.txt.gz
#           susiex.cs_summary.tsv
#           susiex.status                            (PASS/FAIL/NULL marker)
#
#   [CHECKPOINT] fetch_eqtl (per locus)
#       └─→ eqtl_data/manifest.tsv
#           eqtl_data/{dataset_id}.{gene_id}.nominal.tsv
#           eqtl_data/{dataset_id}.{gene_id}.lbf.tsv (if method != ABF)
#           eqtl_data/{dataset_id}.{gene_id}.cs.tsv  (if method != ABF)
#
#   ================== DYNAMIC DAG EVALUATION ==================
#   Snakemake pauses here to read manifest.tsv. It then spawns
#   the downstream jobs for every valid dataset-gene pair found.
#   ============================================================
#
#   run_coloc_pair       (per (locus, dataset_id, gene_id) in manifest)
#       └─→ coloc_results/{dataset_id}.{gene_id}.coloc.RDS
#           coloc_results/{dataset_id}.{gene_id}.summary_row.tsv
#
#   aggregate_locus      (per locus)
#       └─→ coloc_summary.tsv
#
#   locus_report         (per locus)
#       └─→ coloc_report.pdf
#
#   summary_report       (1 global job)
#       └─→ summary_report.pdf  + all_coloc_summary.tsv
#
# =============================================================================

import os
import pandas as pd
from snakemake.shell import shell
PIPELINE_DIR = workflow.basedir

# ─── Config & inputs ──────────────────────────────────────────────────────────
configfile: "config/pipeline.yaml"

PHENOTYPE    = config["phenotype"]
SUSIEX_DIR   = config["susiex_dir"]
GWAS_FILE    = config["gwas_file"]
GWAS_TYPE    = config.get("gwas_type", "GWAMA").upper()
TRAIT_TYPE   = config.get("trait_type", "quant")
SDY          = str(config.get("sdy", "")).strip()
S            = str(config.get("s",   "")).strip()
ANCESTRY     = config["ancestry"]
STUDY        = config["study"]
METHOD       = config.get("coloc_method", "ALL").upper()
QUANT        = config.get("quant_method", "ge")
OUTPUT_DIR = config.get("outdir", "output").removeprefix("./")

BIM_FILE        = config.get("bim_file", "")
METADATA_FILE   = config["metadata_file"]
TABIX_PATHS     = config["tabix_paths"]
LBF_GENE_INDEX  = config.get("lbf_gene_index", "")
LBF_CACHE       = config.get("lbf_cache", "")
KEEP_DIR        = config.get("keep_dir", "")     

# Helper: pick the right GWAS extractor based on format
if GWAS_TYPE == "PLINK":
    EXTRACT_GWAS_SH = f"{PIPELINE_DIR}/scripts/extract_coloc_region_PLINK.sh"
elif GWAS_TYPE == "GWAMA":
    EXTRACT_GWAS_SH = f"{PIPELINE_DIR}/scripts/extract_coloc_region_GWAMA.sh"
else:
    raise ValueError(f"Unknown gwas_type: {GWAS_TYPE} — expected PLINK or GWAMA")

# ─── Load loci.tsv ────────────────────────────────────────────────────────────
loci_df = pd.read_csv("config/loci.tsv", sep="\t", dtype=str).fillna("")
required = ["locus_name", "chr", "lead_snp", "window_kb", "tissues", "gene_ids"]
for col in required:
    if col not in loci_df.columns:
        raise ValueError(f"loci.tsv is missing column: {col}")

LOCI = loci_df["locus_name"].tolist()
LOCUS_META = {row.locus_name: row for row in loci_df.itertuples(index=False)}

def get_gene_ids_param(locus):
    """Returns empty string if user left gene_ids blank or used ALL/NA/NONE."""
    raw = str(LOCUS_META[locus].gene_ids).strip().upper()
    if raw in ["", "NA", "ALL", "NONE", "NAN"]:
        return ""
    return str(LOCUS_META[locus].gene_ids).strip()

print(f"[Snakefile] Initialized {len(LOCI)} loci. eQTL genes will be dynamically discovered via Checkpoint.")

# ─── rule all ─────────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_report.pdf", locus=LOCI),
        f"{OUTPUT_DIR}/{PHENOTYPE}/summary_report.pdf",
        f"{OUTPUT_DIR}/{PHENOTYPE}/all_coloc_summary.tsv"

# ─── extract_gwas_region ──────────────────────────────────────────────────────
rule extract_gwas_region:
    input:
        gwas = GWAS_FILE
    output:
        f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/gwas_region.txt"
    params:
        snp       = lambda w: LOCUS_META[w.locus].lead_snp,
        window_kb = lambda w: LOCUS_META[w.locus].window_kb,
        sh        = EXTRACT_GWAS_SH
    shell:
        """
        bash {params.sh} \
            -f {input.gwas} \
            -s {params.snp} \
            -w {params.window_kb} \
            -o {output}
        """

# ─── compute_ld (1KG, for plots only) ─────────────────────────────────────────
rule compute_ld:
    input:
        gwas_region = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/gwas_region.txt"
    output:
        f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/ld_{ANCESTRY}.RDS"
    params:
        ancestry  = ANCESTRY,
        chr       = lambda w: LOCUS_META[w.locus].chr,
        snp       = lambda w: LOCUS_META[w.locus].lead_snp,
        window_kb = lambda w: LOCUS_META[w.locus].window_kb
    shell:
        """
        # Find the BP of the lead SNP from the GWAS region file
        BP=$(awk -v s={params.snp} 'NR>1 && $1==s {{print $3; exit}}' "{input.gwas_region}")
        
        if [ -z "$BP" ]; then
            echo "[ERROR] Could not find BP for {params.snp} in {input.gwas_region}" >&2
            exit 1
        fi
        
        WIN_KB={params.window_kb}
        HALF=$(( WIN_KB * 1000 / 2 ))
        START=$((BP - HALF)); [ "$START" -lt 1 ] && START=1
        END=$((BP + HALF))
        
        # Assign Snakemake's output to a bash variable first
        OUT_FILE="{output}"
        PREFIX="${{OUT_FILE%.RDS}}"
        
        bash {PIPELINE_DIR}/scripts/query_1kg_ld.sh \
            --ancestry {params.ancestry} \
            --chr {params.chr} \
            --start "$START" \
            --end "$END" \
            --out "$PREFIX"
        """
# ─── parse_susiex ─────────────────────────────────────────────────────────────
rule parse_susiex:
    input:
        gwas_region = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/gwas_region.txt"
    output:
        lbf        = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.lbf_variable.txt.gz",
        cs_summary = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.cs_summary.tsv",
        status     = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.status"
    params:
        sx_dir   = lambda w: f"{SUSIEX_DIR}/{w.locus}",
        sx_name  = lambda w: f"SuSiEx.{PHENOTYPE}.{w.locus}",
        bim      = BIM_FILE,
        prefix   = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/susiex"
    shell:
        """
        Rscript {PIPELINE_DIR}/scripts/parse_susiex_output.R \
            --susiex_dir  {params.sx_dir} \
            --susiex_name {params.sx_name} \
            --gwas        {input.gwas_region} \
            --bim         "{params.bim}" \
            --locus_name  {wildcards.locus} \
            --out_prefix  {params.prefix}
        """

# ─── fetch_eqtl (CHECKPOINT) ──────────────────────────────────────────────────
checkpoint fetch_eqtl:
    input:
        gwas_region = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/gwas_region.txt"
    output:
        manifest = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/eqtl_data/manifest.tsv"
    params:
        chr        = lambda w: LOCUS_META[w.locus].chr,
        snp        = lambda w: LOCUS_META[w.locus].lead_snp,
        window_kb  = lambda w: LOCUS_META[w.locus].window_kb,
        tissues    = lambda w: LOCUS_META[w.locus].tissues,
        gene_ids   = lambda w: get_gene_ids_param(w.locus),
        out_dir    = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/eqtl_data",
        metadata   = METADATA_FILE,
        tabix_paths= TABIX_PATHS,
        lbf_gene_index = LBF_GENE_INDEX,
        lbf_cache  = LBF_CACHE,
        study      = STUDY,
        quant      = QUANT,
        method     = METHOD
    shell:
        """
        WIN_KB={params.window_kb}
        HALF=$(( WIN_KB * 1000 / 2 ))
        
        BP=$(awk -v s={params.snp} 'NR>1 && $1==s {{print $3; exit}}' "{input.gwas_region}")
        
        if [ -z "$BP" ]; then
            echo "[ERROR] Could not find BP for {params.snp} in {input.gwas_region}" >&2
            exit 1
        fi
        START=$((BP - HALF)); [ "$START" -lt 1 ] && START=1
        END=$((BP + HALF))
        
        Rscript {PIPELINE_DIR}/scripts/fetch_eqtl_region.R \
            --chr {params.chr} --start "$START" --end "$END" \
            --metadata {params.metadata} \
            --tabix_paths {params.tabix_paths} \
            --gene_ids "{params.gene_ids}" \
            --study "{params.study}" \
            --tissues "{params.tissues}" \
            --quant {params.quant} \
            --method {params.method} \
            --lbf_gene_index "{params.lbf_gene_index}" \
            --lbf_cache "{params.lbf_cache}" \
            --out {params.out_dir}
        """

# ─── run_coloc_pair ───────────────────────────────────────────────────────────
rule run_coloc_pair:
    input:
        gwas        = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/gwas_region.txt",
        susiex_lbf  = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.lbf_variable.txt.gz",
        susiex_st   = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.status",
        manifest    = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/eqtl_data/manifest.tsv"
    output:
        rds = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_results/{{dataset_id}}.{{gene_id}}.coloc.RDS",
        row = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_results/{{dataset_id}}.{{gene_id}}.summary_row.tsv"
    params:
        eqtl_nominal = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/eqtl_data/{w.dataset_id}.{w.gene_id}.nominal.tsv",
        eqtl_lbf     = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/eqtl_data/{w.dataset_id}.{w.gene_id}.lbf.tsv",
        eqtl_cs      = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/eqtl_data/{w.dataset_id}.{w.gene_id}.cs.tsv",
        method       = METHOD,
        trait_type   = TRAIT_TYPE,
        sdy          = SDY,
        s            = S,
        study        = STUDY,
        quant        = QUANT,
        tissue       = lambda w: LOCUS_META[w.locus].tissues,
        out_prefix   = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/coloc_results/{w.dataset_id}.{w.gene_id}"
    shell:
        """
        Rscript {PIPELINE_DIR}/scripts/run_coloc_one_pair.R \
            --gwas          {input.gwas} \
            --eqtl_nominal  {params.eqtl_nominal} \
            --eqtl_lbf      "{params.eqtl_lbf}" \
            --eqtl_cs       "{params.eqtl_cs}" \
            --susiex_lbf    {input.susiex_lbf} \
            --susiex_status {input.susiex_st} \
            --method        {params.method} \
            --trait_type    {params.trait_type} \
            --sdy           "{params.sdy}" \
            --s             "{params.s}" \
            --locus         {wildcards.locus} \
            --dataset_id    {wildcards.dataset_id} \
            --gene_id       {wildcards.gene_id} \
            --study_label   "{params.study}" \
            --tissue        "{params.tissue}" \
            --quant_method  {params.quant} \
            --out           {params.out_prefix}
        """

# ─── aggregate_locus ──────────────────────────────────────────────────────────
def locus_pair_rds(wildcards):
    """Dynamic evaluation: reads the checkpoint manifest to generate targets."""
    locus = wildcards.locus
    manifest_path = checkpoints.fetch_eqtl.get(locus=locus).output.manifest
    try:
        df = pd.read_csv(manifest_path, sep="\t")
        if df.empty:
            return []
        files = []
        for _, row in df.iterrows():
            d = str(row["dataset_id"]).strip()
            g = str(row["gene_id"]).strip()
            if d and g and d != "nan" and g != "nan":
                files.append(f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{locus}/coloc_results/{d}.{g}.summary_row.tsv")
        return files
    except Exception:
        return []

rule aggregate_locus:
    input:
        rows = locus_pair_rds
    output:
        summary = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_summary.tsv"
    params:
        in_dir = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/coloc_results"
    run:
        # If the checkpoint found ZERO significant genes, bypass R and write an empty summary.
        if not input.rows:
            with open(output.summary, "w") as f:
                f.write("locus\tdataset_id\tstudy_label\ttissue\tcell_type\tgene_id\tmethod\tsusie_detail\tn_gwas_cs\tn_eqtl_cs\tn_snps\tPP4_ABF\tPP4_SuSiE\tbest_PP4\tabf_warnings\n")
        else:
            shell(
                "Rscript {PIPELINE_DIR}/scripts/aggregate_locus_results.R "
                "--in_dir {params.in_dir} "
                "--out    {output.summary}"
            )

# ─── locus_report ─────────────────────────────────────────────────────────────
rule locus_report:
    input:
        summary    = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_summary.tsv",
        cs_summary = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.cs_summary.tsv",
        status     = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/susiex.status",
        ld         = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/ld_{ANCESTRY}.RDS"
    output:
        pdf = f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_report.pdf"
    params:
        results_dir = lambda w: f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{w.locus}/coloc_results",
        lead_snp    = lambda w: LOCUS_META[w.locus].lead_snp,
        window_kb   = lambda w: LOCUS_META[w.locus].window_kb,
        phenotype   = PHENOTYPE,
        gwas_file   = GWAS_FILE,
        gwas_type   = GWAS_TYPE,
        trait_type  = TRAIT_TYPE,
        method      = METHOD,
        ancestry    = ANCESTRY,
        study       = STUDY
    shell:
        """
        Rscript {PIPELINE_DIR}/scripts/generate_locus_report.R \
            --results_dir {params.results_dir} \
            --cs_summary  {input.cs_summary} \
            --status      {input.status} \
            --ld          {input.ld} \
            --locus       {wildcards.locus} \
            --lead_snp    {params.lead_snp} \
            --window_kb   {params.window_kb} \
            --phenotype   "{params.phenotype}" \
            --gwas_file   "{params.gwas_file}" \
            --gwas_type   "{params.gwas_type}" \
            --trait_type  {params.trait_type} \
            --method      {params.method} \
            --ancestry    {params.ancestry} \
            --study       "{params.study}" \
            --out         {output.pdf}
        """

# ─── all_coloc_summary + summary_report ──────────────────────────────────────
rule all_coloc_summary:
    input:
        expand(f"{OUTPUT_DIR}/{PHENOTYPE}/loci/{{locus}}/coloc_summary.tsv", locus=LOCI)
    output:
        f"{OUTPUT_DIR}/{PHENOTYPE}/all_coloc_summary.tsv"
    run:
        dfs = [pd.read_csv(f, sep="\t") for f in input]
        out = pd.concat(dfs, ignore_index=True)
        # Prevent pandas crash if every locus was completely empty
        if not out.empty and "best_PP4" in out.columns:
            out = out.sort_values("best_PP4", ascending=False)
        out.to_csv(output[0], sep="\t", index=False)
        print(f"[all_coloc_summary] {len(out)} rows → {output[0]}")

rule summary_report:
    input:
        summary = f"{OUTPUT_DIR}/{PHENOTYPE}/all_coloc_summary.tsv"
    output:
        pdf = f"{OUTPUT_DIR}/{PHENOTYPE}/summary_report.pdf"
    params:
        phenotype = PHENOTYPE,
        gwas_file = GWAS_FILE,
        method    = METHOD,
        ancestry  = ANCESTRY,
        study     = STUDY
    shell:
        """
        Rscript {PIPELINE_DIR}/scripts/generate_summary_report.R \
            --summary   {input.summary} \
            --phenotype "{params.phenotype}" \
            --gwas_file "{params.gwas_file}" \
            --method    {params.method} \
            --ancestry  {params.ancestry} \
            --study     "{params.study}" \
            --out       {output.pdf}
        """