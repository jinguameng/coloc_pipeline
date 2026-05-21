#!/bin/bash
#SBATCH --job-name=coloc_ctrl
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/slurm_controller_%j.log
#SBATCH --error=logs/slurm_controller_%j.log
#SBATCH --account=h_vmac

# Find where the master pipeline is installed by tracing the symlink of Snakefile
# PIPELINE_ROOT=$(dirname $(readlink -f Snakefile))
PIPELINE_ROOT=/data/h_vmac/zhanm32/coloc_pipeline

# Load Anaconda and activate the shared environment
module load miniconda3/23.9.0-0 2>/dev/null || true
source activate "$PIPELINE_ROOT/env"

# Run Snakemake
snakemake \
    --snakefile "$PIPELINE_ROOT/Snakefile" \
    --profile "$PIPELINE_ROOT/snakemake_slurm_profile" \
    --keep-going \
    "$@"