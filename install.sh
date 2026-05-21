#!/bin/bash
# install.sh - Centralized installation script for coloc_pipeline
set -e

echo "=== Initializing Environment Setup ==="

# 1. Load Cluster Miniconda Module
if command -v conda &>/dev/null; then
    echo "Conda is already loaded."
else
    echo "Loading miniconda3 module..."
    module load miniconda3/23.9.0-0 2>/dev/null || { echo "ERROR: Failed to load miniconda3 module"; exit 1; }
fi

# Ensure shell is configured for conda usage
eval "$(conda shell.bash hook)"

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ENV_PATH="$REPO_DIR/env"


# 2. Build local prefix environment
echo "Creating/Updating isolated Conda environment at: $ENV_PATH"
if command -v mamba &>/dev/null; then
    # If the directory already exists, we use 'update' with '--prune'
    if [ -d "$ENV_PATH" ]; then
        mamba env update -p "$ENV_PATH" -f "$REPO_DIR/environment.yml" --prune
    else
        mamba env create -p "$ENV_PATH" -f "$REPO_DIR/environment.yml" -y
    fi
else
    echo "Mamba not found, falling back to Conda (this may take longer)..."
    if [ -d "$ENV_PATH" ]; then
        conda env update -p "$ENV_PATH" -f "$REPO_DIR/environment.yml" --prune
    else
        # Standard conda env create does not accept '--prune' or '-y'
        conda env create -p "$ENV_PATH" -f "$REPO_DIR/environment.yml"
    fi
fi

echo "=== Installation & Environment Setup Complete ==="