#!/usr/bin/env bash

# Script to generate uv.lock files for both pipelines

set -e

echo "=================================================================="
echo "Generating uv.lock files for RNA-seq Pipelines"
echo "=================================================================="
echo ""

# Check if UV is installed
if ! command -v uv &> /dev/null; then
    echo "UV is not installed. Installing now..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
fi

echo "UV version: $(uv --version)"
echo ""

# Generate lock file for Bulk RNA-seq pipeline
echo "=================================================================="
echo "1. Bulk RNA-seq Pipeline"
echo "=================================================================="
cd /home/kiennm/Nextflow_Study

if [ -f "pyproject.toml" ]; then
    echo "Creating uv.lock for bulk RNA-seq pipeline..."
    uv lock
    
    if [ -f "uv.lock" ]; then
        echo "✅ uv.lock created successfully!"
        echo "   Location: /home/kiennm/Nextflow_Study/uv.lock"
        echo "   Size: $(du -h uv.lock | cut -f1)"
    else
        echo "❌ Failed to create uv.lock"
    fi
else
    echo "⚠️  pyproject.toml not found, skipping..."
fi

echo ""

# Generate lock file for Single-Cell RNA-seq pipeline
echo "=================================================================="
echo "2. Single-Cell RNA-seq Pipeline"
echo "=================================================================="
cd /home/kiennm/Nextflow_Study/scRNAseq

if [ -f "pyproject.toml" ]; then
    echo "Creating uv.lock for single-cell RNA-seq pipeline..."
    uv lock
    
    if [ -f "uv.lock" ]; then
        echo "✅ uv.lock created successfully!"
        echo "   Location: /home/kiennm/Nextflow_Study/scRNAseq/uv.lock"
        echo "   Size: $(du -h uv.lock | cut -f1)"
    else
        echo "❌ Failed to create uv.lock"
    fi
else
    echo "⚠️  pyproject.toml not found, skipping..."
fi

echo ""
echo "=================================================================="
echo "Summary"
echo "=================================================================="
echo ""
echo "✅ Lock files generated!"
echo ""
echo "Next steps:"
echo ""
echo "1. To install bulk RNA-seq pipeline:"
echo "   cd /home/kiennm/Nextflow_Study"
echo "   uv sync"
echo ""
echo "2. To install single-cell RNA-seq pipeline:"
echo "   cd /home/kiennm/Nextflow_Study/scRNAseq"
echo "   uv sync"
echo ""
echo "3. To install with optional features:"
echo "   uv sync --extra dev --extra batch --extra trajectory"
echo ""
echo "4. To reproduce on another machine:"
echo "   uv sync --frozen"
echo ""
echo "=================================================================="
