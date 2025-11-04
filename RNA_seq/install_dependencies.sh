#!/usr/bin/env bash

# RNA-seq Pipeline Installation Script
# This script will install all required dependencies using conda

echo "==================================================================="
echo "RNA-seq Pipeline - Dependency Installation"
echo "==================================================================="

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda is not installed."
    echo "Please install Miniconda or Anaconda first:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Creating conda environment: rnaseq-pipeline"
echo ""

# Create conda environment with all dependencies
conda create -n rnaseq-pipeline -c bioconda -c conda-forge -y \
    nextflow \
    trim-galore \
    fastqc \
    multiqc \
    star \
    samtools \
    subread \
    python=3.9 \
    pandas

if [ $? -eq 0 ]; then
    echo ""
    echo "==================================================================="
    echo "Installation completed successfully!"
    echo "==================================================================="
    echo ""
    echo "To activate the environment, run:"
    echo "    conda activate rnaseq-pipeline"
    echo ""
    echo "To test the pipeline, run:"
    echo "    nextflow run RNAseq.nf --help"
    echo ""
else
    echo ""
    echo "ERROR: Installation failed!"
    echo "Please check the error messages above."
    exit 1
fi
