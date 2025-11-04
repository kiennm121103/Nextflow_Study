#!/usr/bin/env bash

# Single-Cell RNA-seq Pipeline - Dependency Installation Script

echo "==================================================================="
echo "Single-Cell RNA-seq Pipeline - Installation"
echo "==================================================================="

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda is not installed."
    echo "Please install Miniconda or Anaconda first:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Creating conda environment: scrna-pipeline"
echo ""

# Create conda environment with all dependencies
conda create -n scrna-pipeline -c bioconda -c conda-forge -y \
    nextflow \
    kb-python \
    fastqc \
    multiqc \
    scanpy \
    python-igraph \
    leidenalg \
    matplotlib \
    seaborn \
    pandas \
    numpy \
    scipy \
    h5py \
    python=3.9

if [ $? -eq 0 ]; then
    echo ""
    echo "==================================================================="
    echo "Installation completed successfully!"
    echo "==================================================================="
    echo ""
    echo "To activate the environment, run:"
    echo "    conda activate scrna-pipeline"
    echo ""
    echo "To test the installation:"
    echo "    kb --version"
    echo "    python -c 'import scanpy; print(scanpy.__version__)'"
    echo ""
    echo "To run the pipeline:"
    echo "    cd scRNAseq"
    echo "    nextflow run scRNAseq.nf"
    echo ""
else
    echo ""
    echo "ERROR: Installation failed!"
    echo "Please check the error messages above."
    exit 1
fi
