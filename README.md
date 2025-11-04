# RNA-seq Analysis Pipelines

Comprehensive Nextflow pipelines for bulk and single-cell RNA sequencing data analysis.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

This repository contains two production-ready Nextflow pipelines:

### 1. Bulk RNA-seq Pipeline
Complete analysis workflow from raw FASTQ files to differential expression:
- Quality control (FastQC, MultiQC)
- Adapter trimming (Trim Galore)
- Read alignment (STAR)
- Gene quantification (featureCounts)
- Count matrix generation
- Differential expression (DESeq2)

### 2. Single-Cell RNA-seq Pipeline  
Comprehensive scRNA-seq analysis from raw reads to cell type annotation:
- Quality control
- Pseudoalignment (kallisto|bustools)
- Cell/gene filtering (Scanpy)
- Dimensionality reduction (PCA, UMAP)
- Clustering (Leiden algorithm)
- Marker gene discovery
- Cell type annotation

## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥21.04.0)
- [UV](https://github.com/astral-sh/uv) or Conda/Mamba
- Java 11+

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/rna-seq-pipelines.git
cd rna-seq-pipelines

# Option 1: Using UV (recommended - fast!)
curl -LsSf https://astral.sh/uv/install.sh | sh

# For single-cell pipeline
cd scRNAseq
uv sync

# Option 2: Using Conda
cd scRNAseq
./install_dependencies.sh
```

### Run Pipelines

```bash
# Bulk RNA-seq
cd RNA_seq
nextflow run RNAseq.nf --reads "data/*_{1,2}.fastq.gz" --genome_fasta genome.fa --gtf genes.gtf

# Single-cell RNA-seq
cd scRNAseq
nextflow run scRNAseq.nf --reads "data" --genome_fasta genome.fa --gtf genes.gtf --expect_cells 3000
```

## 📁 Repository Structure

```
.
├── RNA_seq/                    # Bulk RNA-seq pipeline
│   ├── RNAseq.nf              # Main pipeline
│   ├── nextflow.config         # Configuration
│   ├── README.md              # Documentation
│   └── deseq2_analysis.R      # DE analysis
│
├── scRNAseq/                   # Single-cell pipeline
│   ├── scRNAseq.nf            # Main pipeline
│   ├── nextflow.config         # Configuration
│   ├── pyproject.toml         # Python dependencies
│   ├── uv.lock                # Lock file
│   ├── README.md              # Documentation
│   ├── QUICKSTART.md          # Quick guide
│   └── analyze_scrna.py       # Analysis script
│
├── UV_GUIDE.md                # UV package manager guide
├── requirements.txt           # Python requirements
└── .gitignore                 # Git ignore rules
```

## 📚 Documentation

- **Bulk RNA-seq**: See [RNA_seq/README.md](RNA_seq/README.md)
- **Single-cell RNA-seq**: See [scRNAseq/README.md](scRNAseq/README.md)
- **Quick Start**: See [scRNAseq/QUICKSTART.md](scRNAseq/QUICKSTART.md)
- **UV Guide**: See [UV_GUIDE.md](UV_GUIDE.md)

## 🔧 Configuration

Both pipelines use `nextflow.config` for configuration. Key parameters:

```groovy
params {
    reads = "/path/to/fastq"
    output = "/path/to/results"
    genome_fasta = "/path/to/genome.fa"
    gtf = "/path/to/genes.gtf"
}
```

## 📊 Outputs

### Bulk RNA-seq
- Trimmed reads
- Aligned BAM files
- Gene count matrix
- MultiQC report
- DESeq2 results

### Single-cell RNA-seq
- Count matrices
- QC plots
- UMAP visualizations
- Marker genes
- Cell type annotations
- AnnData objects (.h5ad)

## 🛠️ Development

### Using UV (Recommended)

```bash
# Install dependencies
uv sync

# Install with optional features
uv sync --extra dev --extra batch --extra trajectory

# Add new package
uv add scanpy

# Update dependencies
uv lock --upgrade
```

### Using Conda

```bash
# Create environment
conda env create -f environment.yml

# Activate
conda activate rnaseq-pipeline
```

## 🧪 Testing

```bash
# Run tests (if available)
pytest tests/

# Test pipeline with small dataset
nextflow run scRNAseq.nf -profile test
```

## 📝 Citation

If you use these pipelines, please cite:

**Tools:**
- Nextflow: Di Tommaso et al. (2017) Nature Biotechnology
- STAR: Dobin et al. (2013) Bioinformatics
- Scanpy: Wolf et al. (2018) Genome Biology
- kallisto|bustools: Melsted et al. (2021) Nature Biotechnology

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

For issues or questions:
- Open an [issue](https://github.com/yourusername/rna-seq-pipelines/issues)
- Check documentation in respective README files
- Review `.nextflow.log` for detailed errors

## 👥 Authors

- Your Name - Initial work

## 🙏 Acknowledgments

- nf-core community for pipeline inspiration
- Scanpy developers
- All tool developers and contributors

## 📈 Roadmap

- [ ] Add automated testing
- [ ] Implement CI/CD pipelines
- [ ] Add Docker containers
- [ ] Support for more technologies
- [ ] Integrated differential expression in single-cell
- [ ] Batch effect correction options
- [ ] Interactive visualization dashboard
