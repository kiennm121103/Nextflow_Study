# Comprehensive RNA-seq Analysis Pipeline

A Nextflow pipeline for RNA-seq data analysis from raw reads to gene expression counts.

## Pipeline Overview

This pipeline performs the following steps:

1. **Quality Control (Raw)** - FastQC on raw reads
2. **Trimming** - Adapter and quality trimming with Trim Galore
3. **Quality Control (Trimmed)** - FastQC on trimmed reads
4. **Genome Indexing** - Generate STAR index (optional, if not pre-built)
5. **Alignment** - Align reads to reference genome with STAR
6. **BAM Indexing** - Index aligned BAM files with samtools
7. **Alignment Statistics** - Generate alignment stats with samtools
8. **Quantification** - Count reads per gene with featureCounts
9. **Merge Counts** - Combine counts from all samples
10. **MultiQC Report** - Aggregate all QC metrics

## Requirements

### Software Dependencies

- Nextflow (>= 21.04.0)
- Java (>= 11)
- Conda/Mamba OR Docker/Singularity

### Bioinformatics Tools

- Trim Galore
- FastQC
- MultiQC
- STAR aligner
- Samtools
- featureCounts (from Subread package)
- Python 3 (with pandas)

## Installation

### Option 1: Using Conda

```bash
# Create conda environment
conda create -n rnaseq-pipeline -c bioconda -c conda-forge \
    nextflow \
    trim-galore \
    fastqc \
    multiqc \
    star \
    samtools \
    subread \
    python \
    pandas

# Activate environment
conda activate rnaseq-pipeline
```

### Option 2: Using Docker

```bash
# Pull the container
docker pull nfcore/rnaseq:3.12.0
```

## Input Data Structure

Your FASTQ files should follow this naming convention:

```
RNA_seq_data/
├── sample1_1.fastq.gz  # Forward reads
├── sample1_2.fastq.gz  # Reverse reads
├── sample2_1.fastq.gz
├── sample2_2.fastq.gz
└── ...
```

The pipeline uses `*_{1,2}.fastq.gz` pattern to automatically pair reads.

## Configuration

Edit the parameters in `nextflow.config` or `RNAseq.nf`:

```groovy
params {
    // Input files (use glob pattern for paired-end reads)
    reads = "/path/to/RNA_seq_data/*_{1,2}.fastq.gz"
    
    // Output directory
    output = "/path/to/results"
    
    // Reference genome files
    genome_fasta = "/path/to/genome.fa"
    gtf = "/path/to/annotation.gtf"
    star_index = "/path/to/STAR_index"
    
    // Set to true if STAR index already exists
    skip_star_index = false
}
```

## Running the Pipeline

### Basic Run

```bash
nextflow run RNAseq.nf
```

### With Custom Parameters

```bash
nextflow run RNAseq.nf \
    --reads "/path/to/data/*_{1,2}.fastq.gz" \
    --output "/path/to/results" \
    --genome_fasta "/path/to/genome.fa" \
    --gtf "/path/to/genes.gtf"
```

### With Existing STAR Index

```bash
nextflow run RNAseq.nf \
    --star_index "/path/to/existing/STAR_index" \
    --skip_star_index true
```

### Using Conda Profile

```bash
nextflow run RNAseq.nf -profile conda
```

### Using Docker Profile

```bash
nextflow run RNAseq.nf -profile docker
```

### Resume Failed Run

```bash
nextflow run RNAseq.nf -resume
```

## Output Structure

```
results/
├── trimmed_reads/           # Trimmed FASTQ files
├── QC_reports/             
│   ├── raw_fastqc/         # FastQC reports for raw reads
│   ├── trimmed_fastqc/     # FastQC reports for trimmed reads
│   ├── multiqc_report.html # Aggregated QC report
│   └── multiqc_data/       # MultiQC data
├── STAR_index/             # STAR genome index (if generated)
├── aligned_reads/          # BAM files and indices
├── alignment_stats/        # Alignment statistics
├── feature_counts/         # Per-sample count files
├── merged_counts/          # Merged gene count matrix
│   └── merged_gene_counts.txt
└── reports/                # Pipeline execution reports
    ├── report.html         # Execution report
    ├── timeline.html       # Timeline
    ├── trace.txt          # Resource usage
    └── dag.svg            # Pipeline DAG
```

## Key Output Files

1. **merged_gene_counts.txt** - Gene expression matrix for all samples (ready for DESeq2/edgeR)
2. **multiqc_report.html** - Comprehensive QC report
3. **BAM files** - Aligned reads for visualization in IGV
4. **featureCounts summaries** - Mapping statistics per sample

## Downstream Analysis

The merged gene count file can be used for differential expression analysis:

### In R with DESeq2:

```R
library(DESeq2)

# Read count matrix
counts <- read.table("results/merged_counts/merged_gene_counts.txt", 
                     header=TRUE, row.names=1, sep="\t")

# Create sample metadata
coldata <- data.frame(
    sample = colnames(counts),
    condition = c("control", "control", "treatment", "treatment")
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
results <- results(dds)
```

## Troubleshooting

### Out of Memory Errors

Increase memory in `nextflow.config`:

```groovy
process {
    withName: STAR_ALIGN {
        memory = '64 GB'
    }
}
```

### STAR Index Issues

If you have memory limitations, reduce `genomeSAindexNbases`:

```groovy
--genomeSAindexNbases 10  // For small genomes
```

### File Not Found Errors

- Check file paths are absolute
- Verify glob patterns match your file naming
- Ensure `checkIfExists: true` warnings are addressed

## Pipeline Customization

### Adding New Processes

You can add additional analysis steps like:
- Salmon/Kallisto for transcript quantification
- RSeQC for RNA quality metrics
- DESeq2/edgeR for differential expression
- Gene set enrichment analysis

### Modifying Parameters

Edit `nextflow.config` to change:
- CPU/memory allocation
- Quality trimming parameters
- STAR alignment parameters
- featureCounts options

## Performance Tips

1. Use `-resume` to restart from last successful step
2. Pre-build STAR index and set `skip_star_index = true`
3. Use SSD storage for temporary files
4. Allocate sufficient memory for STAR (32-64 GB recommended)
5. Use multiple CPUs/threads where possible

## Citation

If you use this pipeline, please cite the tools:

- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology
- **Trim Galore**: https://github.com/FelixKrueger/TrimGalore
- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: Ewels et al. (2016) Bioinformatics
- **STAR**: Dobin et al. (2013) Bioinformatics
- **Samtools**: Li et al. (2009) Bioinformatics
- **featureCounts**: Liao et al. (2014) Bioinformatics

## License

MIT License

## Contact

For questions or issues, please open an issue on GitHub.
