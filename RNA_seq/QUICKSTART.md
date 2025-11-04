# Quick Start Guide - RNA-seq Pipeline

## 1. Setup Environment

```bash
# Make the installation script executable
chmod +x install_dependencies.sh

# Install dependencies
./install_dependencies.sh

# Activate the environment
conda activate rnaseq-pipeline
```

## 2. Prepare Your Data

### File Structure
Organize your FASTQ files with paired-end naming:
```
RNA_seq_data/
├── sample1_1.fastq.gz
├── sample1_2.fastq.gz
├── sample2_1.fastq.gz
├── sample2_2.fastq.gz
└── ...
```

### Reference Genome
Download reference genome and annotation:
```bash
# Example for human genome (GRCh38)
mkdir -p genome

# Download genome FASTA
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome/genome.fa

# Download GTF annotation
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf genome/genes.gtf
```

## 3. Configure Pipeline

Edit `nextflow.config` or use command-line parameters:

```bash
# Option 1: Edit nextflow.config
nano nextflow.config

# Option 2: Use command-line parameters (see below)
```

## 4. Run Pipeline

### First Run (Generate STAR Index)
```bash
nextflow run RNAseq.nf \
    --reads "RNA_seq_data/*_{1,2}.fastq.gz" \
    --genome_fasta "genome/genome.fa" \
    --gtf "genome/genes.gtf" \
    --output "results"
```

### Subsequent Runs (Use Existing Index)
```bash
nextflow run RNAseq.nf \
    --reads "RNA_seq_data/*_{1,2}.fastq.gz" \
    --star_index "results/STAR_index/STAR_index" \
    --gtf "genome/genes.gtf" \
    --output "results" \
    --skip_star_index true
```

### Resume Failed Run
```bash
nextflow run RNAseq.nf -resume
```

## 5. Check Results

### Key Output Files
```bash
# View MultiQC report
firefox results/QC_reports/multiqc_report.html

# Check merged counts
head results/merged_counts/merged_gene_counts.txt

# View pipeline report
firefox results/reports/report.html
```

## 6. Differential Expression Analysis

### Prepare Sample Metadata
Edit `sample_metadata.csv`:
```csv
sample_name,condition,replicate,batch
sample1,control,1,batch1
sample2,control,2,batch1
sample3,treatment,1,batch1
sample4,treatment,2,batch1
```

### Run DESeq2 Analysis
```bash
# Install R packages (first time only)
R -e "install.packages('BiocManager'); BiocManager::install(c('DESeq2', 'ggplot2', 'pheatmap', 'RColorBrewer'))"

# Run analysis
Rscript deseq2_analysis.R \
    results/merged_counts/merged_gene_counts.txt \
    sample_metadata.csv \
    DESeq2_results
```

### View DE Results
```bash
# Check significant genes
head DESeq2_results/deseq2_results.csv

# View plots
firefox DESeq2_results/PCA_plot.pdf
firefox DESeq2_results/volcano_plot.pdf
firefox DESeq2_results/top_genes_heatmap.pdf
```

## 7. Troubleshooting

### Common Issues

**Issue**: Out of memory during STAR alignment
```bash
# Solution: Increase memory in nextflow.config
process {
    withName: STAR_ALIGN {
        memory = '64 GB'
    }
}
```

**Issue**: Files not found
```bash
# Solution: Use absolute paths
nextflow run RNAseq.nf \
    --reads "/full/path/to/RNA_seq_data/*_{1,2}.fastq.gz"
```

**Issue**: STAR index generation fails
```bash
# Solution: Reduce genomeSAindexNbases for small genomes
# Edit RNAseq.nf, line with --genomeSAindexNbases
# Change 12 to a smaller value like 10
```

## 8. Test Run

For a quick test with small data:

```bash
# Create test data directory
mkdir -p test_data

# Download small test dataset (E. coli)
# (Use your own small test dataset or subset of real data)

# Run pipeline
nextflow run RNAseq.nf \
    --reads "test_data/*_{1,2}.fastq.gz" \
    --genome_fasta "test_data/ecoli.fa" \
    --gtf "test_data/ecoli.gtf" \
    --output "test_results"
```

## 9. Pipeline Monitoring

Monitor execution in real-time:

```bash
# In another terminal
tail -f .nextflow.log

# Check resource usage
htop
```

## 10. Clean Up

Remove work directory after successful run:

```bash
# Remove intermediate files (saves disk space)
rm -rf work/

# Keep only the work directory for a specific run
nextflow clean -f
```

## Execution Time Estimates

Typical execution times (per sample):
- FastQC: 5-10 minutes
- Trim Galore: 10-20 minutes
- STAR Index: 30-60 minutes (one-time)
- STAR Align: 20-40 minutes
- featureCounts: 5-10 minutes

Total pipeline time: **2-4 hours** for 4-8 samples (after index is built)

## Resource Requirements

Minimum recommended:
- CPU: 8 cores
- RAM: 32 GB
- Disk: 100 GB free space

Optimal:
- CPU: 16+ cores
- RAM: 64 GB
- Disk: 500 GB+ (SSD preferred)

## Next Steps

1. **Quality Control**: Review MultiQC report for any quality issues
2. **Visualization**: Load BAM files into IGV for visual inspection
3. **Differential Expression**: Run DESeq2 analysis
4. **Gene Set Enrichment**: Perform GO/KEGG enrichment analysis
5. **Biological Interpretation**: Validate findings with RT-qPCR

## Help & Support

For detailed information, see `README.md`

For Nextflow help:
```bash
nextflow run RNAseq.nf --help
```

For issues:
- Check `.nextflow.log` for detailed errors
- Review process logs in `work/` directory
- Consult tool documentation (STAR, featureCounts, etc.)
