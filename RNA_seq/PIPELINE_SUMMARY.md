# RNA-seq Pipeline - Project Summary

## 📋 Overview

This comprehensive RNA-seq analysis pipeline performs end-to-end processing of RNA sequencing data, from raw FASTQ files to gene expression quantification and differential expression analysis.

## 🗂️ Project Structure

```
Nextflow_Study/
├── RNAseq.nf                    # Main Nextflow pipeline
├── nextflow.config              # Pipeline configuration
├── README.md                    # Detailed documentation
├── QUICKSTART.md               # Quick start guide
├── install_dependencies.sh     # Dependency installation script
├── deseq2_analysis.R           # Differential expression R script
└── sample_metadata.csv         # Example metadata template
```

## 🔬 Pipeline Workflow

### Analysis Steps

1. **Quality Control (Raw Reads)**
   - Tool: FastQC
   - Purpose: Assess raw read quality

2. **Adapter Trimming & Quality Filtering**
   - Tool: Trim Galore
   - Purpose: Remove adapters and low-quality bases

3. **Quality Control (Trimmed Reads)**
   - Tool: FastQC
   - Purpose: Verify trimming effectiveness

4. **Genome Index Generation** (optional)
   - Tool: STAR
   - Purpose: Create genome index for alignment

5. **Read Alignment**
   - Tool: STAR
   - Purpose: Align reads to reference genome

6. **BAM Processing**
   - Tools: samtools index, samtools stats
   - Purpose: Index and generate alignment statistics

7. **Gene Quantification**
   - Tool: featureCounts
   - Purpose: Count reads per gene

8. **Count Matrix Generation**
   - Tool: Custom Python script
   - Purpose: Merge counts from all samples

9. **QC Report Aggregation**
   - Tool: MultiQC
   - Purpose: Combine all QC metrics into one report

## 📊 Output Files

### Primary Outputs

| File | Description | Use Case |
|------|-------------|----------|
| `merged_gene_counts.txt` | Gene expression matrix | Differential expression analysis |
| `multiqc_report.html` | Comprehensive QC report | Quality assessment |
| `*.bam` files | Aligned reads | Visualization in IGV |
| `*_counts.txt` | Per-sample counts | Individual sample analysis |

### Report Files

| File | Description |
|------|-------------|
| `report.html` | Pipeline execution report |
| `timeline.html` | Execution timeline |
| `trace.txt` | Resource usage trace |
| `dag.svg` | Pipeline DAG visualization |

## ⚙️ Configuration

### Key Parameters

```groovy
params {
    reads          = "path/to/*_{1,2}.fastq.gz"
    output         = "path/to/results"
    genome_fasta   = "path/to/genome.fa"
    gtf            = "path/to/genes.gtf"
    star_index     = "path/to/STAR_index"
    skip_star_index = false
}
```

### Resource Allocation

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| TRIM_GALORE | 4 | 8 GB | 2h |
| STAR_INDEX | 8 | 32 GB | 4h |
| STAR_ALIGN | 8 | 32 GB | 4h |
| FEATURE_COUNTS | 8 | 8 GB | 2h |

## 🚀 Quick Start Commands

### Install Dependencies
```bash
chmod +x install_dependencies.sh
./install_dependencies.sh
conda activate rnaseq-pipeline
```

### Run Pipeline
```bash
# First run (generate STAR index)
nextflow run RNAseq.nf \
    --reads "data/*_{1,2}.fastq.gz" \
    --genome_fasta "genome/genome.fa" \
    --gtf "genome/genes.gtf" \
    --output "results"

# Subsequent runs (use existing index)
nextflow run RNAseq.nf \
    --skip_star_index true \
    --star_index "results/STAR_index/STAR_index" \
    -resume
```

### Differential Expression Analysis
```bash
Rscript deseq2_analysis.R \
    results/merged_counts/merged_gene_counts.txt \
    sample_metadata.csv \
    DESeq2_results
```

## 📦 Software Dependencies

### Core Tools
- Nextflow (>=21.04.0)
- Java (>=11)
- Python 3 with pandas

### Bioinformatics Tools
- Trim Galore (v0.6.x)
- FastQC (v0.11.x)
- MultiQC (v1.x)
- STAR (v2.7.x)
- samtools (v1.x)
- Subread/featureCounts (v2.x)

### R Packages (for DE analysis)
- DESeq2
- ggplot2
- pheatmap
- RColorBrewer

## 💾 Resource Requirements

### Minimum System
- **CPU**: 8 cores
- **RAM**: 32 GB
- **Storage**: 100 GB free

### Recommended System
- **CPU**: 16+ cores
- **RAM**: 64 GB
- **Storage**: 500 GB+ (SSD)

### Expected Execution Time
- Small dataset (4 samples): ~2-3 hours
- Medium dataset (12 samples): ~6-8 hours
- Large dataset (24+ samples): ~12-16 hours

*Note: Times exclude STAR index generation (30-60 min one-time)*

## 🔧 Pipeline Features

### ✅ Implemented
- [x] Paired-end read support
- [x] Quality control (raw and trimmed)
- [x] Adapter trimming
- [x] Genome alignment with STAR
- [x] Gene quantification
- [x] Count matrix generation
- [x] Comprehensive QC reporting
- [x] Resume capability
- [x] Resource optimization
- [x] Multiple execution profiles (local, conda, docker)

### 🔮 Potential Extensions
- [ ] Single-end read support
- [ ] Transcript-level quantification (Salmon/Kallisto)
- [ ] RNA quality metrics (RSeQC)
- [ ] Integrated differential expression
- [ ] Gene set enrichment analysis
- [ ] Batch effect correction
- [ ] Alternative splicing analysis

## 📈 Downstream Analysis Workflow

```
merged_gene_counts.txt
         ↓
    DESeq2/edgeR
         ↓
Differential Expression Results
         ↓
    GO Enrichment
         ↓
   KEGG Pathways
         ↓
    Visualization
```

## 🐛 Troubleshooting

### Common Issues & Solutions

1. **Out of Memory (STAR)**
   - Increase memory in config: `memory = '64 GB'`
   - Reduce `genomeSAindexNbases` for small genomes

2. **Files Not Found**
   - Use absolute paths
   - Check glob pattern matches filenames
   - Verify `checkIfExists: true` paths

3. **STAR Alignment Fails**
   - Check genome/GTF compatibility
   - Ensure sufficient disk space
   - Verify read file quality

4. **Low Mapping Rate**
   - Check adapter contamination
   - Verify correct reference genome
   - Review FastQC reports

## 📚 References

### Tool Citations
1. **Nextflow**: Di Tommaso et al. (2017) Nat Biotech 35:316
2. **STAR**: Dobin et al. (2013) Bioinformatics 29:15
3. **featureCounts**: Liao et al. (2014) Bioinformatics 30:923
4. **DESeq2**: Love et al. (2014) Genome Biol 15:550
5. **MultiQC**: Ewels et al. (2016) Bioinformatics 32:3047

### Useful Resources
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core/rnaseq](https://nf-co.re/rnaseq)
- [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [DESeq2 Vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## 📝 Best Practices

### Data Organization
- Use consistent naming conventions
- Keep raw data separate from results
- Document sample metadata thoroughly
- Backup critical data

### Pipeline Execution
- Test with small dataset first
- Use `-resume` for interrupted runs
- Monitor resource usage
- Check MultiQC report after completion

### Quality Control
- Review FastQC reports for all samples
- Check alignment rates (target >70%)
- Examine read distribution across features
- Validate with biological replicates

## 🤝 Contributing

To extend or modify this pipeline:
1. Fork/clone the repository
2. Create a new branch for features
3. Test changes thoroughly
4. Document new parameters/processes
5. Update README and QUICKSTART guides

## 📄 License

This pipeline is provided under the MIT License.

## 👤 Contact

For questions, issues, or contributions:
- Review documentation in README.md
- Check QUICKSTART.md for common problems
- Examine `.nextflow.log` for debugging

---

**Last Updated**: November 4, 2025
**Version**: 1.0.0
**Nextflow Version**: >=21.04.0
