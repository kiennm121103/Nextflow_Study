# Single-Cell RNA-seq Pipeline - Complete Summary

## 📋 Overview

A comprehensive Nextflow pipeline for single-cell RNA sequencing analysis using modern tools (kb-python/kallisto|bustools and Scanpy). This pipeline handles 10X Genomics, Drop-seq, inDrop, and other droplet-based scRNA-seq technologies.

## 🗂️ Project Structure

```
scRNAseq/
├── scRNAseq.nf                 # Main Nextflow pipeline
├── nextflow.config             # Configuration file
├── README.md                   # Comprehensive documentation
├── QUICKSTART.md              # Quick start guide
├── install_dependencies.sh    # Installation script
├── analyze_scrna.py           # Standalone analysis script
└── sample_sheet.csv           # Sample metadata template
```

## 🔬 Complete Workflow

### Pipeline Steps

```
Raw 10X FASTQ Files
         ↓
    FastQC (QC)
         ↓
Build Kallisto Index
         ↓
kb-python Count (Pseudoalignment)
         ↓
Scanpy QC (Cell/Gene Filtering)
         ↓
Normalization & HVG Selection
         ↓
PCA → UMAP
         ↓
Leiden Clustering
         ↓
Marker Gene Discovery
         ↓
Cell Type Annotation
         ↓
MultiQC Report
```

### Key Processes

| Process | Tool | Output | Purpose |
|---------|------|--------|---------|
| **FASTQC** | FastQC | QC reports | Quality assessment |
| **BUILD_KALLISTO_INDEX** | kb ref | Index + t2g | Create transcriptome index |
| **KB_COUNT** | kb count | Count matrices | Pseudoalignment & quantification |
| **SCANPY_QC** | Scanpy | Filtered h5ad | Remove low-quality cells |
| **SCANPY_PREPROCESS** | Scanpy | Preprocessed h5ad | Normalization, HVG selection |
| **SCANPY_CLUSTER** | Scanpy | Clustered h5ad | Dimensionality reduction, clustering |
| **FIND_MARKERS** | Scanpy | Marker genes | Identify cluster markers |
| **CELL_TYPE_ANNOTATION** | Scanpy | Annotated h5ad | Assign cell types |

## 🎯 Key Features

### ✅ Technology Support
- **10X Genomics** (3' v1/v2/v3, 5' kits)
- **Drop-seq**
- **inDrop**
- **CEL-Seq2**

### ✅ Analysis Capabilities
- Automated quality control
- Cell and gene filtering
- Normalization (CP10K + log)
- Highly variable gene selection
- PCA, UMAP visualization
- Graph-based clustering (Leiden)
- Differential expression
- Marker gene discovery
- Basic cell type annotation

### ✅ Alternative Workflows
- CellRanger workflow (10X official)
- STARsolo workflow (STAR-based)
- Easy switching between methods

## 📊 Output Files

### Primary Outputs

| File | Format | Description |
|------|--------|-------------|
| `*_clustered.h5ad` | AnnData | Final clustered data with UMAP coords |
| `*_annotated.h5ad` | AnnData | Cell type annotated data |
| `*_marker_genes.csv` | CSV | All marker genes per cluster |
| `*_umap_clusters.png` | PNG | UMAP visualization with clusters |
| `multiqc_report.html` | HTML | Comprehensive QC report |

### Directory Structure

```
results/
├── QC_reports/
│   ├── fastqc/
│   └── multiqc_report.html
├── kallisto_index/
│   ├── transcriptome.idx
│   └── t2g.txt
├── kb_count/
│   └── [sample]/
│       ├── counts_filtered/
│       └── counts_unfiltered/
├── scanpy_qc/
│   └── [sample]/
│       ├── *_filtered.h5ad
│       └── *_qc_*.png
├── scanpy_preprocess/
│   └── [sample]/
│       ├── *_preprocessed.h5ad
│       └── *_hvg_plot.png
├── scanpy_analysis/
│   └── [sample]/
│       ├── *_clustered.h5ad
│       ├── *_pca.png
│       ├── *_umap*.png
│       └── *_cluster_stats.csv
├── marker_genes/
│   └── [sample]/
│       ├── *_marker_genes.csv
│       ├── *_top_markers_heatmap.png
│       └── *_marker_dotplot.png
└── cell_types/
    └── [sample]/
        ├── *_annotated.h5ad
        ├── *_cell_types.png
        └── *_cell_type_composition.csv
```

## ⚙️ Configuration

### Key Parameters

```groovy
params {
    // Paths
    reads = "/path/to/scRNA_seq_data"
    output = "/path/to/results"
    genome_fasta = "/path/to/genome.fa"
    gtf = "/path/to/genes.gtf"
    
    // Technology
    technology = "10X"
    expect_cells = 3000
    
    // QC Thresholds
    min_genes = 200          // Remove cells with < 200 genes
    max_genes = 2500         // Remove doublets
    max_mito_pct = 5         // Remove dying cells
    
    // Analysis
    n_top_genes = 2000       // Number of HVGs
    n_neighbors = 10         // KNN graph
    n_pcs = 40               // Principal components
    leiden_resolution = 0.5  // Clustering resolution
}
```

### Resource Requirements

| Process | CPUs | Memory | Time | Notes |
|---------|------|--------|------|-------|
| FASTQC | 2 | 4 GB | 10m | Per sample |
| BUILD_KALLISTO_INDEX | 8 | 16 GB | 1h | One-time |
| KB_COUNT | 8 | 16 GB | 1h | Main bottleneck |
| SCANPY_QC | 4 | 16 GB | 10m | |
| SCANPY_CLUSTER | 4 | 16 GB | 20m | |
| FIND_MARKERS | 4 | 16 GB | 30m | CPU intensive |

## 🚀 Quick Commands

### Installation
```bash
./install_dependencies.sh
conda activate scrna-pipeline
```

### Run Pipeline
```bash
# Basic run
nextflow run scRNAseq.nf

# Custom parameters
nextflow run scRNAseq.nf \
    --reads "data/*_R{1,2}*.fastq.gz" \
    --genome_fasta "genome.fa" \
    --gtf "genes.gtf" \
    --expect_cells 5000

# With Docker
nextflow run scRNAseq.nf -profile docker

# Resume
nextflow run scRNAseq.nf -resume
```

### Standalone Analysis
```bash
python analyze_scrna.py \
    --input results/kb_count/sample1/counts_filtered \
    --output analysis_results \
    --min-genes 200 \
    --max-genes 2500 \
    --resolution 0.5
```

## 📈 Expected Results

### Quality Metrics (Good Data)

| Metric | Good Range | Warning | Poor |
|--------|-----------|---------|------|
| Cells recovered | 50-90% expected | 30-50% | <30% |
| Median genes/cell | 1000-5000 | 500-1000 | <500 |
| Median UMI/cell | 2000-50000 | 1000-2000 | <1000 |
| Mitochondrial % | <5% | 5-10% | >10% |
| Doublet rate | <5% | 5-10% | >10% |

### Analysis Outputs

**For 3000 cells:**
- Number of clusters: 5-15 (tissue dependent)
- Highly variable genes: ~2000
- Marker genes per cluster: 50-200 significant
- Processing time: 2-3 hours total

## 🔧 Advanced Features

### 1. Batch Integration

```python
import scanpy as sc
import scanpy.external as sce

# Load multiple samples
adata = sc.read_h5ad('clustered.h5ad')

# Harmony integration
sce.pp.harmony_integrate(adata, 'batch')

# Re-cluster
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

### 2. Trajectory Analysis

```python
# PAGA for trajectories
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata)

# Diffusion pseudotime
adata.uns['iroot'] = 0  # Set root cell
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
sc.pl.umap(adata, color='dpt_pseudotime')
```

### 3. Cell Type Annotation (Automated)

```python
# Using celltypist
import celltypist
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model)

# Using scanpy scoring
marker_genes = {
    'T cells': ['CD3D', 'CD3E', 'CD8A'],
    'B cells': ['MS4A1', 'CD79A'],
    'Monocytes': ['CD14', 'LYZ']
}
for cell_type, genes in marker_genes.items():
    sc.tl.score_genes(adata, genes, score_name=f'{cell_type}_score')
```

## 🐛 Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Low cell recovery | Wrong expected_cells | Adjust --expect_cells |
| Empty matrix | Wrong R1/R2 order | R1=barcode, R2=cDNA |
| Memory error | Insufficient RAM | Increase memory in config |
| Too many clusters | High resolution | Lower --leiden_resolution |
| No markers found | Poor separation | Check batch effects |

### Diagnostic Commands

```bash
# Check kb-python output
cat results/kb_count/sample1/run_info.json

# Validate FASTQ files
zcat sample_R1.fastq.gz | head -4  # Check barcodes
zcat sample_R2.fastq.gz | head -4  # Check sequences

# Monitor pipeline
tail -f .nextflow.log

# Check resource usage
nextflow log -f status,name,cpus,memory,duration
```

## 📚 Downstream Analysis

### Differential Expression
```python
# Between conditions
sc.tl.rank_genes_groups(adata, 'condition', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20)

# Between clusters
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
```

### Pathway Analysis
```python
import gseapy as gp

# Get markers
markers = sc.get.rank_genes_groups_df(adata, group='0')
gene_list = markers['names'].tolist()

# Run enrichment
enr = gp.enrichr(gene_list=gene_list,
                gene_sets=['KEGG_2021_Human'],
                organism='Human')
```

### Visualization
```python
# Dotplot of markers
sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# Stacked violin
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden')

# Heatmap
sc.pl.heatmap(adata, marker_genes, groupby='leiden')

# Track plot
sc.pl.tracksplot(adata, marker_genes, groupby='leiden')
```

## 📦 Software Versions

**Tested with:**
- Nextflow: 21.04.0+
- kb-python: 0.27.3
- Scanpy: 1.9.3
- Python: 3.9+
- FastQC: 0.11.9
- MultiQC: 1.13

## 🎓 Learning Resources

### Documentation
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [kb-python Documentation](https://www.kallistobus.tools/)
- [Best Practices Paper](https://www.embopress.org/doi/full/10.15252/msb.20188746)

### Example Datasets
- [10X PBMC Datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets)
- [Single Cell Portal](https://singlecell.broadinstitute.org/)
- [CELLxGENE](https://cellxgene.cziscience.com/)

### Video Tutorials
- [Scanpy Workshop](https://www.youtube.com/watch?v=uvyG9yLuNSE)
- [10X Genomics Webinars](https://www.10xgenomics.com/resources/videos)

## 📊 Performance Benchmarks

### Execution Time (3000 cells)

| Step | Time | Bottleneck |
|------|------|-----------|
| FastQC | 5-10 min | I/O |
| Index Building | 30-60 min | CPU (one-time) |
| kb count | 30-90 min | CPU + I/O |
| QC & Filtering | 5 min | Memory |
| Clustering | 10-20 min | CPU |
| Marker Finding | 10-30 min | CPU |
| **Total** | **~2-3 hours** | |

### Scalability

| Cell Count | Memory | Time | Notes |
|-----------|--------|------|-------|
| 3,000 | 16 GB | 2-3h | Typical |
| 10,000 | 32 GB | 4-6h | Large |
| 50,000 | 64 GB | 8-12h | Very large |
| 100,000+ | 128 GB | 12-24h | Use subsampling |

## 🔬 Comparison with Other Tools

| Tool | Speed | Accuracy | Memory | Ease |
|------|-------|----------|--------|------|
| **kb-python** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| CellRanger | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| STARsolo | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ |
| Alevin | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |

## 📄 Citation

If you use this pipeline, please cite:

**Tools:**
1. **Nextflow**: Di Tommaso et al. (2017) Nat Biotechnol
2. **kallisto|bustools**: Melsted et al. (2021) Nat Biotechnol
3. **Scanpy**: Wolf et al. (2018) Genome Biol
4. **Leiden**: Traag et al. (2019) Sci Rep

**Best Practices:**
- Luecken & Theis (2019) "Current best practices in single-cell RNA-seq" Mol Syst Biol

## 📞 Support

For help:
1. Check `README.md` for detailed documentation
2. Review `QUICKSTART.md` for common workflows
3. Examine `.nextflow.log` for errors
4. Visit [Scanpy Discourse](https://discourse.scverse.org/)
5. Check [kb-python Issues](https://github.com/pachterlab/kb_python/issues)

---

**Version**: 1.0.0  
**Last Updated**: November 4, 2025  
**Author**: Pipeline for single-cell RNA-seq analysis  
**License**: MIT
