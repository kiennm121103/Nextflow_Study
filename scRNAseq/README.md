# Single-Cell RNA-seq Analysis Pipeline

A comprehensive Nextflow pipeline for single-cell RNA sequencing data analysis from raw reads to cell type annotation.

## Pipeline Overview

This pipeline performs the following steps:

1. **Quality Control** - FastQC on raw reads
2. **Index Building** - Generate kallisto/bustools index
3. **Pseudoalignment & Counting** - kb-python (kallisto|bustools)
4. **Quality Control** - Filter cells and genes with Scanpy
5. **Preprocessing** - Normalization and HVG selection
6. **Dimensionality Reduction** - PCA and UMAP
7. **Clustering** - Leiden clustering algorithm
8. **Marker Gene Discovery** - Identify cluster-specific genes
9. **Cell Type Annotation** - Automated cell type identification
10. **Reporting** - MultiQC aggregated report

## Supported Technologies

- **10X Genomics** (3' v1, v2, v3, 5')
- **Drop-seq**
- **inDrop**
- **CEL-Seq2**
- **Smart-seq2** (with modifications)

## Requirements

### Software Dependencies

- Nextflow (>= 21.04.0)
- Java (>= 11)

### Bioinformatics Tools

#### Main Workflow (kb-python)
- kb-python (kallisto|bustools)
- FastQC
- MultiQC
- Scanpy (Python)
- pandas, numpy, matplotlib

#### Alternative Workflows
- **CellRanger** (10X Genomics proprietary)
- **STARsolo** (STAR aligner for scRNA-seq)

## Installation

### Option 1: Using Conda

```bash
# Create conda environment
conda create -n scrna-pipeline -c bioconda -c conda-forge \
    nextflow \
    kb-python \
    fastqc \
    multiqc \
    scanpy \
    python-igraph \
    leidenalg

# Activate environment
conda activate scrna-pipeline
```

### Option 2: Using Docker

```bash
# Pull necessary containers
docker pull quay.io/biocontainers/kb-python:0.27.3--pyhdfd78af_0
docker pull quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0
```

## Input Data Structure

### For 10X Genomics Data

```
scRNA_seq_data/
├── sample1_S1_L001_R1_001.fastq.gz  # Barcode + UMI
├── sample1_S1_L001_R2_001.fastq.gz  # cDNA reads
├── sample2_S1_L001_R1_001.fastq.gz
├── sample2_S1_L001_R2_001.fastq.gz
└── ...
```

**Important**: 
- R1 = Cell barcode (16bp) + UMI (10-12bp)
- R2 = cDNA sequence reads

## Configuration

Edit `nextflow.config` or use command-line parameters:

```groovy
params {
    // Input directory containing FASTQ files
    reads = "/path/to/scRNA_seq_data"
    
    // Output directory
    output = "/path/to/results"
    
    // Reference files
    genome_fasta = "/path/to/genome.fa"
    gtf = "/path/to/genes.gtf"
    
    // Technology
    technology = "10X"        // 10X, DropSeq, InDrop, CELSeq2
    expect_cells = 3000       // Expected number of cells
    
    // QC thresholds
    min_genes = 200           // Min genes per cell
    max_genes = 2500          // Max genes per cell
    max_mito_pct = 5          // Max mitochondrial %
}
```

## Running the Pipeline

### Basic Run

```bash
nextflow run scRNAseq.nf
```

### With Custom Parameters

```bash
nextflow run scRNAseq.nf \
    --reads "/path/to/scRNA_seq_data" \
    --output "/path/to/results" \
    --genome_fasta "/path/to/genome.fa" \
    --gtf "/path/to/genes.gtf" \
    --technology "10X" \
    --expect_cells 5000
```

### Using Conda Profile

```bash
nextflow run scRNAseq.nf -profile conda
```

### Using Docker Profile

```bash
nextflow run scRNAseq.nf -profile docker
```

### Resume Failed Run

```bash
nextflow run scRNAseq.nf -resume
```

## Output Structure

```
results/
├── QC_reports/
│   ├── fastqc/                      # FastQC reports
│   ├── multiqc_report.html          # Aggregated QC
│   └── multiqc_data/
├── kallisto_index/                  # Transcriptome index
│   ├── transcriptome.idx
│   └── t2g.txt
├── kb_count/
│   └── sample1/
│       ├── counts_filtered/          # Filtered count matrix
│       ├── counts_unfiltered/        # Unfiltered counts
│       ├── run_info.json
│       └── inspect.json
├── scanpy_qc/
│   └── sample1/
│       ├── sample1_filtered.h5ad     # Filtered AnnData object
│       ├── sample1_qc_metrics.csv
│       ├── sample1_qc_violin_before.png
│       ├── sample1_qc_violin_after.png
│       └── sample1_qc_scatter.png
├── scanpy_preprocess/
│   └── sample1/
│       ├── sample1_preprocessed.h5ad
│       └── sample1_hvg_plot.png      # Highly variable genes
├── scanpy_analysis/
│   └── sample1/
│       ├── sample1_clustered.h5ad    # Final clustered data
│       ├── sample1_pca.png
│       ├── sample1_umap.png
│       ├── sample1_umap_clusters.png
│       └── sample1_cluster_stats.csv
├── marker_genes/
│   └── sample1/
│       ├── sample1_marker_genes.csv
│       ├── sample1_top_markers_heatmap.png
│       ├── sample1_marker_dotplot.png
│       └── sample1_top10_markers_per_cluster.csv
├── cell_types/
│   └── sample1/
│       ├── sample1_annotated.h5ad
│       ├── sample1_cell_types.png
│       └── sample1_cell_type_composition.csv
└── reports/
    ├── report.html                   # Pipeline execution report
    ├── timeline.html
    ├── trace.txt
    └── dag.svg
```

## Key Output Files

### Primary Analysis Files

| File | Description | Format |
|------|-------------|--------|
| `*_clustered.h5ad` | Clustered single-cell data | AnnData (HDF5) |
| `*_annotated.h5ad` | Cell type annotated data | AnnData (HDF5) |
| `*_marker_genes.csv` | Marker genes per cluster | CSV |
| `*_umap_clusters.png` | UMAP visualization | PNG |

### QC Metrics

| File | Description |
|------|-------------|
| `*_qc_metrics.csv` | Per-cell QC statistics |
| `*_cluster_stats.csv` | Per-cluster statistics |
| `multiqc_report.html` | Comprehensive QC report |

## Downstream Analysis

### Loading Results in Python

```python
import scanpy as sc

# Load annotated data
adata = sc.read_h5ad('results/cell_types/sample1/sample1_annotated.h5ad')

# View basic info
print(adata)
print(adata.obs.head())  # Cell metadata
print(adata.var.head())  # Gene metadata

# Plot UMAP
sc.pl.umap(adata, color='leiden')

# Plot marker genes
sc.pl.umap(adata, color=['CD3D', 'MS4A1', 'CD14'])

# Differential expression between conditions
sc.tl.rank_genes_groups(adata, 'condition', method='wilcoxon')
sc.pl.rank_genes_groups(adata)
```

### Exporting for Seurat (R)

```python
# Save as Seurat-compatible format
adata.write_csvs('seurat_output/', skip_data=False)

# Or convert directly
import anndata2ri
anndata2ri.activate()
%load_ext rpy2.ipython

%%R -i adata
library(Seurat)
seurat_obj <- as.Seurat(adata)
```

### Integration with Seurat (R)

```R
library(Seurat)
library(SeuratDisk)

# Convert h5ad to Seurat
Convert("sample1_annotated.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat_obj <- LoadH5Seurat("sample1_annotated.h5seurat")

# Continue analysis in Seurat
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
DimPlot(seurat_obj, group.by = "leiden")
```

## Quality Control Guidelines

### Expected Metrics for Good Quality Data

| Metric | Good Range | Warning | Action if Outside |
|--------|-----------|---------|-------------------|
| Genes per cell | 500-5000 | <500 or >7000 | Adjust filtering |
| UMI per cell | 1000-50000 | <1000 or >100000 | Check sequencing depth |
| Mitochondrial % | <5% | 5-20% | Stricter filtering |
| Doublet rate | <5% | >10% | Run doublet detection |

### Common Issues

1. **High mitochondrial content**
   - Cause: Dead/dying cells
   - Solution: Lower `max_mito_pct` threshold

2. **Low genes per cell**
   - Cause: Poor library quality or low sequencing depth
   - Solution: Check library prep protocol

3. **High doublet rate**
   - Cause: Cell concentration too high
   - Solution: Use doublet detection tools (Scrublet, DoubletFinder)

4. **Low cell recovery**
   - Cause: Incorrect `expect_cells` parameter
   - Solution: Adjust based on loading concentration

## Advanced Analysis Options

### 1. Batch Effect Correction

```python
import scanpy as sc

# Load data
adata = sc.read_h5ad('results/scanpy_preprocess/sample1/sample1_preprocessed.h5ad')

# Harmony integration
import scanpy.external as sce
sce.pp.harmony_integrate(adata, 'batch')

# Or use scVI
import scvi
scvi.model.SCVI.setup_anndata(adata, batch_key='batch')
vae = scvi.model.SCVI(adata)
vae.train()
adata.obsm['X_scVI'] = vae.get_latent_representation()
```

### 2. Trajectory Analysis

```python
import scanpy as sc

# PAGA trajectory
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, plot=False)
sc.tl.draw_graph(adata, init_pos='paga')

# Or use scVelo for RNA velocity
import scvelo as scv
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
```

### 3. Cell-Cell Communication

```python
import cellphonedb

# Prepare data
cellphonedb.method.statistical_analysis(
    meta_file_path='metadata.txt',
    counts_file_path='counts.txt',
    output_path='cellphonedb_output'
)
```

## Troubleshooting

### Issue: Low alignment rate

```bash
# Check FASTQ quality
fastqc sample_R2_001.fastq.gz

# Verify barcode chemistry
kb count --report sample_R1.fastq.gz sample_R2.fastq.gz

# Try different chemistry
nextflow run scRNAseq.nf --chemistry SC3Pv3
```

### Issue: Memory errors

```groovy
// Increase memory in nextflow.config
process {
    withName: SCANPY_CLUSTER {
        memory = '32 GB'
    }
}
```

### Issue: Empty count matrix

```bash
# Check expected cells parameter
nextflow run scRNAseq.nf --expect_cells 1000

# Verify input file format
# R1 should contain barcodes, R2 should contain cDNA
```

## Performance Optimization

### Execution Time Estimates

| Step | Time (per sample) | Notes |
|------|-------------------|-------|
| FastQC | 5-10 min | Depends on file size |
| Index Building | 30-60 min | One-time process |
| kb count | 30-90 min | Depends on read depth |
| QC & Filtering | 5-10 min | |
| Clustering | 10-20 min | Depends on cell count |
| Marker Finding | 10-30 min | |
| **Total** | **~2-3 hours** | Excluding index building |

### Resource Requirements

**Minimum:**
- CPU: 8 cores
- RAM: 16 GB
- Storage: 50 GB

**Recommended:**
- CPU: 16+ cores
- RAM: 32-64 GB
- Storage: 200 GB (SSD preferred)

## Alternative Pipelines

### Using CellRanger (10X Official)

```bash
# Uncomment CELLRANGER_WORKFLOW in scRNAseq.nf
# Requires CellRanger installation and license

cellranger count \
    --id=sample1 \
    --transcriptome=/path/to/refdata-gex-GRCh38-2020-A \
    --fastqs=/path/to/fastqs \
    --sample=sample1
```

### Using STARsolo

```bash
# Uncomment STARSOLO_WORKFLOW in scRNAseq.nf
# Provides similar results to CellRanger

STAR --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile genes.gtf
```

## Citations

### Tools Used

1. **kb-python**: Melsted et al. (2021) Nat Biotechnol
2. **Scanpy**: Wolf et al. (2018) Genome Biol
3. **STAR**: Dobin et al. (2013) Bioinformatics
4. **CellRanger**: 10X Genomics
5. **MultiQC**: Ewels et al. (2016) Bioinformatics

### Recommended Reading

- Luecken & Theis (2019) "Current best practices in single-cell RNA-seq analysis"
- Hie et al. (2020) "Computational methods for single-cell RNA sequencing"
- Zappia et al. (2018) "Exploring the single-cell RNA-seq analysis landscape"

## Support

For issues or questions:
- Check `.nextflow.log` for detailed errors
- Review Scanpy documentation: https://scanpy.readthedocs.io
- kb-python docs: https://www.kallistobus.tools

## License

MIT License
