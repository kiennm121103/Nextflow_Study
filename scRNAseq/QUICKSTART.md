# Quick Start Guide - Single-Cell RNA-seq Pipeline

## 1. Installation

### Using Conda (Recommended)

```bash
# Create environment
conda create -n scrna-pipeline -c bioconda -c conda-forge \
    nextflow \
    kb-python \
    fastqc \
    multiqc \
    scanpy \
    python-igraph \
    leidenalg \
    matplotlib \
    pandas

# Activate
conda activate scrna-pipeline
```

### Using Docker

```bash
# No installation needed, Docker will pull images automatically
nextflow run scRNAseq.nf -profile docker
```

## 2. Prepare Your Data

### Data Structure for 10X Genomics

Your FASTQ files should look like this:

```
scRNA_seq_data/
├── sample1_S1_L001_R1_001.fastq.gz  # Barcode + UMI (R1)
├── sample1_S1_L001_R2_001.fastq.gz  # cDNA reads (R2)
├── sample2_S1_L001_R1_001.fastq.gz
├── sample2_S1_L001_R2_001.fastq.gz
└── ...
```

**Important**:
- R1 = Cell barcodes (16bp) + UMI (10-12bp)
- R2 = cDNA sequence (biological reads)

### Download Reference Genome

```bash
# Example for human genome
mkdir -p genome

# Download genome FASTA
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome/genome.fa

# Download GTF
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf genome/genes.gtf
```

## 3. Configure Pipeline

Edit `nextflow.config`:

```groovy
params {
    reads = "/full/path/to/scRNA_seq_data"
    output = "/full/path/to/results"
    genome_fasta = "/full/path/to/genome.fa"
    gtf = "/full/path/to/genes.gtf"
    
    technology = "10X"
    expect_cells = 3000    // Adjust based on your experiment
}
```

## 4. Run Pipeline

### First Run

```bash
cd /home/kiennm/Nextflow_Study/scRNAseq

nextflow run scRNAseq.nf \
    --reads "scRNA_seq_data" \
    --genome_fasta "genome/genome.fa" \
    --gtf "genome/genes.gtf" \
    --output "results" \
    --expect_cells 3000
```

### Resume After Failure

```bash
nextflow run scRNAseq.nf -resume
```

### With Different Technology

```bash
# For Drop-seq
nextflow run scRNAseq.nf --technology DropSeq

# For inDrop
nextflow run scRNAseq.nf --technology InDrop
```

## 5. Check Results

### View QC Report

```bash
firefox results/QC_reports/multiqc_report.html
```

### Check Key Outputs

```bash
# QC metrics
cat results/scanpy_qc/sample1/sample1_qc_metrics.csv

# Cluster statistics
cat results/scanpy_analysis/sample1/sample1_cluster_stats.csv

# Cell type composition
cat results/cell_types/sample1/sample1_cell_type_composition.csv
```

### View Visualizations

```bash
# UMAP with clusters
eog results/scanpy_analysis/sample1/sample1_umap_clusters.png

# Marker gene heatmap
eog results/marker_genes/sample1/sample1_top_markers_heatmap.png

# Cell types
eog results/cell_types/sample1/sample1_cell_types.png
```

## 6. Explore Results in Python

### Load and Visualize

```python
import scanpy as sc
import matplotlib.pyplot as plt

# Load annotated data
adata = sc.read_h5ad('results/cell_types/sample1/sample1_annotated.h5ad')

# Basic info
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Clusters: {adata.obs['leiden'].unique()}")

# Plot UMAP
sc.pl.umap(adata, color='leiden', legend_loc='on data')

# Plot QC metrics
sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])

# Plot marker genes
sc.pl.umap(adata, color=['CD3D', 'CD79A', 'CD14'], ncols=3)

# Dotplot of top markers
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
```

### Find Differentially Expressed Genes

```python
# Between clusters
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20)

# Get results as dataframe
result = sc.get.rank_genes_groups_df(adata, group='0')
print(result.head(20))
```

### Filter and Subset

```python
# Get specific cluster
cluster_0 = adata[adata.obs['leiden'] == '0'].copy()

# Get cells with high gene counts
high_quality = adata[adata.obs['n_genes_by_counts'] > 1000].copy()

# Get specific genes
genes_of_interest = ['CD3D', 'CD8A', 'CD4', 'IL7R']
subset = adata[:, genes_of_interest].copy()
```

## 7. Advanced Analysis

### Trajectory Analysis

```python
import scanpy as sc

adata = sc.read_h5ad('results/scanpy_analysis/sample1/sample1_clustered.h5ad')

# PAGA trajectory
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, plot=False)

# Diffusion map
sc.tl.diffmap(adata)
sc.pl.diffmap(adata, color='leiden')

# Draw graph layout
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color='leiden', legend_loc='on data')
```

### Compare Conditions

```python
# If you have multiple conditions
import scanpy as sc

# Assuming you have condition info in metadata
adata.obs['condition'] = ['control', 'treatment', ...] 

# Find DE genes between conditions
sc.tl.rank_genes_groups(adata, 'condition', method='wilcoxon')
sc.pl.rank_genes_groups(adata)

# Volcano plot
from matplotlib import pyplot as plt
result = sc.get.rank_genes_groups_df(adata, group='treatment')
plt.figure(figsize=(10, 8))
plt.scatter(result['logfoldchanges'], -np.log10(result['pvals_adj']), alpha=0.5)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-value')
plt.title('Volcano Plot: Treatment vs Control')
```

### Gene Set Enrichment

```python
import gseapy as gp

# Get marker genes for a cluster
markers = sc.get.rank_genes_groups_df(adata, group='0')
gene_list = markers['names'].tolist()

# Run enrichment analysis
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
    organism='Human',
    outdir='enrichment_results'
)

# View results
print(enr.results.head(10))
```

## 8. Export to Other Tools

### Export to Seurat (R)

```python
# In Python
import scanpy as sc

adata = sc.read_h5ad('results/cell_types/sample1/sample1_annotated.h5ad')
adata.write_csvs('seurat_export/', skip_data=False)
```

```R
# In R
library(Seurat)

# Read exported data
counts <- read.csv('seurat_export/X.csv', row.names=1)
metadata <- read.csv('seurat_export/obs.csv', row.names=1)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = t(counts), meta.data = metadata)

# Continue analysis
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

DimPlot(seurat_obj, group.by = "leiden")
```

## 9. Common Issues & Solutions

### Issue: "No files found matching pattern"

```bash
# Solution: Check file naming
ls scRNA_seq_data/*_R1*.fastq.gz
ls scRNA_seq_data/*_R2*.fastq.gz

# Update pattern if needed
nextflow run scRNAseq.nf --reads "scRNA_seq_data/*_L001_R{1,2}_*.fastq.gz"
```

### Issue: Low cell recovery

```bash
# Solution: Adjust expected cells
nextflow run scRNAseq.nf --expect_cells 5000

# Or check kb-count output
cat results/kb_count/sample1/run_info.json
```

### Issue: Too many/few clusters

Edit filtering parameters:
```bash
nextflow run scRNAseq.nf \
    --min_genes 500 \
    --max_genes 5000 \
    --leiden_resolution 0.8
```

### Issue: Memory errors

Edit `nextflow.config`:
```groovy
process {
    withName: 'SCANPY_.*' {
        memory = '32 GB'
    }
}
```

## 10. Quality Control Checklist

✅ **Before Analysis**
- [ ] FASTQ files correctly named (R1/R2)
- [ ] Expected cells parameter is reasonable
- [ ] Reference genome matches sample species

✅ **After QC**
- [ ] Cells recovered: 50-90% of expected
- [ ] Median genes per cell: 500-5000
- [ ] Median UMI per cell: 1000-50000
- [ ] Mitochondrial %: <5% average

✅ **After Clustering**
- [ ] Number of clusters makes biological sense
- [ ] No obvious batch effects in UMAP
- [ ] Marker genes are cluster-specific
- [ ] Cell type annotations are reasonable

## 11. Next Steps

1. **Refine cell type annotations** - Use marker databases
2. **Perform differential expression** - Between conditions/clusters
3. **Run pathway analysis** - KEGG, GO enrichment
4. **Trajectory analysis** - If differentiation is expected
5. **Integration** - Combine multiple samples/batches
6. **Cell-cell communication** - CellPhoneDB, CellChat
7. **Validate** - qPCR, immunofluorescence

## Resources

- **Scanpy Tutorial**: https://scanpy-tutorials.readthedocs.io
- **kb-python Docs**: https://www.kallistobus.tools
- **Best Practices**: Luecken & Theis (2019) Mol Syst Biol
- **Cell Type Markers**: PanglaoDB, CellMarker

## Execution Time

Typical timeline for 3000 cells:
- Quality Control: 10 min
- Index Building: 30-60 min (one-time)
- kb count: 30-60 min
- Scanpy Analysis: 20-30 min
- **Total: ~2-3 hours**

## Getting Help

```bash
# View pipeline help
nextflow run scRNAseq.nf --help

# Check logs
tail -f .nextflow.log

# Validate config
nextflow config scRNAseq.nf
```

For more detailed information, see `README.md`.
