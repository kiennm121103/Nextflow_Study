#!/usr/bin/env python3
"""
Comprehensive Single-Cell RNA-seq Analysis Script

This script performs complete scRNA-seq analysis from h5ad files:
- Quality control and filtering
- Normalization and preprocessing
- Dimensionality reduction
- Clustering
- Marker gene identification
- Cell type annotation
- Visualization

Usage:
    python analyze_scrna.py --input counts_filtered --output analysis_results
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=False)

def load_data(input_path):
    """Load 10X data from directory or h5ad file"""
    print(f"Loading data from {input_path}...")
    
    if Path(input_path).suffix == '.h5ad':
        adata = sc.read_h5ad(input_path)
    else:
        adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
    
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
    return adata

def quality_control(adata, min_genes=200, max_genes=2500, max_mito_pct=5, output_dir='qc'):
    """Perform quality control and filtering"""
    print("\nPerforming quality control...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Save QC metrics before filtering
    qc_df = pd.DataFrame({
        'n_genes': adata.obs['n_genes_by_counts'],
        'n_counts': adata.obs['total_counts'],
        'pct_mito': adata.obs['pct_counts_mt']
    })
    qc_df.to_csv(output_dir / 'qc_metrics_before_filtering.csv')
    
    # Plot QC metrics before filtering
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # Violin plots
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes[0], show=False)
    
    # Scatter plots
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1,0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1,1], show=False)
    axes[1,2].hist(adata.obs['n_genes_by_counts'], bins=50)
    axes[1,2].set_xlabel('Number of genes')
    axes[1,2].set_ylabel('Number of cells')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'qc_plots_before_filtering.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Filter cells
    print(f"Cells before filtering: {adata.n_obs}")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mito_pct, :]
    print(f"Cells after filtering: {adata.n_obs}")
    
    # Plot QC metrics after filtering
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)
    plt.tight_layout()
    plt.savefig(output_dir / 'qc_plots_after_filtering.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def preprocess(adata, n_top_genes=2000, output_dir='preprocess'):
    """Normalize and identify highly variable genes"""
    print("\nPreprocessing data...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    
    # Plot HVGs
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(output_dir / 'highly_variable_genes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Identified {sum(adata.var.highly_variable)} highly variable genes")
    
    # Save raw data for later
    adata.raw = adata
    
    return adata

def cluster_and_visualize(adata, n_neighbors=10, n_pcs=40, resolution=0.5, output_dir='clustering'):
    """Perform dimensionality reduction and clustering"""
    print("\nClustering cells...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Regress out unwanted variation
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Plot PCA variance
    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    plt.savefig(output_dir / 'pca_variance_ratio.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # PCA plot
    sc.pl.pca(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], show=False)
    plt.savefig(output_dir / 'pca_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Clustering
    sc.tl.leiden(adata, resolution=resolution)
    
    # Plot UMAP
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    sc.pl.umap(adata, color='leiden', legend_loc='on data', ax=axes[0,0], show=False, title='Clusters')
    sc.pl.umap(adata, color='total_counts', ax=axes[0,1], show=False)
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1,0], show=False)
    sc.pl.umap(adata, color='pct_counts_mt', ax=axes[1,1], show=False)
    plt.tight_layout()
    plt.savefig(output_dir / 'umap_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Cluster statistics
    cluster_stats = adata.obs.groupby('leiden').agg({
        'n_genes_by_counts': ['mean', 'std'],
        'total_counts': ['mean', 'std'],
        'pct_counts_mt': ['mean', 'std']
    })
    cluster_stats['n_cells'] = adata.obs.groupby('leiden').size()
    cluster_stats.to_csv(output_dir / 'cluster_statistics.csv')
    
    print(f"Identified {adata.obs['leiden'].nunique()} clusters")
    
    return adata

def find_markers(adata, output_dir='markers'):
    """Find marker genes for each cluster"""
    print("\nFinding marker genes...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find markers
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Save all marker genes
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    marker_df = pd.DataFrame({
        group + '_' + key: result[key][group]
        for group in groups 
        for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    })
    marker_df.to_csv(output_dir / 'all_marker_genes.csv', index=False)
    
    # Save top markers per cluster
    top_markers = pd.DataFrame(result['names']).head(20)
    top_markers.to_csv(output_dir / 'top20_markers_per_cluster.csv', index=False)
    
    # Plot markers
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    plt.savefig(output_dir / 'marker_genes_ranking.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Heatmap
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, show_gene_labels=True, show=False)
    plt.savefig(output_dir / 'marker_genes_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Dotplot
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig(output_dir / 'marker_genes_dotplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Violin plot
    sc.pl.rank_genes_groups_violin(adata, n_genes=8, show=False)
    plt.savefig(output_dir / 'marker_genes_violin.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq analysis')
    parser.add_argument('--input', required=True, help='Input directory or h5ad file')
    parser.add_argument('--output', default='analysis_results', help='Output directory')
    parser.add_argument('--min-genes', type=int, default=200, help='Minimum genes per cell')
    parser.add_argument('--max-genes', type=int, default=2500, help='Maximum genes per cell')
    parser.add_argument('--max-mito', type=float, default=5, help='Maximum mitochondrial percentage')
    parser.add_argument('--n-top-genes', type=int, default=2000, help='Number of HVGs')
    parser.add_argument('--resolution', type=float, default=0.5, help='Clustering resolution')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    adata = load_data(args.input)
    
    # Quality control
    adata = quality_control(adata, args.min_genes, args.max_genes, args.max_mito, 
                           output_dir / 'qc')
    
    # Preprocessing
    adata = preprocess(adata, args.n_top_genes, output_dir / 'preprocess')
    
    # Clustering
    adata = cluster_and_visualize(adata, resolution=args.resolution, 
                                  output_dir=output_dir / 'clustering')
    
    # Find markers
    adata = find_markers(adata, output_dir / 'markers')
    
    # Save final object
    adata.write(output_dir / 'final_analyzed.h5ad')
    
    print(f"\n{'='*60}")
    print("Analysis complete!")
    print(f"{'='*60}")
    print(f"Final dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
    print(f"Results saved to: {output_dir}")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
