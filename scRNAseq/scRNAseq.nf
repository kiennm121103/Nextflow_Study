// This is a Nextflow pipeline for single-cell RNA sequencing data analysis
nextflow.enable.dsl = 2

// Define default parameters for input and output
params.reads = "/home/kiennm/Nextflow_Study/scRNA_seq_data"
params.output = "/home/kiennm/Nextflow_Study/scRNA_results"
params.sample_sheet = "/home/kiennm/Nextflow_Study/scRNAseq/sample_sheet.csv"

// Reference genome parameters
params.genome_fasta = "/home/kiennm/Nextflow_Study/genome/genome.fa"
params.gtf = "/home/kiennm/Nextflow_Study/genome/annotation/genes.gtf"
params.transcriptome_index = "/home/kiennm/Nextflow_Study/genome/kallisto_index"

// Technology and parameters
params.chemistry = "auto"  // 10X chemistry version: auto, SC3Pv2, SC3Pv3, SC5P-PE, SC5P-R2
params.expect_cells = 3000  // Expected number of cells
params.technology = "10X"   // Options: 10X, DropSeq, InDrop, CELSeq2

// Quality control thresholds
params.min_genes = 200      // Minimum genes per cell
params.min_cells = 3        // Minimum cells expressing a gene
params.max_genes = 2500     // Maximum genes per cell (filter doublets)
params.max_mito_pct = 5     // Maximum mitochondrial percentage

// Run FastQC to assess the quality of raw single-cell RNA-seq reads
process FASTQC {
    tag "$sample_id"
    publishDir "${params.output}/QC_reports/fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -t 2 ${reads}
    """
}

// Build transcriptome index for pseudoalignment (Kallisto/bustools)
process BUILD_KALLISTO_INDEX {
    tag "Building Kallisto Index"
    publishDir "${params.output}/kallisto_index", mode: "copy"

    input:
    path genome_fasta
    path gtf

    output:
    path "transcriptome.idx", emit: kallisto_index
    path "t2g.txt", emit: t2g

    script:
    """
    # Extract transcriptome sequences from genome
    kb ref \\
        -i transcriptome.idx \\
        -g t2g.txt \\
        -f1 transcriptome.fa \\
        ${genome_fasta} \\
        ${gtf}
    """
}

// Alternative: Build STAR index for single-cell
process BUILD_STAR_INDEX {
    tag "Building STAR Index"
    publishDir "${params.output}/STAR_index", mode: "copy"

    input:
    path genome_fasta
    path gtf

    output:
    path "STAR_index", emit: star_index

    script:
    """
    mkdir -p STAR_index
    STAR --runMode genomeGenerate \\
         --genomeDir STAR_index \\
         --genomeFastaFiles ${genome_fasta} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN 8 \\
         --genomeSAindexNbases 12 \\
         --sjdbOverhang 100
    """
}

// Run kb-python (kallisto|bustools) for 10X data processing
process KB_COUNT {
    tag "$sample_id"
    publishDir "${params.output}/kb_count/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(reads)
    path kallisto_index
    path t2g

    output:
    tuple val(sample_id), path("counts_unfiltered"), emit: counts_unfiltered
    tuple val(sample_id), path("counts_filtered"), emit: counts_filtered
    path "run_info.json", emit: run_info
    path "inspect.json", emit: inspect

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    kb count \\
        -i ${kallisto_index} \\
        -g ${t2g} \\
        -x ${params.technology} \\
        -t 8 \\
        --h5ad \\
        -o . \\
        ${read1} ${read2}
    """
}

// Alternative: Process 10X data with CellRanger
process CELLRANGER_COUNT {
    tag "$sample_id"
    publishDir "${params.output}/cellranger/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(fastq_dir)
    path transcriptome

    output:
    tuple val(sample_id), path("${sample_id}/outs"), emit: cellranger_output
    path "${sample_id}/outs/web_summary.html", emit: web_summary
    path "${sample_id}/outs/metrics_summary.csv", emit: metrics

    script:
    """
    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id} \\
        --chemistry=${params.chemistry} \\
        --expect-cells=${params.expect_cells} \\
        --localcores=8 \\
        --localmem=64
    """
}

// Run STARsolo for single-cell alignment and quantification
process STARSOLO {
    tag "$sample_id"
    publishDir "${params.output}/starsolo/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(reads)
    path star_index
    path whitelist

    output:
    tuple val(sample_id), path("${sample_id}_Solo.out"), emit: solo_output
    path "${sample_id}_Log.final.out", emit: log_final
    path "${sample_id}_Aligned.sortedByCoord.out.bam", emit: aligned_bam

    script:
    def read1 = reads[0]  // Cell barcode + UMI
    def read2 = reads[1]  // cDNA
    """
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${read2} ${read1} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${sample_id}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN 8 \\
         --soloType CB_UMI_Simple \\
         --soloCBwhitelist ${whitelist} \\
         --soloUMIlen 12 \\
         --soloCBlen 16 \\
         --soloFeatures Gene GeneFull \\
         --soloOutFileNames ${sample_id}_Solo.out/ features.tsv barcodes.tsv matrix.mtx
    """
}

// Quality control and filtering using Scanpy
process SCANPY_QC {
    tag "$sample_id"
    publishDir "${params.output}/scanpy_qc/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(counts_dir)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.h5ad"), emit: filtered_adata
    path "${sample_id}_qc_metrics.csv", emit: qc_metrics
    path "*.png", emit: qc_plots

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Read count matrix
    adata = sc.read_10x_mtx('${counts_dir}', var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Save raw QC metrics
    qc_df = pd.DataFrame({
        'n_genes_by_counts': adata.obs['n_genes_by_counts'],
        'total_counts': adata.obs['total_counts'],
        'pct_counts_mt': adata.obs['pct_counts_mt']
    })
    qc_df.to_csv('${sample_id}_qc_metrics.csv')
    
    # QC plots before filtering
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)
    plt.tight_layout()
    plt.savefig('${sample_id}_qc_violin_before.png', dpi=300)
    plt.close()
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=${params.min_genes})
    sc.pp.filter_genes(adata, min_cells=${params.min_cells})
    adata = adata[adata.obs.n_genes_by_counts < ${params.max_genes}, :]
    adata = adata[adata.obs.pct_counts_mt < ${params.max_mito_pct}, :]
    
    # QC plots after filtering
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)
    plt.tight_layout()
    plt.savefig('${sample_id}_qc_violin_after.png', dpi=300)
    plt.close()
    
    # Scatter plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1], show=False)
    plt.tight_layout()
    plt.savefig('${sample_id}_qc_scatter.png', dpi=300)
    plt.close()
    
    # Save filtered data
    adata.write('${sample_id}_filtered.h5ad')
    
    print(f"Cells before filtering: {qc_df.shape[0]}")
    print(f"Cells after filtering: {adata.n_obs}")
    """
}

// Normalization and preprocessing
process SCANPY_PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output}/scanpy_preprocess/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(filtered_adata)

    output:
    tuple val(sample_id), path("${sample_id}_preprocessed.h5ad"), emit: preprocessed_adata
    path "${sample_id}_hvg_plot.png", emit: hvg_plot

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    # Read filtered data
    adata = sc.read_h5ad('${filtered_adata}')
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig('${sample_id}_hvg_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save preprocessed data
    adata.write('${sample_id}_preprocessed.h5ad')
    
    print(f"Total genes: {adata.n_vars}")
    print(f"Highly variable genes: {sum(adata.var.highly_variable)}")
    """
}

// Dimensionality reduction and clustering
process SCANPY_CLUSTER {
    tag "$sample_id"
    publishDir "${params.output}/scanpy_analysis/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(preprocessed_adata)

    output:
    tuple val(sample_id), path("${sample_id}_clustered.h5ad"), emit: clustered_adata
    path "${sample_id}_pca.png", emit: pca_plot
    path "${sample_id}_umap.png", emit: umap_plot
    path "${sample_id}_umap_clusters.png", emit: umap_clusters
    path "${sample_id}_cluster_stats.csv", emit: cluster_stats

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Read preprocessed data
    adata = sc.read_h5ad('${preprocessed_adata}')
    
    # Regress out effects of total counts and mitochondrial genes
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='total_counts', show=False)
    plt.savefig('${sample_id}_pca.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Clustering with Leiden algorithm
    sc.tl.leiden(adata, resolution=0.5)
    
    # Plot UMAP
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], 
               ncols=3, show=False)
    plt.savefig('${sample_id}_umap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot UMAP with clusters
    sc.pl.umap(adata, color='leiden', legend_loc='on data', show=False)
    plt.savefig('${sample_id}_umap_clusters.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Cluster statistics
    cluster_stats = adata.obs.groupby('leiden').agg({
        'n_genes_by_counts': 'mean',
        'total_counts': 'mean',
        'pct_counts_mt': 'mean'
    })
    cluster_stats['n_cells'] = adata.obs.groupby('leiden').size()
    cluster_stats.to_csv('${sample_id}_cluster_stats.csv')
    
    # Save clustered data
    adata.write('${sample_id}_clustered.h5ad')
    
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
    """
}

// Find marker genes for each cluster
process FIND_MARKERS {
    tag "$sample_id"
    publishDir "${params.output}/marker_genes/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(clustered_adata)

    output:
    path "${sample_id}_marker_genes.csv", emit: marker_genes
    path "${sample_id}_top_markers_heatmap.png", emit: marker_heatmap
    path "${sample_id}_marker_dotplot.png", emit: marker_dotplot
    path "${sample_id}_top10_markers_per_cluster.csv", emit: top_markers

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Read clustered data
    adata = sc.read_h5ad('${clustered_adata}')
    
    # Find marker genes using Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Get marker genes
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker_df = pd.DataFrame({
        group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    })
    marker_df.to_csv('${sample_id}_marker_genes.csv', index=False)
    
    # Get top 10 markers per cluster
    top_markers = pd.DataFrame(result['names']).head(10)
    top_markers.to_csv('${sample_id}_top10_markers_per_cluster.csv', index=False)
    
    # Plot heatmap of top markers
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, show_gene_labels=True, show=False)
    plt.savefig('${sample_id}_top_markers_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Dot plot of top markers
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig('${sample_id}_marker_dotplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Marker gene analysis complete")
    """
}

// Cell type annotation (optional - requires marker database)
process CELL_TYPE_ANNOTATION {
    tag "$sample_id"
    publishDir "${params.output}/cell_types/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(clustered_adata)

    output:
    tuple val(sample_id), path("${sample_id}_annotated.h5ad"), emit: annotated_adata
    path "${sample_id}_cell_types.png", emit: cell_type_plot
    path "${sample_id}_cell_type_composition.csv", emit: cell_type_composition

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Read clustered data
    adata = sc.read_h5ad('${clustered_adata}')
    
    # Simple marker-based annotation (you should customize this)
    # Example markers for common cell types
    marker_genes = {
        'T cells': ['CD3D', 'CD3E', 'CD8A', 'CD4'],
        'B cells': ['MS4A1', 'CD79A', 'CD79B'],
        'NK cells': ['GNLY', 'NKG7', 'NCAM1'],
        'Monocytes': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
        'Dendritic cells': ['FCER1A', 'CST3'],
        'Megakaryocytes': ['PPBP']
    }
    
    # Score each cluster for marker genes
    for cell_type, markers in marker_genes.items():
        # Check which markers exist in dataset
        existing_markers = [m for m in markers if m in adata.var_names]
        if existing_markers:
            sc.tl.score_genes(adata, existing_markers, score_name=f'{cell_type}_score')
    
    # Simple annotation based on scores (customize based on your data)
    # This is a placeholder - real annotation should be done carefully
    adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
    
    # Plot with cell type annotations
    sc.pl.umap(adata, color='cell_type', legend_loc='on data', show=False)
    plt.savefig('${sample_id}_cell_types.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Cell type composition
    composition = adata.obs['cell_type'].value_counts()
    composition_df = pd.DataFrame({
        'cell_type': composition.index,
        'n_cells': composition.values,
        'percentage': (composition.values / composition.sum() * 100).round(2)
    })
    composition_df.to_csv('${sample_id}_cell_type_composition.csv', index=False)
    
    # Save annotated data
    adata.write('${sample_id}_annotated.h5ad')
    
    print("Cell type annotation complete")
    """
}

// MultiQC report aggregation
process MULTIQC {
    publishDir "${params.output}/QC_reports", mode: "copy"

    input:
    path "*"

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}



// Main workflow
workflow {
    
    // Print pipeline information
    log.info """
    ================================================================
    Single-Cell RNA-seq Analysis Pipeline
    ================================================================
    Reads Directory   : ${params.reads}
    Output Directory  : ${params.output}
    Technology        : ${params.technology}
    Expected Cells    : ${params.expect_cells}
    Min Genes/Cell    : ${params.min_genes}
    Max Mito %        : ${params.max_mito_pct}
    ================================================================
    """.stripIndent()
    
    // Step 1: Create input channel from FASTQ files
    // For 10X data: sample_S1_L001_R1_001.fastq.gz and sample_S1_L001_R2_001.fastq.gz
    read_pairs_ch = channel
        .fromFilePairs("${params.reads}/*_R{1,2}*.fastq.gz", size: 2, checkIfExists: true)
    
    // Step 2: Quality control on raw reads
    FASTQC(read_pairs_ch)
    
    // Step 3: Build kallisto/bustools index (if not exists)
    genome_ch = channel.fromPath(params.genome_fasta, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
    BUILD_KALLISTO_INDEX(genome_ch, gtf_ch)
    
    // Step 4: Run kb-python for counting (kallisto|bustools workflow)
    KB_COUNT(
        read_pairs_ch,
        BUILD_KALLISTO_INDEX.out.kallisto_index,
        BUILD_KALLISTO_INDEX.out.t2g
    )
    
    // Step 5: Quality control and filtering with Scanpy
    SCANPY_QC(KB_COUNT.out.counts_filtered)
    
    // Step 6: Preprocessing (normalization, HVG selection)
    SCANPY_PREPROCESS(SCANPY_QC.out.filtered_adata)
    
    // Step 7: Dimensionality reduction and clustering
    SCANPY_CLUSTER(SCANPY_PREPROCESS.out.preprocessed_adata)
    
    // Step 8: Find marker genes
    FIND_MARKERS(SCANPY_CLUSTER.out.clustered_adata)
    
    // Step 9: Cell type annotation
    CELL_TYPE_ANNOTATION(SCANPY_CLUSTER.out.clustered_adata)
    
    // Step 10: Aggregate QC reports
    multiqc_input = FASTQC.out.fastqc_results
        .mix(KB_COUNT.out.run_info)
        .collect()
    
    MULTIQC(multiqc_input)
}

// Alternative workflow using CellRanger (commented out)
/*
workflow CELLRANGER_WORKFLOW {
    // Create sample channel
    samples_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample_id, file(row.fastq_dir)) }
    
    // Build reference
    genome_ch = channel.fromPath(params.genome_fasta)
    gtf_ch = channel.fromPath(params.gtf)
    
    // Run CellRanger
    transcriptome_ch = channel.fromPath(params.transcriptome_index)
    CELLRANGER_COUNT(samples_ch, transcriptome_ch)
    
    // Continue with Scanpy analysis
    SCANPY_QC(CELLRANGER_COUNT.out.cellranger_output)
    SCANPY_PREPROCESS(SCANPY_QC.out.filtered_adata)
    SCANPY_CLUSTER(SCANPY_PREPROCESS.out.preprocessed_adata)
    FIND_MARKERS(SCANPY_CLUSTER.out.clustered_adata)
    CELL_TYPE_ANNOTATION(SCANPY_CLUSTER.out.clustered_adata)
}
*/

// Alternative workflow using STARsolo (commented out)
/*
workflow STARSOLO_WORKFLOW {
    // Prepare inputs
    read_pairs_ch = channel.fromFilePairs("${params.reads}/*_R{1,2}*.fastq.gz", size: 2)
    genome_ch = channel.fromPath(params.genome_fasta)
    gtf_ch = channel.fromPath(params.gtf)
    whitelist_ch = channel.fromPath(params.whitelist)
    
    // Build STAR index
    BUILD_STAR_INDEX(genome_ch, gtf_ch)
    
    // Run STARsolo
    STARSOLO(read_pairs_ch, BUILD_STAR_INDEX.out.star_index, whitelist_ch)
    
    // Continue with Scanpy analysis
    SCANPY_QC(STARSOLO.out.solo_output)
    SCANPY_PREPROCESS(SCANPY_QC.out.filtered_adata)
    SCANPY_CLUSTER(SCANPY_PREPROCESS.out.preprocessed_adata)
    FIND_MARKERS(SCANPY_CLUSTER.out.clustered_adata)
    CELL_TYPE_ANNOTATION(SCANPY_CLUSTER.out.clustered_adata)
}
*/
