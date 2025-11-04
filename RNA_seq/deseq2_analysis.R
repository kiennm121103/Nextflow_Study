#!/usr/bin/env Rscript

# Differential Expression Analysis using DESeq2
# This script takes the merged count matrix and performs DE analysis

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript deseq2_analysis.R <count_matrix> <sample_metadata>")
}

count_file <- args[1]
metadata_file <- args[2]
output_dir <- ifelse(length(args) >= 3, args[3], "DESeq2_results")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
cat("Reading count matrix...\n")
counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

cat("Reading sample metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1)

# Ensure samples match between counts and metadata
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

cat(sprintf("Analyzing %d genes across %d samples\n", nrow(counts), ncol(counts)))

# Create DESeq2 object
cat("Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d genes after filtering\n", sum(keep)))

# Run DESeq2
cat("Running differential expression analysis...\n")
dds <- DESeq(dds)

# Get results
res <- results(dds, alpha = 0.05)
cat(sprintf("Found %d significant genes (padj < 0.05)\n", sum(res$padj < 0.05, na.rm = TRUE)))

# Write results
cat("Writing results...\n")
write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = file.path(output_dir, "normalized_counts.csv"))

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

# PCA plot
cat("Generating PCA plot...\n")
pdf(file.path(output_dir, "PCA_plot.pdf"), width = 8, height = 6)
plotPCA(vsd, intgroup = "condition") +
    theme_bw() +
    ggtitle("PCA of samples")
dev.off()

# Sample distance heatmap
cat("Generating sample distance heatmap...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pdf(file.path(output_dir, "sample_distances.pdf"), width = 10, height = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

# MA plot
cat("Generating MA plot...\n")
pdf(file.path(output_dir, "MA_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
dev.off()

# Volcano plot
cat("Generating volcano plot...\n")
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              "Significant", "Not significant")
res_df$significant[is.na(res_df$significant)] <- "Not significant"

pdf(file.path(output_dir, "volcano_plot.pdf"), width = 10, height = 8)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey")) +
    theme_bw() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme(legend.position = "top")
dev.off()

# Heatmap of top differentially expressed genes
cat("Generating heatmap of top DE genes...\n")
top_genes <- head(order(res$padj), 50)
pdf(file.path(output_dir, "top_genes_heatmap.pdf"), width = 10, height = 12)
pheatmap(assay(vsd)[top_genes, ],
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         annotation_col = metadata[, "condition", drop = FALSE],
         scale = "row",
         main = "Top 50 Differentially Expressed Genes")
dev.off()

# Summary statistics
cat("\n=================================================================\n")
cat("Differential Expression Analysis Summary\n")
cat("=================================================================\n")
cat(sprintf("Total genes analyzed: %d\n", nrow(res)))
cat(sprintf("Significant genes (padj < 0.05): %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
cat(sprintf("Upregulated genes (log2FC > 1, padj < 0.05): %d\n", 
            sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)))
cat(sprintf("Downregulated genes (log2FC < -1, padj < 0.05): %d\n", 
            sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)))
cat(sprintf("\nResults saved to: %s\n", output_dir))
cat("=================================================================\n")

cat("\nAnalysis complete!\n")
