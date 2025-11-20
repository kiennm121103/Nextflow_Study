#!/usr/bin/env Rscript

## PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE.

## Load required libraries
library(optparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

## Define command line options
option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=3       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}

## Read in counts file
count.table <- read.delim(file=opt$count_file, header=TRUE, row.names=NULL)
rownames(count.table) <- count.table[,opt$id_col]
count.table <- count.table[,opt$count_col:ncol(count.table), drop = FALSE]
colnames(count.table) <- gsub(paste0(pattern="\\.$", replacement =  '', colnames(count.table))

## Run DESeq2
if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive = TRUE)
}
setwd(opt$outdir)

samples.vec <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components <- length(name_components[[1]])

## Plot QC
##' Generate all required information to plot PCA from a DESeq2 ojbject
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca$x)
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

PlotFile <- paste(opt$outprefix, ".plots.pdf", seq="")

pdf(file=PlotFile, onefile=TRUE, width=7, height=7)

## PCA





