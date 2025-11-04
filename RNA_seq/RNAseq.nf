// This is a Nextflow pipeline for RNA sequencing data analysis
nextflow.enable.dsl = 2

// Define default parameters input and output
params.reads = "/home/kiennm/Nextflow_Study/RNA_seq_data/*_{1,2}.fastq.gz"
params.output = "/home/kiennm/Nextflow_Study/results"
params.genome_fasta = "/home/kiennm/Nextflow_Study/genome/genome.fa"
params.gtf = "/home/kiennm/Nextflow_Study/genome/annotation/genes.gtf"
params.star_index = "/home/kiennm/Nextflow_Study/genome/STAR_index"
params.skip_star_index = false  // Set to true if STAR index already exists

// Run Trim_Galore to  trim adapters and low - quality bases from raw RNA- seq reads
process TRIM_GALORE {

tag "$sample_id"
publishDir "${params.output}/trimmed_reads", mode: "copy"

input:
    tuple val(sample_id), path(reads)

output:
    tuple val(sample_id), path("*_val_*.fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: reports

script:
    """
    trim_galore --paired -q 20 --gzip ${reads[0]} ${reads[1]}
    """
}

// Run FastQC on raw reads
process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.output}/QC_reports/raw_fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -t 2 ${reads}
    """
}

// Run FastQC on trimmed reads
process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.output}/QC_reports/trimmed_fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -t 2 ${trimmed_reads}
    """
}

// Run MultiQC to aggregate all QC reports
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

// Generate STAR genome index
process STAR_INDEX {
    tag "Generating STAR Index"
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
         --genomeSAindexNbases 12
    """
}

// Align trimmed reads to the reference genome using STAR
process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.output}/aligned_reads", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_reads)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    path "${sample_id}_Log.final.out", emit: log_final
    path "${sample_id}_Log.out", emit: log_out
    path "${sample_id}_SJ.out.tab", emit: splice_junctions

    script:
    """
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${sample_id}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMunmapped Within \\
         --outSAMattributes Standard \\
         --runThreadN 8 \\
         --quantMode GeneCounts
    """
}

// Index BAM files for downstream analysis
process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.output}/aligned_reads", mode: "copy"

    input:
    tuple val(sample_id), path(aligned_bam)

    output:
    tuple val(sample_id), path(aligned_bam), path("${aligned_bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index ${aligned_bam}
    """
}

// Generate alignment statistics
process SAMTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.output}/alignment_stats", mode: "copy"

    input:
    tuple val(sample_id), path(aligned_bam), path(bam_index)

    output:
    path "${sample_id}_stats.txt", emit: stats
    path "${sample_id}_flagstat.txt", emit: flagstat

    script:
    """
    samtools stats ${aligned_bam} > ${sample_id}_stats.txt
    samtools flagstat ${aligned_bam} > ${sample_id}_flagstat.txt
    """
}

// Quantify gene expression levels using featureCounts
process FEATURE_COUNTS {
    tag "$sample_id"
    publishDir "${params.output}/feature_counts", mode: "copy"

    input:
    tuple val(sample_id), path(aligned_bam), path(bam_index)
    path gtf

    output:
    path "${sample_id}_counts.txt", emit: counts
    path "${sample_id}_counts.txt.summary", emit: summary

    script:
    """
    featureCounts -a ${gtf} \\
                  -o ${sample_id}_counts.txt \\
                  -T 8 \\
                  -p \\
                  -B \\
                  -C \\
                  -Q 10 \\
                  -g gene_id \\
                  -t exon \\
                  ${aligned_bam}
    """
}

// Merge all count files for differential expression analysis
process MERGE_COUNTS {
    publishDir "${params.output}/merged_counts", mode: "copy"

    input:
    path counts

    output:
    path "merged_gene_counts.txt", emit: merged_counts

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    
    count_files = glob.glob("*_counts.txt")
    
    # Read first file to get gene info
    first_df = pd.read_csv(count_files[0], sep='\\t', comment='#')
    merged = first_df[['Geneid']].copy()
    
    # Add counts from all samples
    for count_file in count_files:
        sample_name = count_file.replace('_counts.txt', '')
        df = pd.read_csv(count_file, sep='\\t', comment='#')
        merged[sample_name] = df.iloc[:, -1]  # Last column is the count
    
    merged.to_csv('merged_gene_counts.txt', sep='\\t', index=False)
    """
}

// Main workflow
workflow {
    
    // Print pipeline information
    log.info """
    ================================================================
    RNA-seq Analysis Pipeline
    ================================================================
    Reads            : ${params.reads}
    Genome FASTA     : ${params.genome_fasta}
    GTF Annotation   : ${params.gtf}
    Output Directory : ${params.output}
    Skip STAR Index  : ${params.skip_star_index}
    ================================================================
    """.stripIndent()
    
    // Step 1: Create input channel from paired-end FASTQ files
    read_pairs_ch = channel
        .fromFilePairs(params.reads, size: 2, checkIfExists: true)
    
    // Step 2: Run FastQC on raw reads
    FASTQC_RAW(read_pairs_ch)
    
    // Step 3: Trim adapters and low-quality bases
    TRIM_GALORE(read_pairs_ch)
    
    // Step 4: Run FastQC on trimmed reads
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads)
    
    // Step 5: Generate or use existing STAR index
    if (params.skip_star_index) {
        star_index_ch = channel.fromPath(params.star_index, checkIfExists: true)
    } else {
        genome_ch = channel.fromPath(params.genome_fasta, checkIfExists: true)
        gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
        STAR_INDEX(genome_ch, gtf_ch)
        star_index_ch = STAR_INDEX.out.star_index
    }
    
    // Step 6: Align reads to reference genome
    STAR_ALIGN(TRIM_GALORE.out.trimmed_reads, star_index_ch)
    
    // Step 7: Index BAM files
    SAMTOOLS_INDEX(STAR_ALIGN.out.aligned_bam)
    
    // Step 8: Generate alignment statistics
    SAMTOOLS_STATS(SAMTOOLS_INDEX.out.indexed_bam)
    
    // Step 9: Quantify gene expression
    gtf_for_counts = channel.fromPath(params.gtf, checkIfExists: true)
    FEATURE_COUNTS(SAMTOOLS_INDEX.out.indexed_bam, gtf_for_counts)
    
    // Step 10: Merge count files from all samples
    MERGE_COUNTS(FEATURE_COUNTS.out.counts.collect())
    
    // Step 11: Run MultiQC to aggregate all QC reports
    multiqc_input = FASTQC_RAW.out.fastqc_results
        .mix(FASTQC_TRIMMED.out.fastqc_results)
        .mix(TRIM_GALORE.out.reports)
        .mix(STAR_ALIGN.out.log_final)
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(SAMTOOLS_STATS.out.flagstat)
        .mix(FEATURE_COUNTS.out.summary)
        .collect()
    
    MULTIQC(multiqc_input)
}