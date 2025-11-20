process TRIM_GALORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::trim-galore=0.6.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz")    , emit: reads
    tuple val(meta), path("*_trimming_report.txt"), emit: log
    path "*.html"                                , emit: html, optional: true
    path "*.zip"                                 , emit: zip, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Check if single-end or paired-end
    if (meta.single_end) {
        """
        trim_galore \\
            --cores $task.cpus \\
            --gzip \\
            --fastqc \\
            $args \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim-galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        trim_galore \\
            --cores $task.cpus \\
            --gzip \\
            --fastqc \\
            --paired \\
            $args \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim-galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}_trimmed.fq.gz
        touch ${prefix}_trimming_report.txt
        touch ${prefix}_fastqc.html
        touch ${prefix}_fastqc.zip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim-galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_1_val_1.fq.gz
        touch ${prefix}_2_val_2.fq.gz
        touch ${prefix}_1.fastq.gz_trimming_report.txt
        touch ${prefix}_2.fastq.gz_trimming_report.txt
        touch ${prefix}_1_fastqc.html
        touch ${prefix}_2_fastqc.html
        touch ${prefix}_1_fastqc.zip
        touch ${prefix}_2_fastqc.zip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim-galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}