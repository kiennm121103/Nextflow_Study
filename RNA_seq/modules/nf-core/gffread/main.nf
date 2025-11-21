process GFFREAD {
    tag "$meta.id"
    label "process_low"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff)
    path fasta

    output:
    tuple  val(meta), path("*.gtf") , emit: gtf, optional: true
    tuple  val(meta), path("*.gff3"), emit: gffread_gff, optional: true
    tuple  val(meta), path("*.fasta"), emit: gffread_fasta, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    scrip:

    """
    gffread \\
        $gff \\
        $fasta_arg \\
        $args_sorted \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --versions 2>&1)
    END_VERSIONS
    """
    
    
    }