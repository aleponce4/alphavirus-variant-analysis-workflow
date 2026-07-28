process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_low'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("*.fai"), emit: fai
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
