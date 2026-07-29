process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.stats"), path("*.flagstat"), path("*.idxstats"), emit: stats
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools stats ${bam} > ${prefix}.stats
    samtools flagstat ${bam} > ${prefix}.flagstat
    samtools idxstats ${bam} > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    touch ${prefix}.flagstat
    touch ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
