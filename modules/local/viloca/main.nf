process VILOCA {
    tag "$meta.id"
    label 'process_high'

    container "${params.container_shorah}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shorah shotgun -b ${bam} -f ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shorah: \$(shorah -v 2>&1 | sed 's/ShoRAH //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shorah: \$(shorah -v 2>&1 | sed 's/ShoRAH //')
    END_VERSIONS
    """
}
