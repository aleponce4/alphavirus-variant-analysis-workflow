process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_ivar}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.consensus.fa"), emit: consensus
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools mpileup -aa -A -d 0 -Q 0 --reference ${fasta} ${bam} | \\
        ivar consensus -p ${prefix}.consensus -m ${params.ivar_consensus_min_cov} -t ${params.ivar_consensus_threshold} -q ${params.ivar_min_bq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep iVar | sed 's/iVar version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep iVar | sed 's/iVar version //')
    END_VERSIONS
    """
}
