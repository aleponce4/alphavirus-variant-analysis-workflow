process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_ivar}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta), path(fai)
    tuple val(meta_gff), path(gff)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools mpileup -aa -A -d 0 -Q 0 -r ${params.viral_contig} --reference ${fasta} ${bam} | \\
        ivar variants -p ${prefix} -q ${params.ivar_min_bq} -t ${params.ivar_min_freq} -m ${params.ivar_min_depth} -r ${fasta} -g ${gff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep iVar | sed 's/iVar version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep iVar | sed 's/iVar version //')
    END_VERSIONS
    """
}
