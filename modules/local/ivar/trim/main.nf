process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_ivar}"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed

    output:
    tuple val(meta), path("*.trimmed.sorted.bam"), path("*.trimmed.sorted.bam.bai"), emit: bam
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar trim -i ${bam} -b ${primer_bed} -p ${prefix}.trimmed -q ${params.ivar_min_bq} -m 30 -e
    samtools sort -@ ${task.cpus} -o ${prefix}.trimmed.sorted.bam ${prefix}.trimmed.bam
    samtools index ${prefix}.trimmed.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep -i version | sed 's/iVar version //' | sed 's/iVar //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trimmed.sorted.bam
    touch ${prefix}.trimmed.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep -i version | sed 's/iVar version //' | sed 's/iVar //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
