process EXTRACT_VIRAL_BAM {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam), path(bai)
    val contig

    output:
    tuple val(meta), path("*.viral_only.bam"), path("*.viral_only.bam.bai"), emit: bam
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -b ${bam} ${contig} | samtools sort -o ${prefix}.viral_only.bam -
    samtools index ${prefix}.viral_only.bam

    READ_COUNT=\$(samtools view -c ${prefix}.viral_only.bam)
    echo "Extracted \${READ_COUNT} reads for ${contig} in ${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.viral_only.bam
    touch ${prefix}.viral_only.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
