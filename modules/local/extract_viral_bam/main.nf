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
    TARGET_CONTIG="${contig != null && contig != 'null' ? contig : ''}"
    if [ -z "\$TARGET_CONTIG" ]; then
        TARGET_CONTIG=\$(samtools view -H ${bam} | grep '^@SQ' | head -n1 | cut -f2 | sed 's/SN://')
    fi

    samtools view -b ${bam} "\$TARGET_CONTIG" | samtools sort -o ${prefix}.viral_only.bam -
    samtools index ${prefix}.viral_only.bam

    READ_COUNT=\$(samtools view -c ${prefix}.viral_only.bam)
    echo "Extracted \${READ_COUNT} reads for \${TARGET_CONTIG} in ${prefix}"

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
