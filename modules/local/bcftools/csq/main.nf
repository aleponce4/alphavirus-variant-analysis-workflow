process BCFTOOLS_CSQ {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_bcftools}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta_ref), path(fasta), path(fai)
    tuple val(meta_gff), path(gff)

    output:
    tuple val(meta), path("*.csq.vcf"), emit: vcf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools csq -f ${fasta} -g ${gff} --local-csq ${vcf} -o ${prefix}.csq.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csq.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
