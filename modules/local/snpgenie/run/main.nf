process SNPGENIE_RUN {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_perl_bioperl}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta_ref), path(fasta), path(fai)
    tuple val(meta_gtf), path(gtf)

    output:
    tuple val(meta), path("*.tsv"), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    perl ${projectDir}/assets/snpgenie/snpgenie.pl --vcfformat=2 --snpreport=${vcf} --fastafile=${fasta} --gtffile=${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpgenie: 1.0.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_snpgenie_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpgenie: 1.0.0
    END_VERSIONS
    """
}
