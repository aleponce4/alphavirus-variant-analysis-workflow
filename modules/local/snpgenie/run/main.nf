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
    VCF_INPUT="${vcf}"
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > ${prefix}.uncompressed.vcf
        VCF_INPUT="${prefix}.uncompressed.vcf"
    fi

    perl ${projectDir}/assets/snpgenie/snpgenie.pl --vcfformat=2 --minfreq=0 --snpreport=\${VCF_INPUT} --fastafile=${fasta} --gtffile=${gtf}
    mv population_summary.txt ${prefix}_population_summary.tsv 2>/dev/null || touch ${prefix}_population_summary.tsv
    mv product_results.txt ${prefix}_product_results.tsv 2>/dev/null || touch ${prefix}_product_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpgenie: 1.0.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_population_summary.tsv
    touch ${prefix}_product_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpgenie: 1.0.0
    END_VERSIONS
    """
}
