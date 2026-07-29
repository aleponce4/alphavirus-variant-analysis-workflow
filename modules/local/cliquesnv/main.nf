process CLIQUESNV {
    tag "$meta.id"
    label 'process_high'

    container "${params.container_cliquesnv}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.fasta"), emit: haplotypes
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : "-Xmx8g"
    """
    cliquesnv ${avail_mem} -m snv-illumina -in ${bam} -outDir . -threads ${task.cpus}
    mv *.fasta ${prefix}.fasta 2>/dev/null || touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cliquesnv: \$(cliquesnv -version 2>&1 | grep CliqueSNV | sed 's/CliqueSNV v//' || echo "2.0.3")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cliquesnv: \$(cliquesnv -version 2>&1 | grep CliqueSNV | sed 's/CliqueSNV v//' || echo "2.0.3")
    END_VERSIONS
    """
}
