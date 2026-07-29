process BWA_MEM {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_bwa}"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    tuple val(meta_ref), path(fasta), path(index_dir)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.bwa.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bwa mem -t ${task.cpus} ${index_dir}/${fasta.name} ${fastq_1} ${fastq_2} 2> ${prefix}.bwa.log > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam
    touch ${prefix}.bwa.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}
