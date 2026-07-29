process BWA_INDEX {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_bwa}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("bwa"), emit: index
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir bwa
    bwa index -p bwa/${fasta.name} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa
    touch bwa/${fasta.name}.amb
    touch bwa/${fasta.name}.ann
    touch bwa/${fasta.name}.bwt
    touch bwa/${fasta.name}.pac
    touch bwa/${fasta.name}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}
