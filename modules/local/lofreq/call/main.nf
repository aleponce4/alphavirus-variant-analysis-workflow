process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_lofreq}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_bam = bam
    """
    if [ "${params.lofreq_enable_baq}" = "true" ]; then
        echo "WARNING: BAQ enabled for LoFreq (may crash on STAR spliced alignment BAMs)"
        lofreq viterbi -f ${fasta} ${bam} | samtools sort -o ${prefix}.viterbi.bam -
        samtools index ${prefix}.viterbi.bam
        input_bam="${prefix}.viterbi.bam"
    fi

    if [ "${params.lofreq_enable_indelqual}" = "true" ]; then
        echo "WARNING: Indel quality calculation enabled for LoFreq"
        lofreq indelqual --dindel -f ${fasta} \$input_bam -o ${prefix}.indelqual.bam
        samtools index ${prefix}.indelqual.bam
        input_bam="${prefix}.indelqual.bam"
    fi

    lofreq call-parallel --pp-threads ${task.cpus} -f ${fasta} --min-cov ${params.lofreq_min_depth} --min-bq ${params.lofreq_min_bq} --min-alt-bq ${params.lofreq_min_bq} --min-mq ${params.lofreq_min_mq} --sig ${params.lofreq_sig} -o ${prefix}.vcf ${input_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version | grep version | sed 's/version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version | grep version | sed 's/version: //')
    END_VERSIONS
    """
}
