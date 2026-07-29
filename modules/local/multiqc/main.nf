process MULTIQC {
    label 'process_low'

    container "${params.container_multiqc}"

    input:
    path multiqc_files
    path multiqc_config

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc --force --config ${multiqc_config} .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")
    END_VERSIONS
    """
}
