process REPORT_RUN_SUMMARY {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path results_dir

    output:
    path "*.tsv"        , emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_run_summary.py --results-dir ${results_dir} --output-tsv consolidated_sample_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch consolidated_sample_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
