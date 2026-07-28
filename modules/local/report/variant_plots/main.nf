process REPORT_VARIANT_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path vcf_dir

    output:
    path "*.png"        , emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_variants.py --vcf-dir ${vcf_dir} --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch variant_plot_stub.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
