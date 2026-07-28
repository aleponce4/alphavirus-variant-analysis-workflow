process REPORT_HAPLOTYPE_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path haplotype_dir

    output:
    path "*.png"        , emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_haplotypes.py --haplotype-dir ${haplotype_dir} --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch haplotype_plot_stub.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
