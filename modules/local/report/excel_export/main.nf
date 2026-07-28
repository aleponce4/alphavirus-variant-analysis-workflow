process REPORT_EXCEL_EXPORT {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path input_dir

    output:
    path "*.xlsx"       , emit: excel
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_excel_variants.py --input-dir ${input_dir} --output-xlsx consolidated_variants.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch consolidated_variants.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
