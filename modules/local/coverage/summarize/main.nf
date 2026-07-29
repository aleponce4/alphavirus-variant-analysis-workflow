process COVERAGE_SUMMARIZE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(depth)

    output:
    tuple val(meta), path("*.coverage_summary.tsv"), emit: summary
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 \$(which summarize_coverage.py || echo ${projectDir}/bin/summarize_coverage.py) --depth-file ${depth} --sample-id ${meta.id} --output-tsv ${prefix}.coverage_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.coverage_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
