process REPORT_RUN_SUMMARY {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("inputs/*")

    output:
    path "consolidated_sample_summary.tsv", emit: summary
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 \$(which generate_run_summary.py || echo ${projectDir}/bin/generate_run_summary.py) --results-dir . --outdir . --dataset ${params.dataset}

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

process REPORT_COVERAGE_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("depth_files/*")

    output:
    path "*.png"        , emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 \$(which generate_coverage_plots.py || echo ${projectDir}/bin/generate_coverage_plots.py) --coverage-dir depth_files --outdir . --dataset ${params.dataset}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.dataset}_coverage_summary.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process REPORT_VARIANT_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("vcf_files/*")

    output:
    path "*.png"                           , emit: plots
    path "variant_frequency_summary_pct.tsv", emit: summary, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 \$(which generate_variant_plots.py || echo ${projectDir}/bin/generate_variant_plots.py) --vcf-dir vcf_files --outdir . --dataset ${params.dataset}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.dataset}_variant_density.png
    touch ${params.dataset}_allele_frequency_spectrum.png
    touch variant_frequency_summary_pct.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
