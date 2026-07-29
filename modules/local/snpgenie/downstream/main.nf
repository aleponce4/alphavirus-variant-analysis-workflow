process SNPGENIE_SUMMARIZE {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("inputs/*")

    output:
    path "merged/population_summary_all_samples.tsv", emit: population_summary
    path "merged/product_results_all_samples.tsv"   , emit: product_results
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    summarize_snpgenie.py --input-dir inputs --output-dir merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p merged
    touch merged/population_summary_all_samples.tsv
    touch merged/product_results_all_samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process SELECTION_DELTA {
    label 'process_medium'

    container "${params.container_python_reporting}"

    input:
    path product_summary
    path manifest

    output:
    path "delta_per_sample.tsv"       , emit: delta_per_sample
    path "delta_kruskal_by_gene.tsv"  , emit: kruskal_by_gene
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    analyze_delta_selection.py --input ${product_summary} --manifest ${manifest} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch delta_per_sample.tsv
    touch delta_kruskal_by_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process SELECTION_LIMMA {
    label 'process_medium'

    container "${params.container_bioconductor_limma}"

    input:
    path delta_per_sample
    path manifest

    output:
    path "limma_overall_by_gene_all_thresholds.tsv"  , emit: limma_overall
    path "limma_contrasts_by_gene_all_thresholds.tsv", emit: limma_contrasts
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    analyze_delta_limma.R --input-dir . --output-dir . --manifest ${manifest}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -n1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    """
    touch limma_overall_by_gene_all_thresholds.tsv
    touch limma_contrasts_by_gene_all_thresholds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -n1 | sed 's/R version //')
    END_VERSIONS
    """
}

process SELECTION_BUILD_TABLES {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path delta_kruskal
    path limma_overall
    path limma_contrasts
    path product_summary
    path manifest

    output:
    path "selection_gene_summary.tsv"   , emit: gene_summary
    path "selection_gene_contrasts.tsv" , emit: gene_contrasts
    path "selection_gene_key_table.tsv" , emit: key_table
    path "methods_key_parameters.tsv"   , emit: methods_report
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_selection_tables.py --base . --product-summary ${product_summary} --manifest ${manifest} --report-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch selection_gene_summary.tsv
    touch selection_gene_contrasts.tsv
    touch selection_gene_key_table.tsv
    touch methods_key_parameters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
