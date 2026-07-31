process HAPLOTYPE_TABLES {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("cliquesnv/*")
    path fasta
    path gff
    path samplesheet

    output:
    path "haplotype_frequency_by_sample.csv", emit: frequency_table
    path "haplotype_summary.csv"            , emit: summary_table
    path "haplotype_sequences.fasta"        , emit: sequences_fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samplesheet_arg = samplesheet ? "--samplesheet ${samplesheet}" : ""
    """
    python3 \$(which build_haplotype_tables.py || echo ${projectDir}/bin/build_haplotype_tables.py) \\
        --cliquesnv-input cliquesnv \\
        --reference-fasta ${fasta} \\
        ${samplesheet_arg} \\
        --out-dir . \\
        --tier-report ${params.haplotype_report_min_freq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch haplotype_frequency_by_sample.csv
    touch haplotype_summary.csv
    touch haplotype_sequences.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process HAPLOTYPE_LINKAGE {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path ("viloca/*")
    path samplesheet

    output:
    path "linked_mutations_long.csv"     , emit: long_table
    path "linked_mutations_recurrent.csv", emit: recurrent_table
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samplesheet_arg = samplesheet ? "--samplesheet ${samplesheet}" : ""
    """
    python3 \$(which summarize_linked_mutations.py || echo ${projectDir}/bin/summarize_linked_mutations.py) \\
        --viloca-input viloca \\
        ${samplesheet_arg} \\
        --out-dir . \\
        --min-reads ${params.viloca_min_reads} \\
        --min-support ${params.viloca_min_pair_support} \\
        --min-samples ${params.viloca_min_pair_samples}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch linked_mutations_long.csv
    touch linked_mutations_recurrent.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process HAPLOTYPE_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path frequency_csv
    path sequences_fasta
    path samplesheet

    output:
    path "*.png"        , emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samplesheet_arg = samplesheet ? "--samplesheet ${samplesheet}" : ""
    def freq_arg = frequency_csv ? "--frequency-csv ${frequency_csv}" : ""
    def seq_arg = sequences_fasta ? "--sequences-fasta ${sequences_fasta}" : ""
    """
    python3 \$(which plot_haplotypes.py || echo ${projectDir}/bin/plot_haplotypes.py) \\
        ${freq_arg} \\
        ${seq_arg} \\
        ${samplesheet_arg} \\
        --dataset ${params.dataset} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.dataset}_haplotype_frequencies.png
    touch ${params.dataset}_haplotype_network.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
