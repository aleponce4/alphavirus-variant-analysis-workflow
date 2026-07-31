process VILOCA {
    tag "$meta.id"
    label 'process_high'

    container "${params.container_viloca}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta), path(fai)

    output:
    tuple val(meta), path("haplotypes/*.fasta")       , emit: local_haplotypes, optional: true
    tuple val(meta), path("${meta.id}_cooccurring_mutations.csv"), emit: cooccurrence
    tuple val(meta), path("coverage.txt")             , emit: coverage, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    viloca run \\
        -b ${bam} \\
        -f ${fasta} \\
        -w ${params.viloca_window} \\
        -t ${task.cpus} \\
        --mode use_quality_scores \\
        --win_min_ext 0.85 \\
        --min_windows_coverage 10 \\
        --exclude_non_var_pos_threshold 0.005

    mkdir -p haplotypes
    if [ -d "snv/haplotypes" ]; then
        mv snv/haplotypes/*.fasta haplotypes/ 2>/dev/null || true
    fi
    if [ -f "snv/cooccurring_mutations.csv" ]; then
        mv snv/cooccurring_mutations.csv ${prefix}_cooccurring_mutations.csv
    elif [ -f "cooccurring_mutations.csv" ]; then
        mv cooccurring_mutations.csv ${prefix}_cooccurring_mutations.csv
    else
        touch ${prefix}_cooccurring_mutations.csv
    fi
    if [ -f "snv/coverage.txt" ]; then
        mv snv/coverage.txt coverage.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viloca: \$(viloca --version 2>&1 | sed 's/viloca //' || echo "1.1.0")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p haplotypes
    touch haplotypes/${prefix}_haplotype_stub.fasta
    touch ${prefix}_cooccurring_mutations.csv
    touch coverage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viloca: \$(viloca --version 2>&1 | sed 's/viloca //' || echo "1.1.0")
    END_VERSIONS
    """
}
