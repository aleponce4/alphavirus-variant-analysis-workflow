/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SELECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_GFF3_TO_GTF                                } from '../../modules/local/snpgenie/prep_gtf/main'
include { SNPGENIE_RUN                                       } from '../../modules/local/snpgenie/run/main'
include { SNPGENIE_SUMMARIZE; SELECTION_DELTA; SELECTION_LIMMA; SELECTION_BUILD_TABLES } from '../../modules/local/snpgenie/downstream/main'

workflow SELECTION {
    take:
    ch_vcf             // channel: [ val(meta), path(vcf) ]
    ch_fasta           // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gff             // channel: [ val(meta), path(gff) ]
    ch_samplesheet     // path: samplesheet.csv
    ch_treatment_groups// channel: [ treatment: treatment, n_replicates: size ]

    main:
    ch_versions = Channel.empty()
    ch_selection_tables = Channel.empty()

    if (params.run_snpgenie) {
        ch_treatment_groups.subscribe { log.info "SELECTION treatment group: ${it.treatment} (n=${it.n_replicates})" }

        CONVERT_GFF3_TO_GTF(ch_gff)
        ch_versions = ch_versions.mix(CONVERT_GFF3_TO_GTF.out.versions)

        SNPGENIE_RUN(ch_vcf, ch_fasta, CONVERT_GFF3_TO_GTF.out.gtf)
        ch_versions = ch_versions.mix(SNPGENIE_RUN.out.versions)

        ch_snpgenie_files = SNPGENIE_RUN.out.results
            .map { meta, files -> files }
            .flatten()
            .collect()

        SNPGENIE_SUMMARIZE(ch_snpgenie_files)
        ch_versions = ch_versions.mix(SNPGENIE_SUMMARIZE.out.versions)

        SELECTION_DELTA(SNPGENIE_SUMMARIZE.out.product_results, ch_samplesheet)
        ch_versions = ch_versions.mix(SELECTION_DELTA.out.versions)

        SELECTION_LIMMA(SELECTION_DELTA.out.delta_per_sample, ch_samplesheet)
        ch_versions = ch_versions.mix(SELECTION_LIMMA.out.versions)

        SELECTION_BUILD_TABLES(
            SELECTION_DELTA.out.kruskal_by_gene,
            SELECTION_LIMMA.out.limma_overall,
            SELECTION_LIMMA.out.limma_contrasts,
            SNPGENIE_SUMMARIZE.out.product_results,
            ch_samplesheet
        )
        ch_selection_tables = SELECTION_BUILD_TABLES.out.key_table
        ch_versions         = ch_versions.mix(SELECTION_BUILD_TABLES.out.versions)
    }

    emit:
    selection_tables = ch_selection_tables
    versions         = ch_versions
}
