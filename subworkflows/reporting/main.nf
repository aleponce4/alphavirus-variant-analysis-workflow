/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REPORT_RUN_SUMMARY; REPORT_COVERAGE_PLOTS; REPORT_VARIANT_PLOTS } from '../../modules/local/reporting/main'
include { EXECUTIVE_REPORT                                                 } from '../../modules/local/report/executive_report/main'
include { HAPLOTYPE_TABLES; HAPLOTYPE_LINKAGE; HAPLOTYPE_PLOTS             } from '../../modules/local/report/haplotype_report/main'

workflow REPORTING {
    take:
    ch_lofreq_qc           // channel: path(qc_stats)
    ch_coverage_depth      // channel: path(depth)
    ch_vcf                 // channel: path(vcf)
    ch_selection_tables    // channel: path(key_table)
    ch_cliquesnv_fasta     // channel: [meta, fasta]
    ch_viloca_cooccurrence // channel: [meta, csv]
    ch_fasta               // channel: [meta, fasta, fai]
    ch_gff                 // channel: [meta, gff]
    ch_samplesheet         // channel: path(samplesheet)
    ch_software_versions   // channel: path(software_versions.yml)

    main:
    ch_versions = Channel.empty()

    REPORT_RUN_SUMMARY(ch_lofreq_qc.collect())
    ch_versions = ch_versions.mix(REPORT_RUN_SUMMARY.out.versions)

    REPORT_COVERAGE_PLOTS(ch_coverage_depth.collect())
    ch_versions = ch_versions.mix(REPORT_COVERAGE_PLOTS.out.versions)

    REPORT_VARIANT_PLOTS(ch_vcf.collect())
    ch_versions = ch_versions.mix(REPORT_VARIANT_PLOTS.out.versions)

    // Haplotype Reporting
    ch_cliquesnv_files = ch_cliquesnv_fasta.map { meta, fasta -> fasta }.collect().ifEmpty([])
    ch_viloca_files = ch_viloca_cooccurrence.map { meta, csv -> csv }.collect().ifEmpty([])
    ch_ref_fasta = ch_fasta.map { meta, fasta, fai -> fasta }
    ch_ref_gff = ch_gff.map { meta, gff -> gff }

    HAPLOTYPE_TABLES(
        ch_cliquesnv_files,
        ch_ref_fasta,
        ch_ref_gff,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(HAPLOTYPE_TABLES.out.versions)

    HAPLOTYPE_LINKAGE(
        ch_viloca_files,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(HAPLOTYPE_LINKAGE.out.versions)

    HAPLOTYPE_PLOTS(
        HAPLOTYPE_TABLES.out.frequency_table,
        HAPLOTYPE_TABLES.out.sequences_fasta,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(HAPLOTYPE_PLOTS.out.versions)

    ch_haplotype_artifacts = Channel.empty()
        .mix(HAPLOTYPE_PLOTS.out.plots)
        .mix(HAPLOTYPE_TABLES.out.summary_table)
        .mix(HAPLOTYPE_LINKAGE.out.recurrent_table)
        .collect()

    EXECUTIVE_REPORT(
        REPORT_RUN_SUMMARY.out.summary,
        REPORT_COVERAGE_PLOTS.out.plots,
        REPORT_VARIANT_PLOTS.out.plots,
        ch_selection_tables.ifEmpty([]),
        ch_haplotype_artifacts.ifEmpty([]),
        ch_software_versions,
        file("${projectDir}/assets/executive_report.qmd")
    )
    ch_versions = ch_versions.mix(EXECUTIVE_REPORT.out.versions)

    emit:
    summary  = REPORT_RUN_SUMMARY.out.summary
    report   = EXECUTIVE_REPORT.out.html
    versions = ch_versions
}
