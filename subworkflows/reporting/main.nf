/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REPORT_RUN_SUMMARY; REPORT_COVERAGE_PLOTS; REPORT_VARIANT_PLOTS } from '../../modules/local/reporting/main'
include { EXECUTIVE_REPORT                                                 } from '../../modules/local/report/executive_report/main'

workflow REPORTING {
    take:
    ch_lofreq_qc         // channel: path(qc_stats)
    ch_coverage_depth    // channel: path(depth)
    ch_vcf               // channel: path(vcf)
    ch_selection_tables  // channel: path(key_table)
    ch_haplotypes        // channel: path(haplotypes)
    ch_software_versions // channel: path(software_versions.yml)

    main:
    ch_versions = Channel.empty()

    REPORT_RUN_SUMMARY(ch_lofreq_qc.collect())
    ch_versions = ch_versions.mix(REPORT_RUN_SUMMARY.out.versions)

    REPORT_COVERAGE_PLOTS(ch_coverage_depth.collect())
    ch_versions = ch_versions.mix(REPORT_COVERAGE_PLOTS.out.versions)

    REPORT_VARIANT_PLOTS(ch_vcf.collect())
    ch_versions = ch_versions.mix(REPORT_VARIANT_PLOTS.out.versions)

    EXECUTIVE_REPORT(
        REPORT_RUN_SUMMARY.out.summary,
        REPORT_COVERAGE_PLOTS.out.plots,
        REPORT_VARIANT_PLOTS.out.plots,
        ch_selection_tables.ifEmpty([]),
        ch_haplotypes.ifEmpty([]),
        ch_software_versions,
        file("${projectDir}/assets/executive_report.qmd")
    )
    ch_versions = ch_versions.mix(EXECUTIVE_REPORT.out.versions)

    emit:
    summary  = REPORT_RUN_SUMMARY.out.summary
    report   = EXECUTIVE_REPORT.out.html
    versions = ch_versions
}
