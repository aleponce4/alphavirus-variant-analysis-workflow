/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REPORT_RUN_SUMMARY    } from '../../modules/local/report/run_summary/main'
include { REPORT_VARIANT_PLOTS  } from '../../modules/local/report/variant_plots/main'
include { REPORT_COVERAGE_PLOTS } from '../../modules/local/report/coverage_plots/main'
include { REPORT_HAPLOTYPE_PLOTS} from '../../modules/local/report/haplotype_plots/main'
include { REPORT_EXCEL_EXPORT   } from '../../modules/local/report/excel_export/main'

workflow REPORTING {
    take:
    ch_results_dir // channel: path(outdir)

    main:
    ch_versions = Channel.empty()

    REPORT_RUN_SUMMARY(ch_results_dir)
    ch_versions = ch_versions.mix(REPORT_RUN_SUMMARY.out.versions)

    REPORT_VARIANT_PLOTS(ch_results_dir)
    ch_versions = ch_versions.mix(REPORT_VARIANT_PLOTS.out.versions)

    REPORT_COVERAGE_PLOTS(ch_results_dir)
    ch_versions = ch_versions.mix(REPORT_COVERAGE_PLOTS.out.versions)

    REPORT_HAPLOTYPE_PLOTS(ch_results_dir)
    ch_versions = ch_versions.mix(REPORT_HAPLOTYPE_PLOTS.out.versions)

    REPORT_EXCEL_EXPORT(ch_results_dir)
    ch_versions = ch_versions.mix(REPORT_EXCEL_EXPORT.out.versions)

    emit:
    versions = ch_versions
}
