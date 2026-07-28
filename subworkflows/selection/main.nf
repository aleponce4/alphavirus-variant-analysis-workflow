/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SELECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPGENIE_RUN } from '../../modules/local/snpgenie/run/main'

workflow SELECTION {
    take:
    ch_lofreq_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gff        // channel: [ val(meta), path(gff) ]

    main:
    ch_versions = Channel.empty()
    ch_selection_results = Channel.empty()

    if (params.run_snpgenie) {
        ch_vcf_input = ch_lofreq_vcf.map { meta, vcf, tbi -> [ meta, vcf ] }

        // GTF prep channel stub
        ch_gtf = ch_gff.map { meta, gff -> [ meta, gff ] }

        SNPGENIE_RUN(ch_vcf_input, ch_fasta, ch_gtf)
        ch_selection_results = SNPGENIE_RUN.out.results
        ch_versions          = ch_versions.mix(SNPGENIE_RUN.out.versions)
    }

    emit:
    results  = ch_selection_results
    versions = ch_versions
}
