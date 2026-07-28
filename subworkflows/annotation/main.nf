/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_CSQ } from '../../modules/local/bcftools/csq/main'

workflow ANNOTATION {
    take:
    ch_lofreq_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gff        // channel: [ val(meta), path(gff) ]

    main:
    ch_versions = Channel.empty()
    ch_annotated_vcfs = Channel.empty()

    if (params.run_annotation) {
        // Strip index for BCFTOOLS_CSQ input tuple
        ch_vcf_input = ch_lofreq_vcf.map { meta, vcf, tbi -> [ meta, vcf ] }

        BCFTOOLS_CSQ(ch_vcf_input, ch_fasta, ch_gff)
        ch_annotated_vcfs = BCFTOOLS_CSQ.out.vcf
        ch_versions       = ch_versions.mix(BCFTOOLS_CSQ.out.versions)
    }

    emit:
    annotated_vcfs = ch_annotated_vcfs
    versions       = ch_versions
}
