/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IVAR_TSV_TO_VCF } from '../../modules/local/ivar/tsv_to_vcf/main'
include { BCFTOOLS_CSQ    } from '../../modules/local/bcftools/csq/main'

workflow ANNOTATION {
    take:
    ch_ivar_tsv   // channel: [ val(meta), path(tsv) ]
    ch_lofreq_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gff        // channel: [ val(meta), path(gff) ]

    main:
    ch_versions = Channel.empty()
    ch_annotated_vcf = Channel.empty()

    if (params.run_annotation) {
        // Convert iVar TSV to VCF
        ch_ivar_vcf = Channel.empty()
        if (params.run_ivar) {
            IVAR_TSV_TO_VCF(ch_ivar_tsv, ch_fasta)
            ch_ivar_vcf = IVAR_TSV_TO_VCF.out.vcf
                .map { meta, vcf -> [ meta + [caller: 'ivar'], vcf ] }
            ch_versions = ch_versions.mix(IVAR_TSV_TO_VCF.out.versions)
        }

        // Prepare LoFreq VCF stream
        ch_lofreq_prep = Channel.empty()
        if (params.run_lofreq) {
            ch_lofreq_prep = ch_lofreq_vcf
                .map { meta, vcf, tbi -> [ meta + [caller: 'lofreq'], vcf ] }
        }

        ch_vcf_to_annotate = ch_ivar_vcf.mix(ch_lofreq_prep)

        BCFTOOLS_CSQ(ch_vcf_to_annotate, ch_fasta, ch_gff)
        ch_annotated_vcf = BCFTOOLS_CSQ.out.vcf
        ch_versions      = ch_versions.mix(BCFTOOLS_CSQ.out.versions)
    }

    emit:
    annotated_vcf = ch_annotated_vcf
    versions      = ch_versions
}
