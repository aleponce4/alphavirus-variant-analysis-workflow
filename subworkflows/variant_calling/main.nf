/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT_CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_VIRAL_BAM } from '../../modules/local/extract_viral_bam/main'
include { IVAR_VARIANTS     } from '../../modules/local/ivar/variants/main'
include { IVAR_CONSENSUS    } from '../../modules/local/ivar/consensus/main'
include { LOFREQ_CALL       } from '../../modules/local/lofreq/call/main'
include { LOFREQ_FILTER     } from '../../modules/local/lofreq/filter/main'

workflow VARIANT_CALLING {
    take:
    ch_samples // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta   // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gff     // channel: [ val(meta), path(gff) ]

    main:
    ch_versions = Channel.empty()

    // 1. Extract viral BAMs
    EXTRACT_VIRAL_BAM(ch_samples, params.viral_contig ?: '')
    ch_viral_bams = EXTRACT_VIRAL_BAM.out.bam
    ch_versions   = ch_versions.mix(EXTRACT_VIRAL_BAM.out.versions)

    // 2. iVar variants branch
    ch_ivar_tsv       = Channel.empty()
    ch_ivar_consensus = Channel.empty()
    if (params.run_ivar) {
        IVAR_VARIANTS(ch_viral_bams, ch_fasta.first(), ch_gff.first())
        ch_ivar_tsv = IVAR_VARIANTS.out.tsv
        ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions)

        IVAR_CONSENSUS(ch_viral_bams, ch_fasta.first())
        ch_ivar_consensus = IVAR_CONSENSUS.out.consensus
        ch_versions       = ch_versions.mix(IVAR_CONSENSUS.out.versions)
    }

    // 3. LoFreq branch
    ch_lofreq_vcf     = Channel.empty()
    ch_lofreq_qc      = Channel.empty()
    if (params.run_lofreq) {
        LOFREQ_CALL(ch_viral_bams, ch_fasta.first())
        ch_versions = ch_versions.mix(LOFREQ_CALL.out.versions)

        LOFREQ_FILTER(LOFREQ_CALL.out.vcf)
        ch_lofreq_vcf = LOFREQ_FILTER.out.vcf
        ch_lofreq_qc  = LOFREQ_FILTER.out.qc_stats
        ch_versions   = ch_versions.mix(LOFREQ_FILTER.out.versions)
    }

    emit:
    viral_bams     = ch_viral_bams
    ivar_tsv       = ch_ivar_tsv
    ivar_consensus = ch_ivar_consensus
    lofreq_vcf     = ch_lofreq_vcf
    lofreq_qc      = ch_lofreq_qc
    versions       = ch_versions
}
