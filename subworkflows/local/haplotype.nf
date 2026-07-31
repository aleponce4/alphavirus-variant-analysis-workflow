/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: HAPLOTYPE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLIQUESNV } from '../../modules/local/cliquesnv/main'
include { VILOCA    } from '../../modules/local/viloca/main'

workflow HAPLOTYPE {
    take:
    ch_viral_bams // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta      // channel: [ val(meta), path(fasta), path(fai) ]

    main:
    ch_versions                = Channel.empty()
    ch_cliquesnv_fasta         = Channel.empty()
    ch_viloca_cooccurrence     = Channel.empty()
    ch_viloca_coverage         = Channel.empty()
    ch_viloca_local_haplotypes = Channel.empty()

    if (params.run_haplotype) {
        CLIQUESNV(ch_viral_bams)
        ch_versions = ch_versions.mix(CLIQUESNV.out.versions)
        ch_cliquesnv_fasta = CLIQUESNV.out.haplotypes

        VILOCA(ch_viral_bams, ch_fasta.first())
        ch_versions = ch_versions.mix(VILOCA.out.versions)
        ch_viloca_cooccurrence = VILOCA.out.cooccurrence
        ch_viloca_coverage = VILOCA.out.coverage
        ch_viloca_local_haplotypes = VILOCA.out.local_haplotypes
    }

    emit:
    cliquesnv_fasta         = ch_cliquesnv_fasta
    viloca_cooccurrence     = ch_viloca_cooccurrence
    viloca_coverage         = ch_viloca_coverage
    viloca_local_haplotypes = ch_viloca_local_haplotypes
    versions                = ch_versions
}
