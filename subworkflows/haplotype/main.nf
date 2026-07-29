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
    ch_versions = Channel.empty()
    ch_haplotypes = Channel.empty()

    if (params.run_haplotype) {
        CLIQUESNV(ch_viral_bams)
        ch_versions = ch_versions.mix(CLIQUESNV.out.versions)

        VILOCA(ch_viral_bams, ch_fasta)
        ch_versions = ch_versions.mix(VILOCA.out.versions)

        ch_haplotypes = CLIQUESNV.out.haplotypes.mix(VILOCA.out.csv)
    }

    emit:
    haplotypes = ch_haplotypes
    versions   = ch_versions
}
