include { FASTQC         } from '../../modules/local/fastqc/main'
include { BWA_INDEX      } from '../../modules/local/bwa/index/main'
include { BWA_MEM        } from '../../modules/local/bwa/mem/main'
include { SAMTOOLS_SORT  } from '../../modules/local/samtools/sort/main'
include { IVAR_TRIM      } from '../../modules/local/ivar/trim/main'
include { SAMTOOLS_STATS } from '../../modules/local/samtools/stats/main'

workflow READ_PREPROCESSING {
    take:
    ch_fastqs   // [ val(meta), path(fastq_1), path(fastq_2) ]
    ch_fasta    // [ val(meta_ref), path(fasta), path(fai) ]

    main:
    ch_versions = Channel.empty()

    FASTQC(ch_fastqs)
    BWA_INDEX(ch_fasta.map { meta, fasta, fai -> [ meta, fasta ] })
    BWA_MEM(ch_fastqs, BWA_INDEX.out.index)
    SAMTOOLS_SORT(BWA_MEM.out.sam)

    ch_bams = SAMTOOLS_SORT.out.bam

    if (params.protocol == 'amplicon') {
        IVAR_TRIM(SAMTOOLS_SORT.out.bam, file(params.primer_bed, checkIfExists: true))
        ch_bams = IVAR_TRIM.out.bam
        ch_versions = ch_versions.mix(IVAR_TRIM.out.versions)
    }

    SAMTOOLS_STATS(ch_bams)

    ch_versions = ch_versions
        .mix(FASTQC.out.versions)
        .mix(BWA_INDEX.out.versions)
        .mix(BWA_MEM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_STATS.out.versions)

    emit:
    bams        = ch_bams                  // [ meta, bam, bai ]
    fastqc_zip  = FASTQC.out.zip
    bwa_log     = BWA_MEM.out.log
    stats       = SAMTOOLS_STATS.out.stats
    versions    = ch_versions
}
