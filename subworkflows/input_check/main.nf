/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: INPUT_CHECK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX } from '../../modules/local/samtools/faidx/main'

workflow INPUT_CHECK {
    take:
    samplesheet // path: samplesheet.csv

    main:
    if (!params.fasta) {
        error "Parameter 'params.fasta' must be specified."
    }
    if (!params.gff) {
        error "Parameter 'params.gff' must be specified."
    }

    def gff_file = file(params.gff, checkIfExists: true)
    def has_cds = false
    gff_file.eachLine { line ->
        if (!line.startsWith('#') && line.contains('\tCDS\t')) {
            has_cds = true
        }
    }
    if (!has_cds) {
        error "Provided GFF file ${params.gff} does not contain CDS features required for annotation and selection analysis."
    }

    ch_samplesheet = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [:]
            meta.id        = row.sample
            meta.condition = row.condition ?: 'infected'
            meta.dpi       = row.dpi ?: '0'

            def bam_path = row.bam.startsWith('/') ? row.bam : "${workflow.projectDir}/${row.bam}"
            def bai_path = row.bai.startsWith('/') ? row.bai : "${workflow.projectDir}/${row.bai}"

            def bam = file(bam_path, checkIfExists: true)
            def bai = file(bai_path, checkIfExists: true)

            return [ meta, bam, bai ]
        }

    ch_raw_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [ [id: 'reference'], fasta ] }

    SAMTOOLS_FAIDX(ch_raw_fasta)
    ch_indexed_fasta = SAMTOOLS_FAIDX.out.fai

    ch_gff = Channel
        .fromPath(params.gff, checkIfExists: true)
        .map { gff -> [ [id: 'annotation'], gff ] }

    emit:
    samples = ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ]
    fasta   = ch_indexed_fasta// channel: [ val(meta), path(fasta), path(fai) ]
    gff     = ch_gff         // channel: [ val(meta), path(gff) ]
}
