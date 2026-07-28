/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: INPUT_CHECK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_CHECK {
    take:
    samplesheet // path: samplesheet.csv

    main:
    // Parse samplesheet CSV
    ch_samplesheet = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [:]
            meta.id        = row.sample
            meta.condition = row.condition ?: 'unknown'
            meta.dpi       = row.dpi ?: '0'

            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)

            return [ meta, bam, bai ]
        }

    ch_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [ [id: 'reference'], fasta ] }

    ch_gff = Channel
        .fromPath(params.gff, checkIfExists: true)
        .map { gff -> [ [id: 'annotation'], gff ] }

    emit:
    samples = ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ]
    fasta   = ch_fasta       // channel: [ val(meta), path(fasta) ]
    gff     = ch_gff         // channel: [ val(meta), path(gff) ]
}
