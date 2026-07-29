/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: INPUT_CHECK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'
include { SAMTOOLS_FAIDX    } from '../../modules/local/samtools/faidx/main'

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
    if (params.protocol == 'amplicon' && !params.primer_bed) {
        error "Parameter --primer_bed is required when --protocol amplicon."
    }
    if (params.protocol == 'metagenomic' && params.primer_bed) {
        log.warn "--primer_bed supplied but --protocol is 'metagenomic'; primer trimming is OFF."
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

    def root_dir = workflow.projectDir.toString().split('/.nf-test')[0].split('\\.nf-test')[0]

    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row ->
            def meta = [:]
            meta.id        = row[0]
            meta.treatment = row[3]
            meta.condition = row[3]
            def f1 = row[1].toString()
            def f2 = row[2].toString()
            def file1 = f1.startsWith('/') ? file(f1, checkIfExists: true) : (file(f1).exists() ? file(f1) : file("${root_dir}/${f1}", checkIfExists: true))
            def file2 = f2.startsWith('/') ? file(f2, checkIfExists: true) : (file(f2).exists() ? file(f2) : file("${root_dir}/${f2}", checkIfExists: true))
            return [ meta, file1, file2 ]
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
    fastqs   = ch_samplesheet // channel: [ val(meta), path(fastq_1), path(fastq_2) ]
    fasta    = ch_indexed_fasta// channel: [ val(meta), path(fasta), path(fai) ]
    gff      = ch_gff         // channel: [ val(meta), path(gff) ]
    versions = SAMTOOLS_FAIDX.out.versions
}
