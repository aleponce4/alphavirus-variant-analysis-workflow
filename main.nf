#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alphavirus-variant-analysis Main Workflow Entrypoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { INPUT_CHECK            } from './subworkflows/input_check/main'
include { READ_PREPROCESSING     } from './subworkflows/local/read_preprocessing'
include { VARIANT_CALLING        } from './subworkflows/variant_calling/main'
include { ANNOTATION             } from './subworkflows/annotation/main'
include { COVERAGE_QC            } from './subworkflows/coverage_qc/main'
include { SELECTION              } from './subworkflows/local/selection'
include { HAPLOTYPE              } from './subworkflows/local/haplotype'
include { REPORTING              } from './subworkflows/reporting/main'
include { DUMP_SOFTWARE_VERSIONS } from './modules/local/dumpsoftwareversions/main'
include { MULTIQC                } from './modules/local/multiqc/main'

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)

    log.info """
    ================================================================
    A L P H A V I R U S   V A R I A N T   A N A L Y S I S
    ================================================================
    Dataset         : ${params.dataset}
    Input           : ${params.input}
    FASTA           : ${params.fasta}
    GFF             : ${params.gff}
    Contig          : ${params.viral_contig}
    Protocol        : ${params.protocol}
    Primer BED      : ${params.primer_bed}
    Outdir          : ${params.outdir}
    Run iVar        : ${params.run_ivar}
    Run LoFreq      : ${params.run_lofreq}
    Run Annotation  : ${params.run_annotation}
    Run Coverage    : ${params.run_coverage}
    Run SNPGenie    : ${params.run_snpgenie}
    Run Haplotype   : ${params.run_haplotype}
    ================================================================
    """

    // 1. Input Validation
    INPUT_CHECK(file(params.input, checkIfExists: true))

    // 2. Read Preprocessing & Alignment
    READ_PREPROCESSING(INPUT_CHECK.out.fastqs, INPUT_CHECK.out.fasta)

    // 3. Variant Calling
    VARIANT_CALLING(READ_PREPROCESSING.out.bams, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 4. Annotation
    ANNOTATION(VARIANT_CALLING.out.ivar_tsv, VARIANT_CALLING.out.lofreq_vcf, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 5. Coverage QC
    COVERAGE_QC(VARIANT_CALLING.out.viral_bams)

    // 6. Selection Analysis (optional, --run_snpgenie)
    ch_snpgenie_vcfs = VARIANT_CALLING.out.lofreq_vcf
        .map { meta, vcf, tbi -> [ meta, vcf ] }

    ch_treatment_groups = INPUT_CHECK.out.fastqs
        .map { meta, fastq_1, fastq_2 -> [ meta.treatment, meta.id ] }
        .groupTuple()
        .map { treatment, ids -> [ treatment: treatment, n_replicates: ids.size() ] }

    SELECTION(
        ch_snpgenie_vcfs,
        INPUT_CHECK.out.fasta,
        INPUT_CHECK.out.gff,
        file(params.input),
        ch_treatment_groups
    )

    // 7. Haplotype Reconstruction (optional, --run_haplotype)
    HAPLOTYPE(VARIANT_CALLING.out.viral_bams, INPUT_CHECK.out.fasta)

    // 8. Software Versions Aggregation
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(READ_PREPROCESSING.out.versions)
        .mix(VARIANT_CALLING.out.versions)
        .mix(ANNOTATION.out.versions)
        .mix(COVERAGE_QC.out.versions)
        .mix(SELECTION.out.versions)
        .mix(HAPLOTYPE.out.versions)

    DUMP_SOFTWARE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml', newLine: true))

    // 9. Executive HTML & Biological Reporting
    REPORTING(
        VARIANT_CALLING.out.lofreq_qc.map { meta, qc -> qc },
        COVERAGE_QC.out.coverage_summary.map { meta, summary -> summary },
        VARIANT_CALLING.out.lofreq_vcf.map { meta, vcf, tbi -> vcf },
        SELECTION.out.selection_tables,
        HAPLOTYPE.out.haplotypes.map { meta, fasta -> fasta },
        DUMP_SOFTWARE_VERSIONS.out.yml
    )

    // 10. MultiQC Terminal Dashboard
    ch_multiqc_files = Channel.empty()
        .mix(READ_PREPROCESSING.out.fastqc_zip.map { meta, zip -> zip })
        .mix(READ_PREPROCESSING.out.bwa_log.map    { meta, log -> log })
        .mix(READ_PREPROCESSING.out.stats.map      { meta, stats, flagstat, idxstats -> [ stats, flagstat, idxstats ] })
        .mix(VARIANT_CALLING.out.lofreq_qc.map     { meta, qc -> qc })
        .mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml)
        .collect()

    MULTIQC(ch_multiqc_files, file("${projectDir}/assets/multiqc_config.yml"))
}
