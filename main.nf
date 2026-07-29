#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alphavirus-variant-analysis Main Workflow Entrypoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK     } from './subworkflows/input_check/main'
include { VARIANT_CALLING } from './subworkflows/variant_calling/main'
include { ANNOTATION      } from './subworkflows/annotation/main'
include { COVERAGE_QC     } from './subworkflows/coverage_qc/main'
include { SELECTION       } from './subworkflows/selection/main'
include { HAPLOTYPE       } from './subworkflows/haplotype/main'
include { REPORTING       } from './subworkflows/reporting/main'

workflow {
    log.info """
    ================================================================
    A L P H A V I R U S   V A R I A N T   A N A L Y S I S
    ================================================================
    Dataset         : ${params.dataset}
    Input           : ${params.input}
    FASTA           : ${params.fasta}
    GFF             : ${params.gff}
    Contig          : ${params.viral_contig}
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

    // 2. Variant Calling
    VARIANT_CALLING(INPUT_CHECK.out.samples, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 3. Annotation
    ANNOTATION(VARIANT_CALLING.out.ivar_tsv, VARIANT_CALLING.out.lofreq_vcf, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 4. Coverage QC
    COVERAGE_QC(VARIANT_CALLING.out.viral_bams)

    // 5. Selection Analysis (optional)
    SELECTION(VARIANT_CALLING.out.lofreq_vcf.map { meta, vcf, tbi -> [ meta, vcf ] }, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff, file(params.input))

    // 6. Haplotype Reconstruction (optional)
    HAPLOTYPE(VARIANT_CALLING.out.viral_bams, INPUT_CHECK.out.fasta)

    // 7. Reporting
    REPORTING(
        VARIANT_CALLING.out.lofreq_qc.map { meta, qc -> qc },
        COVERAGE_QC.out.coverage_summary.map { meta, summary -> summary },
        VARIANT_CALLING.out.lofreq_vcf.map { meta, vcf, tbi -> vcf }
    )
}
