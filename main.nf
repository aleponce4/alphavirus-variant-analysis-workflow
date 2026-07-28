#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alphavirus-variant-analysis Main Workflow Entrypoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// Subworkflow imports (uncomment as subworkflows are activated)
/*
include { INPUT_CHECK     } from './subworkflows/input_check/main'
include { VARIANT_CALLING } from './subworkflows/variant_calling/main'
include { ANNOTATION      } from './subworkflows/annotation/main'
include { COVERAGE_QC     } from './subworkflows/coverage_qc/main'
include { SELECTION       } from './subworkflows/selection/main'
include { HAPLOTYPE       } from './subworkflows/haplotype/main'
include { REPORTING       } from './subworkflows/reporting/main'
*/

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

    // Phase 1 stub workflow execution complete
}
