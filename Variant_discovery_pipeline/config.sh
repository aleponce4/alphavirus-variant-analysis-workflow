#!/bin/bash

# Pipeline Configuration for iVar and LoFreq
# Modify these values as needed

# Input mode controls pipeline entrypoint
# Options:
#   bam  - start from pre-aligned BAM files in Input/BAMs
#   fastq - start from raw FASTQ reads in INPUT_FASTQ_DIR and map with bwa-mem2
PIPELINE_INPUT_MODE="${PIPELINE_INPUT_MODE:-bam}"

# Input FASTQ path (used when PIPELINE_INPUT_MODE=fastq)
INPUT_FASTQ_DIR="${INPUT_FASTQ_DIR:-Input/FASTQ}"
STAGE1_DIR="${STAGE1_DIR:-work/stage1}"
STAGE1_MANIFEST="${STAGE1_MANIFEST:-${STAGE1_DIR}/manifest/sample_inputs.tsv}"

# iVar settings
IVAR_MIN_VARIANT_DEPTH=5000      # require >=5000x total coverage at a site before evaluating variants
IVAR_MIN_VARIANT_FREQ=0.001      # report alleles present at >=0.1% frequency
IVAR_MIN_BASE_QUALITY=30         # ignore bases with Phred quality <30 (~0.1% error rate)
IVAR_MIN_CONSENSUS_COVERAGE=10   # require >=10x depth to emit consensus bases
IVAR_PRIMER_TRIM_QUALITY=20      # trimming primer regions, drop bases with Q<20
IVAR_PRIMER_BASE_QUALITY=30      # bases used in primer-trimmed regions must have Q>=30
IVAR_PRIMER_WINDOW_SIZE=4        # sliding window size for primer trimming

# LoFreq settings
LOFREQ_MIN_VARIANT_DEPTH=5000    # require >=5000x total coverage at a site before evaluating variants
LOFREQ_MIN_VARIANT_FREQ=0.001    # report alleles >=0.1% frequency after statistical testing
LOFREQ_MIN_BASE_QUALITY=30       # require alt-supporting bases to have Phred Q>=30 (~0.1% error rate)
LOFREQ_MIN_MAP_QUALITY=60        # use only uniquely mapped reads with MAPQ 60 (99.9999% confident in mapping)
LOFREQ_CALL_SIG=0.01             # raw p-value cutoff from LoFreq's binomial test
LOFREQ_ENABLE_INDELQUAL=1        # correct base alignment quality around indels
LOFREQ_BAQ=1                     # enable BAQ adjustment, suppresses false SNP calls caused by alignment artifacts

# General settings
THREADS="${THREADS:-20}"         # number of threads to use for parallelizable steps
BAM_OUTPUT_DIR="${BAM_OUTPUT_DIR:-Input/BAMs}"      # default BAM output folder (used by mapping + variant stages)

# Optional fixed reference path for pipeline-level mapping + annotation scripts
REFERENCE_FASTA="${REFERENCE_FASTA:-Input/Reference/inh.fasta}"
ANNOTATION_GFF="${ANNOTATION_GFF:-Input/Reference/target_reference.gff3}"

# Mapping (bwa-mem2) settings
BWAMEM2_THREADS="${THREADS}"     # thread count for bwa-mem2 and post-processing
BWAMEM2_EXTRA_ARGS="${BWAMEM2_EXTRA_ARGS:-}"  # append additional bwa-mem2 args here, e.g. "-k 19"
FORCE_REMAP="${FORCE_REMAP:-0}"   # set to 1 to rebuild BAM outputs even if already present
BWAMEM2_ENV_NAME="${BWAMEM2_ENV_NAME:-bwa_mem2_env}"  # conda environment name containing bwa-mem2

# BAM indexing can be expensive when many samples are present.
# 1 = skip re-mapping existing BAMs, 0 = re-map all BAMs every run.
SKIP_EXISTING_BAM_MAP="${SKIP_EXISTING_BAM_MAP:-1}"

# Target contig restriction
TARGET_CONTIG="${TARGET_CONTIG:-target_contig}"   # Restrict to Target reads before processing

