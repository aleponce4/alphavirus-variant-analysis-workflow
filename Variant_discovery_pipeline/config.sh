#!/bin/bash

# Pipeline Configuration for iVar and LoFreq
# Modify these values as needed

# iVar settings
IVAR_MIN_VARIANT_DEPTH=5000      # require ≥5000× total coverage at a site before evaluating variants
IVAR_MIN_VARIANT_FREQ=0.001      # report alleles present at ≥0.1% frequency
IVAR_MIN_BASE_QUALITY=30         # ignore bases with Phred quality <30 (~0.1% error rate)
IVAR_MIN_CONSENSUS_COVERAGE=10   # require ≥10× depth to emit consensus bases
IVAR_PRIMER_TRIM_QUALITY=20      # trimming primer regions, drop bases with Q<20
IVAR_PRIMER_BASE_QUALITY=30      # bases used in primer-trimmed regions must have Q≥30
IVAR_PRIMER_WINDOW_SIZE=4        # sliding window size for primer trimming

# LoFreq settings
LOFREQ_MIN_VARIANT_DEPTH=5000    # require ≥5000× total coverage at a site before evaluating variants
LOFREQ_MIN_VARIANT_FREQ=0.001    # report alleles ≥0.1% frequency after statistical testing
LOFREQ_MIN_BASE_QUALITY=30       # require alt-supporting bases to have Phred Q≥30 (~0.1% error rate)
LOFREQ_MIN_MAP_QUALITY=60        # use only uniquely mapped reads with MAPQ 60 (99.9999% confident in mapping)
LOFREQ_CALL_SIG=0.01             # raw p-value cutoff from LoFreq's binomial test
LOFREQ_ENABLE_INDELQUAL=1        # correct base alignment quality around indels
LOFREQ_BAQ=1                     # enable BAQ adjustment, suppresses false SNP calls caused by alignment artifacts

# General settings
THREADS=8                       # number of threads to use for parallelizable steps

# Viral contig restriction
VIRAL_CONTIG="${VIRAL_CONTIG:-VEEV_INH}"   # Restrict to viral reads before processing
