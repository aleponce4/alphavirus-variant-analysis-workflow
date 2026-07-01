#!/bin/bash

# TEST CONFIGURATION FILE
# Use these settings for quick pipeline validation (runs in ~minutes instead of hours)
# Source this file instead of config.sh to override with smaller thresholds

# ========== iVar PARAMETERS (Frequency-based) ==========
IVAR_MIN_CONSENSUS_COVERAGE=100        # [TEST] Reduced from 10000 - calls consensus faster
IVAR_MIN_BASE_QUALITY=30               # Quality threshold for bases
IVAR_MIN_VARIANT_DEPTH=100             # [TEST] Reduced from 5000 - detects variants faster
IVAR_MIN_VARIANT_FREQ=0.001            # Minimum 0.1% frequency (same as prod)
IVAR_PRIMER_TRIM_QUALITY=30            # Base quality for primer trimming
IVAR_PRIMER_BASE_QUALITY=20            # Quality for primer bases
IVAR_PRIMER_WINDOW_SIZE=4              # Sliding window size for primer trimming

# ========== LoFreq PARAMETERS (Statistical binomial test) ==========
LOFREQ_MIN_VARIANT_DEPTH=100           # [TEST] Reduced from 5000 - faster calls
LOFREQ_MIN_BASE_QUALITY=30             # Min quality for alt alleles
LOFREQ_MIN_MAP_QUALITY=60              # Only uniquely mapped reads (MAPQ 60)
LOFREQ_CALL_SIG=0.05                   # Significance threshold for calling
LOFREQ_ENABLE_INDELQUAL=1              # Correct indel alignment quality
LOFREQ_BAQ=1                           # Viterbi realignment (fixes alignment artifacts)

# ========== GENERAL SETTINGS ==========
TARGET_CONTIG="target_contig"                # Restrict to Target genome only (exclude host)
THREADS=4                              # [TEST] Reduced from 20 for quick testing

# ========== TEST-SPECIFIC NOTES ==========
# This configuration is designed to:
# 1. Use a subset of BAM files (manually select in test_pipeline.sh)
# 2. Call variants at lower coverage (100× instead of 5000×) for speed
# 3. Still use reasonable quality thresholds to catch real variants
# 4. Produce valid VCF output for annotation testing
#
# Expected runtime: 2-5 minutes for 2 samples
# Output: Same structure as production (Ivar/, LoFreq/, Annotated_variants/)

