#!/bin/bash
# =====================================================================
# config.sh — Adapted configuration for nf-core RNAseq BAM inputs
# 
# Key differences from original workflow:
#   - No primer trimming (RNA-seq, not amplicon)
#   - VIRAL_CONTIG is dataset-dependent (KP282671.1 or EEEV_FL93)
#   - BAMs are markdup.sorted from STAR (already sorted+indexed)
#   - Coverage may be uneven (RNA-seq bias)
# =====================================================================

# ========== DATASET SELECTION ==========
# Set DATASET before sourcing this file, or default to mouse_veev
# Options: mouse_veev | mouse_eeev | rat_veev
DATASET="${DATASET:-mouse_veev}"

case "$DATASET" in
  mouse_veev|rat_veev)
    VIRAL_CONTIG="KP282671.1"
    VIRUS="VEEV"
    ;;
  mouse_eeev)
    VIRAL_CONTIG="EEEV_FL93"
    VIRUS="EEEV"
    ;;
  *)
    echo "ERROR: Unknown DATASET '$DATASET'"
    echo "Set DATASET to one of: mouse_veev, mouse_eeev, rat_veev"
    exit 1
    ;;
esac

# ========== PATHS ==========
# config.sh is at variant_analysis/workflow/config.sh
# PROJECT_ROOT = 045_048_nftcore/ (parent of variant_analysis/)
VARIANT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="$(cd "$VARIANT_ROOT/.." && pwd)"
INPUT_DIR="${VARIANT_ROOT}/input/${DATASET}/BAMs"
REFERENCE="${VARIANT_ROOT}/references/${VIRUS}_INH/viral_only.fasta"
ANNOTATION="${VARIANT_ROOT}/references/${VIRUS}_INH/viral_only.gff3"

# Results go to WSL ext4 (SSD) to avoid 9p filesystem race conditions with parallel writes
# Falls back to VARIANT_ROOT if WSL_HOME is not available
if [ -d "$HOME" ] && [ "$HOME" != "/root" ]; then
  RESULTS_DIR="${HOME}/variant_results/${DATASET}"
else
  RESULTS_DIR="${VARIANT_ROOT}/results/${DATASET}"
fi
mkdir -p "$RESULTS_DIR"

# Note: For EEEV, the reference dir is EEEV_FL93, not EEEV_INH
if [ "$VIRUS" = "EEEV" ]; then
  REFERENCE="${VARIANT_ROOT}/references/EEEV_FL93/viral_only.fasta"
  ANNOTATION="${VARIANT_ROOT}/references/EEEV_FL93/viral_only.gff3"
fi

# ========== iVar SETTINGS ==========
IVAR_MIN_VARIANT_DEPTH=1000      # Require ≥1000× coverage at variant sites
IVAR_MIN_VARIANT_FREQ=0.01       # Report alleles ≥1% frequency
IVAR_MIN_BASE_QUALITY=30         # Ignore bases with Phred quality <30
IVAR_MIN_CONSENSUS_COVERAGE=10   # Require ≥10× depth for consensus bases
IVAR_PRIMER_TRIM_QUALITY=20      # (unused for RNA-seq, no primers)
IVAR_PRIMER_BASE_QUALITY=30      # (unused)
IVAR_PRIMER_WINDOW_SIZE=4        # (unused)

# ========== LoFreq SETTINGS ==========
LOFREQ_MIN_VARIANT_DEPTH=1000    # Require ≥1000× coverage at variant sites
LOFREQ_MIN_VARIANT_FREQ=0.01     # Report alleles ≥1% frequency
LOFREQ_MIN_BASE_QUALITY=30       # Require alt-supporting bases Q≥30
LOFREQ_MIN_MAP_QUALITY=60        # Only uniquely mapped reads (MAPQ 60)
LOFREQ_CALL_SIG=0.01             # Raw p-value cutoff from binomial test
LOFREQ_ENABLE_INDELQUAL=0        # Disabled for RNA-seq (STAR spliced alignments crash lofreq indelqual)
LOFREQ_BAQ=0                     # Disabled for RNA-seq (STAR spliced alignments crash lofreq viterbi)

# ========== GENERAL SETTINGS ==========
# For parallel execution on 16-core/32-thread Threadripper:
#   6 parallel jobs × 3 threads = 18 threads (leaves headroom)
# For single-sample runs, set THREADS=20 for max speed
THREADS="${THREADS:-3}"           # Per-job threads (override with env var)
MAX_JOBS="${MAX_JOBS:-6}"        # Max parallel sample jobs

# ========== OUTPUT DIRECTORIES ==========
LOFREQ_OUT="${RESULTS_DIR}/LoFreq"
IVAR_OUT="${RESULTS_DIR}/Ivar"
ANNOTATED_OUT="${RESULTS_DIR}/Annotated_variants"
COVERAGE_OUT="${RESULTS_DIR}/Coverage"

# ========== DISPLAY ==========
echo "═══════════════════════════════════════════════════════════"
echo "  Dataset:       $DATASET"
echo "  Virus:         $VIRUS"
echo "  Viral contig:  $VIRAL_CONTIG"
echo "  Reference:     $REFERENCE"
echo "  Annotation:    $ANNOTATION"
echo "  Input BAMs:    $INPUT_DIR"
echo "  Results dir:   $RESULTS_DIR"
echo "═══════════════════════════════════════════════════════════"
