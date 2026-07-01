#!/bin/bash
# ============================================================
# run_viloca.sh
# Local haplotype reconstruction with VILOCA.
#
# Per sample it:
#   1. Filters target_only.bam to primary alignments only
#   2. Runs VILOCA in uniform-window mode
#
# Usage:
#   cd Variant_discovery_pipeline
#   ../Haplotype/scripts/run_viloca.sh                       # all samples
#   ../Haplotype/scripts/run_viloca.sh INH_3_DPI_R1_A3       # one sample
# ============================================================

# ── Settings ────────────────────────────────────────────────
WINDOW_SIZE=141             # ~read length; must be divisible by 3 (win_shifts)
MODE="use_quality_scores"   # recommended: uses Phred scores from BAM
WIN_MIN_EXT=0.85            # min fraction of window a read must cover
NON_VAR_POS_THRESHOLD=0.005  # low-frequency noise cutoff (0.5%)
THREADS=24                  # parallelism
MIN_PRIMARY_READS=10000     # skip samples with too few usable reads
MIN_WINDOWS_COVERAGE=10     # omit windows with low depth for stability
# ────────────────────────────────────────────────────────────

# ── Paths ───────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../Variant_discovery_pipeline" && pwd)"
OUTDIR="$(cd "$SCRIPT_DIR/.." && pwd)/viloca"
LOGDIR="$OUTDIR/logs"

REF="$PIPELINE_DIR/Input/Reference/target_only.fasta"
LOFREQ_DIR="$PIPELINE_DIR/LoFreq"

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

LOG_FILE="$LOGDIR/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ── Conda ───────────────────────────────────────────────────
eval "$(conda shell.bash hook)"
conda activate env_viloca 2>/dev/null \
  || { echo "ERROR: could not activate env_viloca"; exit 1; }

# Also need samtools — check
for tool in viloca samtools; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not found"; exit 1; }
done

set -euo pipefail

# ── Decide which samples to process ────────────────────────
if [ $# -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for d in "$LOFREQ_DIR"/*/; do
        SAMPLES+=("$(basename "$d")")
    done
fi

echo "================================================================"
echo "  VILOCA haplotype reconstruction"
echo "  Window: ${WINDOW_SIZE} bp | Mode: ${MODE}"
echo "  Threads: ${THREADS}"
echo "  Noise cutoff: ${NON_VAR_POS_THRESHOLD}"
echo "  Min primary reads: ${MIN_PRIMARY_READS}"
echo "  Min window coverage: ${MIN_WINDOWS_COVERAGE}"
echo "  Samples: ${#SAMPLES[@]}"
echo "  Output : $OUTDIR"
echo "  Log    : $LOG_FILE"
echo "================================================================"
echo ""

SUCCESS_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0
FAILED_SAMPLES=()

# ── Process each sample ────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    echo "── $SAMPLE ──────────────────────────────────"

    SRC_BAM="$LOFREQ_DIR/$SAMPLE/target_only.bam"
    SAMPLE_OUT="$OUTDIR/$SAMPLE"

    # --- check input ---
    if [ ! -f "$SRC_BAM" ]; then
        echo "  SKIP: no target_only.bam found"
        echo ""
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    mkdir -p "$SAMPLE_OUT"

    # 1. Filter to primary alignments only + single-contig header
    #    VILOCA requires: primary-only reads AND single reference in BAM header
    #    -F 0x904 = exclude unmapped, secondary, supplementary
    PRIMARY_BAM="$SAMPLE_OUT/primary_only.bam"
    echo "  1) Filtering to primary alignments (single-contig BAM) ..."

    # Build a minimal header with only the Target contig
    TARGET_CONTIG="${TARGET_CONTIG:-}"
    {
        samtools view -H "$SRC_BAM" | grep "^@HD"
        samtools view -H "$SRC_BAM" | grep "^@SQ" | grep "SN:${TARGET_CONTIG}"
        samtools view -H "$SRC_BAM" | grep "^@RG"
        samtools view -H "$SRC_BAM" | grep "^@PG"
    } > "$SAMPLE_OUT/target_header.sam"

    samtools view -b -F 0x904 "$SRC_BAM" "$TARGET_CONTIG" \
        | samtools reheader "$SAMPLE_OUT/target_header.sam" - \
        | samtools sort -o "$PRIMARY_BAM" -
    samtools index "$PRIMARY_BAM"

    NREADS=$(samtools view -c "$PRIMARY_BAM")
    echo "     → $NREADS primary reads"

    if [ "$NREADS" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few primary reads (${NREADS} < ${MIN_PRIMARY_READS})"
        echo ""
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # 2. Run VILOCA
    echo "  2) Running VILOCA (window=${WINDOW_SIZE}, mode=${MODE}) ..."

    # Remove stale run artifacts from previous interrupted/failed runs
    # that can cause shutil.move collisions (e.g. work/inference, work/raw_reads)
    rm -rf "$SAMPLE_OUT/inference" "$SAMPLE_OUT/raw_reads" "$SAMPLE_OUT/work"

    pushd "$SAMPLE_OUT" > /dev/null

    if viloca run \
        -b "$PRIMARY_BAM" \
        -f "$REF" \
        -w "$WINDOW_SIZE" \
        -t "$THREADS" \
        --mode "$MODE" \
        --win_min_ext "$WIN_MIN_EXT" \
        --min_windows_coverage "$MIN_WINDOWS_COVERAGE" \
        --exclude_non_var_pos_threshold "$NON_VAR_POS_THRESHOLD"; then
        popd > /dev/null
        echo "  ✓  Done: $SAMPLE_OUT"
        echo ""
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        popd > /dev/null
        echo "  ✗  FAILED: $SAMPLE (see $LOG_FILE)"
        echo ""
        FAIL_COUNT=$((FAIL_COUNT + 1))
        FAILED_SAMPLES+=("$SAMPLE")
        continue
    fi
done

echo "================================================================"
echo "Done! Haplotype results are in: $OUTDIR"
echo "Summary: success=$SUCCESS_COUNT, skipped=$SKIP_COUNT, failed=$FAIL_COUNT"
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "Failed samples: ${FAILED_SAMPLES[*]}"
fi
echo ""
echo "Key outputs per sample:"
echo "  haplotypes/                → local haplotype FASTAs"
echo "  cooccurring_mutations.csv  → mutation linkage"
echo "  coverage.txt               → per-window read counts"
echo "================================================================"

