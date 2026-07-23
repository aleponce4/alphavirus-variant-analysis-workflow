#!/bin/bash
# =====================================================================
# run_viloca.sh — Run VILOCA haplotype reconstruction (dataset-aware)
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Uses dataset-aware config.sh paths and variables (dynamic VIRAL_CONTIG)
#   - Safe 3-step BAM extraction with verification and read-count filtering
#   - Controlled concurrency (MAX_JOBS=3) to prevent I/O and RAM thrashing
#   - Saves all intermediate and final outputs inside results/$DATASET/VILOCA
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Paths
VILOCA_OUT="${RESULTS_DIR}/VILOCA"
LOGDIR="${VILOCA_OUT}/logs"

mkdir -p "$VILOCA_OUT" "$LOGDIR"

LOG_FILE="$LOGDIR/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Settings
WINDOW_SIZE=141             # ~read length; must be divisible by 3 (win_shifts)
MODE="use_quality_scores"   # recommended: uses Phred scores from BAM
WIN_MIN_EXT=0.85            # min fraction of window a read must cover
NON_VAR_POS_THRESHOLD=0.005  # low-frequency noise cutoff (0.5%)
MIN_PRIMARY_READS=10000     # skip samples with too few usable reads
MAX_PRIMARY_READS=500000    # post-dedup safety cap (~5000x depth) to prevent VILOCA OOM
MIN_WINDOWS_COVERAGE=10     # omit windows with low depth for stability
THREADS_VILOCA="${THREADS:-4}"
MAX_VILOCA_JOBS="${MAX_JOBS_VILOCA:-1}" # Default 1 job at a time to prevent RAM exhaustion

# Conda activation
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
eval "$(conda shell.bash hook)"
if conda activate env_viloca 2>/dev/null; then
    echo "Activated conda environment: env_viloca"
else
    echo "ERROR: Conda environment 'env_viloca' not found!"
    exit 1
fi

# Check tools
for tool in viloca samtools; do
    command -v "$tool" >/dev/null || { echo "ERROR: required tool $tool not found in environment"; exit 1; }
done

process_sample() {
    local sample="$1"
    local status_file="${VILOCA_OUT}/${sample}/status.txt"
    mkdir -p "${VILOCA_OUT}/${sample}"

    echo "── Processing VILOCA: $sample ──"

    local SRC_BAM="${LOFREQ_OUT}/${sample}/viral_only.bam"
    local SAMPLE_OUT="${VILOCA_OUT}/${sample}"

    if [ ! -f "$SRC_BAM" ]; then
        echo "  SKIP: no viral_only.bam found"
        printf "SKIP\tno_viral_only_bam\n" > "$status_file"
        return 0
    fi

    # Pre-check raw read count
    local raw_read_count
    raw_read_count=$(samtools view -c "$SRC_BAM" 2>/dev/null || echo 0)
    if [ "$raw_read_count" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few raw reads ($raw_read_count < $MIN_PRIMARY_READS)"
        printf "SKIP\ttoo_few_raw_reads\n" > "$status_file"
        return 0
    fi

    local RAW_PRIMARY="$SAMPLE_OUT/raw_primary.bam"
    local PRIMARY_BAM="$SAMPLE_OUT/primary_only.bam"
    echo "  1) Filtering to primary alignments & deduplicating ..."

    # Safe 3-step BAM extraction to prevent pipe truncation
    local HEADER_SAM="$SAMPLE_OUT/header.sam"
    local READS_SAM="$SAMPLE_OUT/reads.sam"

    samtools view -H "$SRC_BAM" | grep -E "^@HD|^@RG|^@PG|^@CO|${VIRAL_CONTIG}" > "$HEADER_SAM"
    samtools view -F 0x904 -q 30 "$SRC_BAM" "$VIRAL_CONTIG" | awk '$6 !~ /N/' > "$READS_SAM"
    
    cat "$HEADER_SAM" "$READS_SAM" | samtools view -b - | samtools sort -o "$RAW_PRIMARY" -
    rm -f "$HEADER_SAM" "$READS_SAM"

    if [ ! -f "$RAW_PRIMARY" ]; then
        echo "  ✗  FAILED: $sample (BAM extraction failed)"
        printf "FAIL\tbam_creation_failed\n" > "$status_file"
        return 1
    fi

    # Collapse PCR duplicates (markdup)
    local DEDUP_BAM="$SAMPLE_OUT/dedup.bam"
    if samtools collate -O -u "$RAW_PRIMARY" 2>/dev/null | samtools fixmate -m -u - - 2>/dev/null | samtools sort -u - 2>/dev/null | samtools markdup -r - "$DEDUP_BAM" 2>/dev/null; then
        mv "$DEDUP_BAM" "$PRIMARY_BAM"
        rm -f "$RAW_PRIMARY"
    else
        mv "$RAW_PRIMARY" "$PRIMARY_BAM"
    fi

    local NREADS
    NREADS=$(samtools view -c "$PRIMARY_BAM" 2>/dev/null || echo 0)
    echo "     → $NREADS deduplicated primary reads"

    if [ "$NREADS" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few deduplicated primary reads (${NREADS} < ${MIN_PRIMARY_READS})"
        printf "SKIP\ttoo_few_primary_reads\n" > "$status_file"
        rm -f "$PRIMARY_BAM"
        return 0
    fi

    # Apply safety cap ONLY if post-dedup coverage is still ultra-high (> 500k reads) to prevent OOM
    if [ "$NREADS" -gt "$MAX_PRIMARY_READS" ]; then
        echo "  Notice: post-dedup reads (${NREADS}) exceed safety limit (${MAX_PRIMARY_READS}). Capping to prevent VILOCA OOM."
        local FRAC
        FRAC=$(awk -v r="$NREADS" -v max="$MAX_PRIMARY_READS" 'BEGIN {printf "%.4f", max/r}')
        local CAPPED_BAM="$SAMPLE_OUT/capped.bam"
        samtools view -s "42.${FRAC#0.}" -b "$PRIMARY_BAM" > "$CAPPED_BAM"
        mv "$CAPPED_BAM" "$PRIMARY_BAM"
        NREADS=$(samtools view -c "$PRIMARY_BAM" 2>/dev/null || echo 0)
        echo "     → Capped to $NREADS primary reads"
    fi

    samtools index "$PRIMARY_BAM"

    # Run VILOCA
    echo "  2) Running VILOCA (window=${WINDOW_SIZE}, mode=${MODE}) ..."
    rm -rf "$SAMPLE_OUT/inference" "$SAMPLE_OUT/raw_reads" "$SAMPLE_OUT/work"

    pushd "$SAMPLE_OUT" > /dev/null
    set +e
    viloca run \
        -b "$PRIMARY_BAM" \
        -f "$REFERENCE" \
        -w "$WINDOW_SIZE" \
        -t "$THREADS_VILOCA" \
        --mode "$MODE" \
        --win_min_ext "$WIN_MIN_EXT" \
        --min_windows_coverage "$MIN_WINDOWS_COVERAGE" \
        --exclude_non_var_pos_threshold "$NON_VAR_POS_THRESHOLD"
    local code=$?
    set -e
    popd > /dev/null

    # Clean intermediate files if successful
    if [ $code -eq 0 ]; then
        echo "  ✓  Done: $SAMPLE_OUT"
        printf "SUCCESS\tok\n" > "$status_file"
        return 0
    fi

    echo "  ✗  FAILED: $sample (exit code $code)"
    printf "FAIL\texit_code_${code}\n" > "$status_file"
    return 1
}

# Run all samples in parallel using the MAX_VILOCA_JOBS pool
echo "Starting VILOCA runs in parallel (MAX_JOBS=$MAX_VILOCA_JOBS)..."
for sample_dir in "$LOFREQ_OUT"/*/; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    sample=$(basename "$sample_dir")

    # Check if sample had 0 reads (pre-screened)
    no_virus_marker="${VARIANT_ROOT}/input/${DATASET}/viral_bams/${sample}.no_virus"
    if [ -f "$no_virus_marker" ]; then
        echo "  ⊘ $sample: 0 viral reads (pre-screened), skipping Haplotypes"
        mkdir -p "${VILOCA_OUT}/${sample}"
        printf "SKIP\tno_viral_reads\n" > "${VILOCA_OUT}/${sample}/status.txt"
        continue
    fi

    SAMPLE_OUT="${VILOCA_OUT}/${sample}"
    if [ "${FORCE_RECALL:-false}" = "false" ] && [ -f "$SAMPLE_OUT/cooccurring_mutations.csv" ]; then
        echo "  ✓ $sample: VILOCA already run, skipping"
        printf "SUCCESS\tok\n" > "$SAMPLE_OUT/status.txt"
        continue
    fi

    process_sample "$sample" &

    # Queue control: wait if we have MAX_VILOCA_JOBS running
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_VILOCA_JOBS" ]; do
        sleep 1
    done
done

# Wait for all remaining background jobs
wait

# Summarize statistics from status files
SUCCESS_COUNT=$(find "$VILOCA_OUT" -name "status.txt" -exec grep -q "^SUCCESS" {} \; -print 2>/dev/null | wc -l)
SKIP_COUNT=$(find "$VILOCA_OUT" -name "status.txt" -exec grep -q "^SKIP" {} \; -print 2>/dev/null | wc -l)
FAIL_COUNT=$(find "$VILOCA_OUT" -name "status.txt" -exec grep -q "^FAIL" {} \; -print 2>/dev/null | wc -l)

FAILED_SAMPLES=()
for f in $(find "$VILOCA_OUT" -name "status.txt" -exec grep -l "^FAIL" {} \; 2>/dev/null); do
    FAILED_SAMPLES+=($(basename "$(dirname "$f")"))
done

echo "================================================================"
echo "Done! VILOCA Haplotype results are in: $VILOCA_OUT"
echo "Summary: success=$SUCCESS_COUNT, skipped=$SKIP_COUNT, failed=$FAIL_COUNT"
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "Failed samples: ${FAILED_SAMPLES[*]}"
fi
echo "================================================================"
