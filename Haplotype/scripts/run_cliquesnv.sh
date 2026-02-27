#!/bin/bash
# ============================================================
# run_cliquesnv.sh
# Local haplotype reconstruction with CliqueSNV (Illumina mode).
#
# Per sample it:
#   1. Filters viral_only.bam to primary alignments only
#   2. Converts BAM -> SAM (CliqueSNV Illumina expects SAM input)
#   3. Runs CliqueSNV for multiple tf thresholds
#   4. Writes brief per-sample summary table
#
# Usage:
#   cd Variant_discovery_pipeline
#   ../Haplotype/scripts/run_cliquesnv.sh                       # all samples
#   ../Haplotype/scripts/run_cliquesnv.sh INH_3_DPI_R1_A3       # one sample
# ============================================================

# ── Settings ────────────────────────────────────────────────
THREADS="${THREADS:-0}"
HEAP_SIZE_GB="${HEAP_SIZE_GB:-0}"
MAX_PARALLEL_JOBS="${MAX_PARALLEL_JOBS:-0}"
TARGET_THREADS_PER_JOB="${TARGET_THREADS_PER_JOB:-16}"
TARGET_HEAP_GB_PER_JOB="${TARGET_HEAP_GB_PER_JOB:-40}"
RESERVE_MEM_GB="${RESERVE_MEM_GB:-8}"

MIN_PRIMARY_READS=10000
TIME_LIMIT_SECONDS=10800
TF_VALUES=(0.01)
CLIQUESNV_TF_VALUES="${CLIQUESNV_TF_VALUES:-}"
HAP_TABLE_THRESHOLD_LABEL="${HAP_TABLE_THRESHOLD_LABEL:-}"

if [ -n "$CLIQUESNV_TF_VALUES" ]; then
    TF_VALUES=()
    for token in ${CLIQUESNV_TF_VALUES//,/ }; do
        [ -n "$token" ] && TF_VALUES+=("$token")
    done
fi

if [ ${#TF_VALUES[@]} -eq 0 ]; then
    echo "ERROR: no tf thresholds configured (TF_VALUES)"
    exit 1
fi

if [ -z "$HAP_TABLE_THRESHOLD_LABEL" ]; then
    if [ ${#TF_VALUES[@]} -eq 1 ]; then
        HAP_TABLE_THRESHOLD_LABEL="tf_${TF_VALUES[0]//./p}"
    else
        HAP_TABLE_THRESHOLD_LABEL="all"
    fi
fi

CLIQUESNV_MODE="snv-illumina"
CLIQUESNV_T="${CLIQUESNV_T:-10}"
CLIQUESNV_MIN_MAPQ="${CLIQUESNV_MIN_MAPQ:-30}"
CLIQUESNV_MAX_NM="${CLIQUESNV_MAX_NM:-5}"
CLIQUESNV_ENV="${CLIQUESNV_ENV:-env_cliquesnv}"
VIRAL_CONTIG="${VIRAL_CONTIG:-VEEV_INH}"
# ────────────────────────────────────────────────────────────

# ── Paths ───────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../Variant_discovery_pipeline" && pwd)"
OUTDIR="$(cd "$SCRIPT_DIR/.." && pwd)/cliquesnv"
LOGDIR="$OUTDIR/logs"
ANALYSIS_DIR="$OUTDIR/Analysis"

LOFREQ_DIR="$PIPELINE_DIR/LoFreq"
SUMMARY_SCRIPT="$SCRIPT_DIR/summarize_cliquesnv_brief.py"
HAP_TABLE_SCRIPT="$SCRIPT_DIR/build_cliquesnv_haplotype_tables.py"

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"
mkdir -p "$ANALYSIS_DIR"

LOG_FILE="$LOGDIR/run_$(date +%Y%m%d_%H%M%S).log"
if [ -t 1 ]; then
    exec > >(tee -a "$LOG_FILE") 2>&1
else
    exec >> "$LOG_FILE" 2>&1
fi

# ── Conda ──────────────────────────────────────────────────
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$CLIQUESNV_ENV" 2>/dev/null \
        || { echo "ERROR: could not activate $CLIQUESNV_ENV"; exit 1; }
    echo "Activated conda environment: $CLIQUESNV_ENV"
else
    echo "ERROR: conda not found (required for $CLIQUESNV_ENV)"
    exit 1
fi

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: CONDA_PREFIX is not set after activating $CLIQUESNV_ENV"
    exit 1
fi

CLIQUESNV_BIN="$CONDA_PREFIX/bin/cliquesnv"
SAMTOOLS_BIN="$CONDA_PREFIX/bin/samtools"
PYTHON_BIN="$CONDA_PREFIX/bin/python"

for tool in "$CLIQUESNV_BIN" "$SAMTOOLS_BIN" "$PYTHON_BIN"; do
    [ -x "$tool" ] || { echo "ERROR: required tool not found in env: $tool"; exit 1; }
done

"$PYTHON_BIN" -c "import pandas" >/dev/null 2>&1 || {
    echo "ERROR: pandas not available in $CLIQUESNV_ENV"
    echo "Install with: conda activate $CLIQUESNV_ENV && conda install -y pandas"
    exit 1
}

if [ ! -f "$SUMMARY_SCRIPT" ]; then
    echo "ERROR: summary script not found: $SUMMARY_SCRIPT"
    exit 1
fi

if [ ! -f "$HAP_TABLE_SCRIPT" ]; then
    echo "ERROR: haplotype table script not found: $HAP_TABLE_SCRIPT"
    exit 1
fi

set -euo pipefail

format_tf_dir() {
    local tf="$1"
    printf "tf_%s" "${tf//./p}"
}

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
echo "  CliqueSNV haplotype reconstruction"
echo "  Mode: ${CLIQUESNV_MODE}"
echo "  Time limit (s): ${TIME_LIMIT_SECONDS}"
echo "  tf thresholds: ${TF_VALUES[*]}"
echo "  Haplotype summary threshold label: ${HAP_TABLE_THRESHOLD_LABEL}"
echo "  Read filter: MAPQ>=${CLIQUESNV_MIN_MAPQ}, NM<=${CLIQUESNV_MAX_NM}"
echo "  Min primary reads: ${MIN_PRIMARY_READS}"
echo "  Samples: ${#SAMPLES[@]}"
echo "  Output : $OUTDIR"
echo "  Log    : $LOG_FILE"
echo "================================================================"
echo ""

CPU_TOTAL=$(nproc)
MEM_AVAILABLE_GB=$(awk '/MemAvailable:/ {printf "%d", $2/1024/1024}' /proc/meminfo)

SAMPLE_COUNT=${#SAMPLES[@]}
if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "No samples found. Exiting."
    exit 0
fi

if [ "$MAX_PARALLEL_JOBS" -gt 0 ]; then
    PARALLEL_JOBS="$MAX_PARALLEL_JOBS"
else
    JOBS_CPU=$(( CPU_TOTAL / TARGET_THREADS_PER_JOB ))
    if [ "$JOBS_CPU" -lt 1 ]; then
        JOBS_CPU=1
    fi

    MEM_FOR_JOBS_GB=$(( MEM_AVAILABLE_GB - RESERVE_MEM_GB ))
    if [ "$MEM_FOR_JOBS_GB" -lt "$TARGET_HEAP_GB_PER_JOB" ]; then
        JOBS_MEM=1
    else
        JOBS_MEM=$(( MEM_FOR_JOBS_GB / TARGET_HEAP_GB_PER_JOB ))
        if [ "$JOBS_MEM" -lt 1 ]; then
            JOBS_MEM=1
        fi
    fi

    PARALLEL_JOBS="$JOBS_CPU"
    if [ "$JOBS_MEM" -lt "$PARALLEL_JOBS" ]; then
        PARALLEL_JOBS="$JOBS_MEM"
    fi
fi

if [ "$PARALLEL_JOBS" -gt "$SAMPLE_COUNT" ]; then
    PARALLEL_JOBS="$SAMPLE_COUNT"
fi
if [ "$PARALLEL_JOBS" -lt 1 ]; then
    PARALLEL_JOBS=1
fi

if [ "$THREADS" -gt 0 ]; then
    THREADS_PER_JOB="$THREADS"
else
    THREADS_PER_JOB=$(( CPU_TOTAL / PARALLEL_JOBS ))
    if [ "$THREADS_PER_JOB" -lt 1 ]; then
        THREADS_PER_JOB=1
    fi
fi

if [ "$HEAP_SIZE_GB" -gt 0 ]; then
    HEAP_GB_PER_JOB="$HEAP_SIZE_GB"
else
    MEM_FOR_JOBS_GB=$(( MEM_AVAILABLE_GB - RESERVE_MEM_GB ))
    if [ "$MEM_FOR_JOBS_GB" -lt 8 ]; then
        MEM_FOR_JOBS_GB=8
    fi
    HEAP_GB_PER_JOB=$(( MEM_FOR_JOBS_GB / PARALLEL_JOBS ))
    if [ "$HEAP_GB_PER_JOB" -lt 8 ]; then
        HEAP_GB_PER_JOB=8
    fi
fi

HEAP_SIZE="${HEAP_GB_PER_JOB}g"

echo "Auto scheduling:"
echo "  CPU total: ${CPU_TOTAL}"
echo "  Mem available: ${MEM_AVAILABLE_GB} GiB"
echo "  Parallel jobs: ${PARALLEL_JOBS}"
echo "  Threads/job: ${THREADS_PER_JOB}"
echo "  Heap/job: ${HEAP_SIZE}"
echo ""

STATUS_DIR="$LOGDIR/status_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$STATUS_DIR"

process_sample() {
    local SAMPLE="$1"
    local STATUS_FILE="$STATUS_DIR/${SAMPLE}.status"

    echo "── $SAMPLE ──────────────────────────────────"

    local SRC_BAM="$LOFREQ_DIR/$SAMPLE/viral_only.bam"
    local SAMPLE_OUT="$OUTDIR/$SAMPLE"

    if [ ! -f "$SRC_BAM" ]; then
        echo "  SKIP: no viral_only.bam found"
        printf "SKIP\tno_viral_only_bam\n" > "$STATUS_FILE"
        return 0
    fi

    mkdir -p "$SAMPLE_OUT"

    local PRIMARY_BAM="$SAMPLE_OUT/primary_only.bam"
    local PRIMARY_SAM="$SAMPLE_OUT/primary_only.sam"
    local HEADER_SAM="$SAMPLE_OUT/viral_header.sam"

    echo "  1) Filtering to high-quality primary alignments ..."
    {
        "$SAMTOOLS_BIN" view -H "$SRC_BAM" | grep "^@HD"
        "$SAMTOOLS_BIN" view -H "$SRC_BAM" | grep "^@SQ" | grep "SN:${VIRAL_CONTIG}"
        "$SAMTOOLS_BIN" view -H "$SRC_BAM" | grep "^@RG" || true
        "$SAMTOOLS_BIN" view -H "$SRC_BAM" | grep "^@PG" || true
    } > "$HEADER_SAM"

    "$SAMTOOLS_BIN" view -b -F 0x904 -q "$CLIQUESNV_MIN_MAPQ" -e "[NM] <= $CLIQUESNV_MAX_NM" "$SRC_BAM" "$VIRAL_CONTIG" \
        | "$SAMTOOLS_BIN" reheader "$HEADER_SAM" - \
        | "$SAMTOOLS_BIN" sort -o "$PRIMARY_BAM" -
    "$SAMTOOLS_BIN" index "$PRIMARY_BAM"

    local NREADS
    NREADS=$("$SAMTOOLS_BIN" view -c "$PRIMARY_BAM")
    echo "     → $NREADS primary reads"

    if [ "$NREADS" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few primary reads (${NREADS} < ${MIN_PRIMARY_READS})"
        printf "SKIP\ttoo_few_primary_reads\n" > "$STATUS_FILE"
        return 0
    fi

    echo "  2) Converting BAM to SAM for CliqueSNV ..."
    "$SAMTOOLS_BIN" view -h "$PRIMARY_BAM" > "$PRIMARY_SAM"

    local AVG_DEPTH
    AVG_DEPTH=$("$SAMTOOLS_BIN" depth -aa -r "$VIRAL_CONTIG" "$PRIMARY_BAM" \
        | awk '{sum += $3; n += 1} END {if (n > 0) printf "%.2f", sum / n; else print "0"}')

    echo "  2b) Mean depth: ${AVG_DEPTH}x"
    echo "      CliqueSNV -t: ${CLIQUESNV_T} (fixed absolute floor) | threads=${THREADS_PER_JOB} | heap=${HEAP_SIZE}"

    local SAMPLE_FAILED=0
    local TF
    for TF in "${TF_VALUES[@]}"; do
        local TF_DIR_NAME
        TF_DIR_NAME="$(format_tf_dir "$TF")"
        local RUN_OUT="$SAMPLE_OUT/$TF_DIR_NAME"
        mkdir -p "$RUN_OUT"

        echo "  3) Running CliqueSNV (tf=${TF}) ..."
        if "$CLIQUESNV_BIN" \
            -Xmx"$HEAP_SIZE" \
            -m "$CLIQUESNV_MODE" \
            -in "$PRIMARY_SAM" \
            -outDir "$RUN_OUT" \
            -threads "$THREADS_PER_JOB" \
            -t "$CLIQUESNV_T" \
            -tf "$TF" \
            -tl "$TIME_LIMIT_SECONDS" \
            -log; then
            echo "     ✓ tf=${TF} done"
        else
            echo "     ✗ tf=${TF} failed"
            SAMPLE_FAILED=1
        fi
    done

    if [ "$SAMPLE_FAILED" -eq 1 ]; then
        printf "FAIL\tcliquesnv_threshold_failed\n" > "$STATUS_FILE"
        echo "  ✗ FAILED: $SAMPLE"
        return 0
    fi

    printf "SUCCESS\tok\n" > "$STATUS_FILE"
    echo "  ✓ Done: $SAMPLE_OUT"
    return 0
}

ACTIVE_JOBS=0
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_LOG="$LOGDIR/sample_${SAMPLE}.log"
    process_sample "$SAMPLE" > "$SAMPLE_LOG" 2>&1 &
    ACTIVE_JOBS=$((ACTIVE_JOBS + 1))

    if [ "$ACTIVE_JOBS" -ge "$PARALLEL_JOBS" ]; then
        wait -n || true
        ACTIVE_JOBS=$((ACTIVE_JOBS - 1))
    fi
done

wait || true

SUCCESS_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0
FAILED_SAMPLES=()

for STATUS_FILE in "$STATUS_DIR"/*.status; do
    [ -e "$STATUS_FILE" ] || continue
    SAMPLE_NAME="$(basename "$STATUS_FILE" .status)"
    STATUS_CODE="$(awk '{print $1}' "$STATUS_FILE")"
    case "$STATUS_CODE" in
        SUCCESS)
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            ;;
        SKIP)
            SKIP_COUNT=$((SKIP_COUNT + 1))
            ;;
        FAIL)
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAILED_SAMPLES+=("$SAMPLE_NAME")
            ;;
        *)
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAILED_SAMPLES+=("$SAMPLE_NAME")
            ;;
    esac
done

echo "================================================================"
echo "Generating brief summary table ..."
"$PYTHON_BIN" "$SUMMARY_SCRIPT" \
    --cliquesnv-dir "$OUTDIR" \
    --out-dir "$ANALYSIS_DIR"

echo "Generating haplotype frequency tables ..."
"$PYTHON_BIN" "$HAP_TABLE_SCRIPT" \
    --cliquesnv-dir "$OUTDIR" \
    --out-dir "$ANALYSIS_DIR" \
    --threshold-label "$HAP_TABLE_THRESHOLD_LABEL"
echo "================================================================"

echo "Done! CliqueSNV results are in: $OUTDIR"
echo "Summary: success=$SUCCESS_COUNT, skipped=$SKIP_COUNT, failed=$FAIL_COUNT"
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "Failed samples: ${FAILED_SAMPLES[*]}"
fi
echo ""
echo "Key outputs per sample:" 
echo "  primary_only.sam                         → CliqueSNV input"
echo "  tf_0p01/                                 → threshold-specific results"
echo "  Analysis/cliquesnv_brief_per_sample.tsv  → concise summary"
echo "  Analysis/cliquesnv_haplotype_*.csv       → haplotype summaries"
echo "================================================================"
