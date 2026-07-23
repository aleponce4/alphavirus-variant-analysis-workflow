#!/bin/bash
# =====================================================================
# run_cliquesnv.sh — Run CliqueSNV haplotype reconstruction (dataset-aware)
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Uses dataset-aware config.sh paths and variables
#   - Processes BAMs from results/$DATASET/LoFreq/
#   - Saves all intermediate and final outputs inside results/$DATASET/CliqueSNV
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Paths
CLIQUESNV_OUT="${RESULTS_DIR}/CliqueSNV"
LOGDIR="${CLIQUESNV_OUT}/logs"
ANALYSIS_DIR="${CLIQUESNV_OUT}/Analysis"

mkdir -p "$CLIQUESNV_OUT" "$LOGDIR" "$ANALYSIS_DIR"

LOG_FILE="$LOGDIR/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Settings
CLIQUESNV_ENV="env_cliquesnv"
MIN_PRIMARY_READS=10000
TIME_LIMIT_SECONDS=10800
TF_VALUES=(0.01)
CLIQUESNV_MODE="snv-illumina"
CLIQUESNV_T=10
CLIQUESNV_MIN_MAPQ=30
CLIQUESNV_MAX_NM=5
THREADS_CLIQUESNV="${THREADS:-4}"


# Conda activation
eval "$(conda shell.bash hook)"
if conda activate "$CLIQUESNV_ENV" 2>/dev/null; then
    echo "Activated conda environment: $CLIQUESNV_ENV"
else
    echo "ERROR: Conda environment '$CLIQUESNV_ENV' not found!"
    echo "  Please create it first by running:"
    echo "  wsl -d Ubuntu -- bash -lc \"conda env create -f ${VARIANT_ROOT}/workflow/original_repo/envs/env_cliquesnv.yml\""
    exit 1
fi

CLIQUESNV_BIN="${CONDA_PREFIX}/bin/cliquesnv"
SAMTOOLS_BIN="/usr/bin/samtools"
PYTHON_BIN="${CONDA_PREFIX}/bin/python"

# Function to process one sample
process_sample() {
    local sample="$1"
    local status_file="$LOGDIR/${sample}.status"
    local sample_out="${CLIQUESNV_OUT}/${sample}"

    echo "── Processing Haplotypes: $sample ──"

    # Input BAM is from LoFreq output folder
    local src_bam="${LOFREQ_OUT}/${sample}/viral_only.bam"

    if [ ! -f "$src_bam" ]; then
        echo "  SKIP: no viral_only.bam found for $sample"
        printf "SKIP\tno_viral_only_bam\n" > "$status_file"
        return 0
    fi

    # Check if the input BAM has at least some alignments to avoid samtools sort crashes
    local raw_read_count
    raw_read_count=$("$SAMTOOLS_BIN" view -c "$src_bam" 2>/dev/null || echo 0)
    if [ "$raw_read_count" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few raw reads ($raw_read_count < $MIN_PRIMARY_READS)"
        printf "SKIP\ttoo_few_raw_reads\n" > "$status_file"
        return 0
    fi

    # Check if we already processed it
    if [ "${FORCE_RECALL:-false}" = "false" ] && [ -f "$sample_out/tf_0p01/primary_only.fasta" ]; then
        echo "  ✓ $sample: CliqueSNV already run, skipping"
        printf "SUCCESS\tok\n" > "$status_file"
        return 0
    fi

    mkdir -p "$sample_out"

    local primary_bam="$sample_out/primary_only.bam"
    local primary_sam="$sample_out/primary_only.sam"

    # Filter to primary alignments
    echo "  1) Filtering to high-quality primary alignments ..."
    "$SAMTOOLS_BIN" view -h -F 0x904 -q "$CLIQUESNV_MIN_MAPQ" -e "[NM] <= $CLIQUESNV_MAX_NM" "$src_bam" "$VIRAL_CONTIG" \
        | grep -E "^@HD|^@RG|^@PG|^@CO|${VIRAL_CONTIG}" \
        | "$SAMTOOLS_BIN" view -b - \
        | "$SAMTOOLS_BIN" sort -o "$primary_bam" -
    "$SAMTOOLS_BIN" index "$primary_bam"

    local nreads
    nreads=$("$SAMTOOLS_BIN" view -c "$primary_bam")
    echo "     → $nreads primary reads"

    if [ "$nreads" -lt "$MIN_PRIMARY_READS" ]; then
        echo "  SKIP: too few primary reads ($nreads < $MIN_PRIMARY_READS)"
        printf "SKIP\ttoo_few_primary_reads\n" > "$status_file"
        return 0
    fi

    echo "  2) Converting BAM to SAM for CliqueSNV ..."
    "$SAMTOOLS_BIN" view -h "$primary_bam" > "$primary_sam"

    # Run for each threshold (tf=0.01 default)
    local sample_failed=0
    for tf in "${TF_VALUES[@]}"; do
        local tf_dir_name="tf_$(echo $tf | sed 's/\./p/')"
        local run_out="$sample_out/$tf_dir_name"
        mkdir -p "$run_out"

        echo "  3) Running CliqueSNV (tf=${tf}) ..."
        set +e
        "$CLIQUESNV_BIN" \
            -Xmx8g \
            -m "$CLIQUESNV_MODE" \
            -in "$primary_sam" \
            -outDir "$run_out" \
            -threads "$THREADS_CLIQUESNV" \
            -t "$CLIQUESNV_T" \
            -tf "$tf" \
            -tl "$TIME_LIMIT_SECONDS" \
            -log
        code=$?
        set -e

        if [ $code -eq 0 ] && [ -d "$run_out/SNPGenie_Results" ] || [ -f "$run_out/cliquesnv.fasta" ] || [ -n "$(ls -A "$run_out" 2>/dev/null)" ]; then
            echo "     ✓ tf=${tf} done"
        else
            echo "     ✗ tf=${tf} failed"
            sample_failed=1
        fi
    done

    # Clean intermediate SAM to save space
    rm -f "$primary_sam"

    if [ "$sample_failed" -eq 1 ]; then
        printf "FAIL\tcliquesnv_threshold_failed\n" > "$status_file"
        echo "  ✗ FAILED: $sample"
        return 1
    fi

    printf "SUCCESS\tok\n" > "$status_file"
    echo "  ✓ Done: $sample"
    return 0
}

# Run all samples in parallel using the MAX_JOBS pool
echo "Starting CliqueSNV runs in parallel (MAX_JOBS=$MAX_JOBS)..."
for sample_dir in "$LOFREQ_OUT"/*/; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    sample=$(basename "$sample_dir")
    
    # Check if sample had 0 reads (pre-screened)
    no_virus_marker="${VARIANT_ROOT}/input/${DATASET}/viral_bams/${sample}.no_virus"
    if [ -f "$no_virus_marker" ]; then
        echo "  ⊘ $sample: 0 viral reads (pre-screened), skipping Haplotypes"
        continue
    fi

    # Skip check outside background job to avoid slot wastage
    sample_out="$CLIQUESNV_OUT/$sample"
    if [ "${FORCE_RECALL:-false}" = "false" ] && [ -f "$sample_out/tf_0p01/primary_only.fasta" ]; then
        echo "  ✓ $sample: CliqueSNV already run, skipping"
        continue
    fi

    process_sample "$sample" &

    # Queue control: wait if we have MAX_JOBS running
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
done

# Wait for all remaining background jobs
wait

# Summarize results
echo ""
echo "Generating Haplotype summary tables..."
"$PYTHON_BIN" "${VARIANT_ROOT}/workflow/original_repo/Haplotype/scripts/summarize_cliquesnv_brief.py" \
    --cliquesnv-dir "$CLIQUESNV_OUT" \
    --out-dir "$ANALYSIS_DIR"

"$PYTHON_BIN" "${VARIANT_ROOT}/workflow/original_repo/Haplotype/scripts/build_cliquesnv_haplotype_tables.py" \
    --cliquesnv-dir "$CLIQUESNV_OUT" \
    --out-dir "$ANALYSIS_DIR" \
    --threshold-label "tf_0p01" \
    --reference-fasta "$REFERENCE" \
    --annotation-gff "$ANNOTATION"

echo "═══════════════════════════════════════════════════════════"
echo "  CliqueSNV haplotype analysis complete"
echo "  Results: $CLIQUESNV_OUT"
echo "═══════════════════════════════════════════════════════════"
