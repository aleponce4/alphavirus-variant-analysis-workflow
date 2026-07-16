#!/bin/bash
# =====================================================================
# run_lofreq.sh — LoFreq variant calling (adapted for nf-core RNAseq BAMs)
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Uses DATASET-aware config.sh (VIRAL_CONTIG per dataset)
#   - Reads from input/<dataset>/BAMs/ (symlinked star_salmon BAMs)
#   - No primer handling (RNA-seq)
#   - Results in results/<dataset>/LoFreq/
#
# Usage:
#   DATASET=mouse_veev bash run_lofreq.sh
#   DATASET=mouse_veev bash run_lofreq.sh s101   # single sample
# =====================================================================
set -uo pipefail
# Note: no 'set -e' — background jobs may return non-zero, which would kill the script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Setup logging
mkdir -p "$LOFREQ_OUT"
LOG_FILE="$LOFREQ_OUT/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Activate conda
eval "$(conda shell.bash hook)"
conda activate lofreq-env

echo "Starting LoFreq pipeline for: $DATASET"

# Index reference
echo "Indexing reference file..."
samtools faidx "$REFERENCE"

# Function to process a single sample
process_sample() {
    local bam_file="$1"
    if [ ! -f "$bam_file" ]; then
        echo "BAM file not found: $bam_file"
        return 1
    fi

    local sample_name
    sample_name=$(basename "$bam_file" .bam)
    local out_dir="$LOFREQ_OUT/${sample_name}"

    # Skip if already processed and FORCE_RECALL is false
    FORCE_RECALL="${FORCE_RECALL:-false}"
    if [ "$FORCE_RECALL" = "false" ] && [ -f "$out_dir/variants.filtered.vcf.gz" ] && [ -f "$out_dir/variants.filtered.vcf.gz.tbi" ]; then
        echo "  ✓ $sample_name: variants already called (LoFreq), skipping"
        return 0
    fi

    mkdir -p "$out_dir"

    echo "── Processing: $sample_name ──"

    # Step 1: Extract viral reads only
    # Use pre-extracted viral BAM if available (from extract_viral_bams.sh)
    local preextracted="${VARIANT_ROOT}/input/${DATASET}/viral_bams/${sample_name}.viral_only.bam"
    local no_virus_marker="${VARIANT_ROOT}/input/${DATASET}/viral_bams/${sample_name}.no_virus"

    if [ -f "$no_virus_marker" ]; then
        echo "  ⊘ No viral reads (pre-screened), skipping"
        echo "Sample: $sample_name" > "$out_dir/qc_stats.txt"
        echo "Viral reads: 0" >> "$out_dir/qc_stats.txt"
        echo "Variants (filtered): 0" >> "$out_dir/qc_stats.txt"
        return 0
    elif [ -f "$preextracted" ]; then
        echo "  Using pre-extracted viral BAM (SSD)"
        cp -f "$preextracted" "$out_dir/viral_only.bam"
        cp -f "${preextracted}.bai" "$out_dir/viral_only.bam.bai" 2>/dev/null || samtools index "$out_dir/viral_only.bam"
    else
        echo "  Extracting viral reads from contig: $VIRAL_CONTIG"
        samtools view -b "$bam_file" "$VIRAL_CONTIG" | \
            samtools sort -o "$out_dir/viral_only.bam" -
        samtools index "$out_dir/viral_only.bam"
    fi
    local input_bam="$out_dir/viral_only.bam"

    # Step 2: Viterbi realignment (BAQ)
    if [ "$LOFREQ_BAQ" -eq 1 ]; then
        echo "  Viterbi realignment..."
        lofreq viterbi -f "$REFERENCE" "$input_bam" | \
            samtools sort -o "$out_dir/realigned.bam" -
        samtools index "$out_dir/realigned.bam"
        input_bam="$out_dir/realigned.bam"
    fi

    # Step 3: Add indel qualities
    if [ "$LOFREQ_ENABLE_INDELQUAL" -eq 1 ]; then
        echo "  Adding indel qualities..."
        lofreq indelqual --dindel -f "$REFERENCE" \
            -o "$out_dir/indelqual.bam" "$input_bam"
        samtools index "$out_dir/indelqual.bam"
        input_bam="$out_dir/indelqual.bam"
    fi

    # Step 4: Call variants
    echo "  Calling variants..."
    lofreq call-parallel --pp-threads "$THREADS" -f "$REFERENCE" \
        --min-cov "$LOFREQ_MIN_VARIANT_DEPTH" \
        --min-bq "$LOFREQ_MIN_BASE_QUALITY" \
        --min-alt-bq "$LOFREQ_MIN_BASE_QUALITY" \
        --min-mq "$LOFREQ_MIN_MAP_QUALITY" \
        --sig "$LOFREQ_CALL_SIG" \
        -o "$out_dir/variants.vcf" \
        "$input_bam"

    # Step 5: Filter variants
    echo "  Filtering variants..."
    lofreq filter -i "$out_dir/variants.vcf" \
        -o "$out_dir/variants.filtered.vcf" \
        --snvqual-thresh 20 --indelqual-thresh 20

    # Step 6: Compress and index
    bgzip -f "$out_dir/variants.filtered.vcf"
    tabix -f -p vcf "$out_dir/variants.filtered.vcf.gz"

    # Step 7: QC stats
    echo "  Generating QC stats..."
    local viral_reads
    viral_reads=$(samtools view -c "$out_dir/viral_only.bam")
    local variant_count
    variant_count=$(zcat "$out_dir/variants.filtered.vcf.gz" | grep -v "^#" | wc -l)

    {
        echo "Sample: $sample_name"
        echo "Viral reads: $viral_reads"
        echo "Variants (filtered): $variant_count"
        echo "Reference: $REFERENCE"
        echo "Contig: $VIRAL_CONTIG"
    } > "$out_dir/qc_stats.txt"

    echo "  ✓ Complete: $sample_name ($variant_count variants)"
    echo ""
}

# Process single sample or all
if [ -n "${1:-}" ]; then
    # Single sample
    process_sample "$INPUT_DIR/${1}.bam"
else
    # All samples
    FORCE_RECALL="${FORCE_RECALL:-false}"
    if [ "$FORCE_RECALL" = "true" ]; then
        echo "FORCE_RECALL=true: Cleaning previous LoFreq results..."
        rm -rf "$LOFREQ_OUT"/*
    else
        echo "FORCE_RECALL=false: Keeping previous LoFreq results, skipping already processed samples."
    fi

    job_count=0
    pids=()

    for bam_file in "$INPUT_DIR"/*.bam; do
        if [ ! -f "$bam_file" ]; then
            echo "No BAM files found in $INPUT_DIR"
            break
        fi

        process_sample "$bam_file" &
        pids+=($!)
        job_count=$((job_count + 1))

        if [ "$job_count" -ge "$MAX_JOBS" ]; then
            echo "Reached max concurrent jobs ($MAX_JOBS). Waiting..."
            wait "${pids[@]}"
            pids=()
            job_count=0
        fi
    done

    # Wait for remaining jobs
    if [ "${#pids[@]}" -gt 0 ]; then
        wait "${pids[@]}"
    fi
fi

echo "═══════════════════════════════════════════════════════════"
echo "  LoFreq pipeline complete for: $DATASET"
echo "  Results: $LOFREQ_OUT"
echo "═══════════════════════════════════════════════════════════"
