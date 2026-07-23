#!/bin/bash
# =====================================================================
# run_ivar.sh — iVar variant calling (adapted for nf-core RNAseq BAMs)
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Uses DATASET-aware config.sh (VIRAL_CONTIG per dataset)
#   - No primer trimming (RNA-seq, primer_bed is always empty)
#   - Results in results/<dataset>/Ivar/
#
# Usage:
#   DATASET=mouse_veev bash run_ivar.sh
#   DATASET=mouse_veev bash run_ivar.sh s101   # single sample
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Setup logging
mkdir -p "$IVAR_OUT"
LOG_FILE="$IVAR_OUT/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Activate conda
eval "$(conda shell.bash hook)"
conda activate ivar_env

echo "Starting iVar pipeline for: $DATASET"

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
    local out_dir="$IVAR_OUT/${sample_name}"

    # Skip if already processed and FORCE_RECALL is false
    FORCE_RECALL="${FORCE_RECALL:-false}"
    if [ "$FORCE_RECALL" = "false" ] && [ -f "$out_dir/variants.tsv" ] && [ -f "$out_dir/consensus.fa" ]; then
        echo "  ✓ $sample_name: variants already called (iVar), skipping"
        return 0
    fi

    mkdir -p "$out_dir"

    echo "── Processing: $sample_name ──"

    # Step 1: Extract viral reads only
    echo "  Extracting viral reads from contig: $VIRAL_CONTIG"
    samtools view -@ 2 -b "$bam_file" "$VIRAL_CONTIG" | \
        samtools sort -@ 2 -o "$out_dir/viral_only.bam" -
    samtools index "$out_dir/viral_only.bam"
    local viral_bam="$out_dir/viral_only.bam"

    # Step 2: No primer trimming (RNA-seq)
    # (Original workflow trims primers here; skipped for RNA-seq)
    local input_bam="$viral_bam"

    # Step 3: Generate consensus
    echo "  Creating consensus..."
    samtools mpileup -aa -A -d 0 -Q 0 -r "$VIRAL_CONTIG" \
        --reference "$REFERENCE" "$input_bam" | \
        ivar consensus -p "$out_dir/consensus" \
            -m "$IVAR_MIN_CONSENSUS_COVERAGE" -t 0.5 -q "$IVAR_MIN_BASE_QUALITY"

    # Step 4: Call variants
    echo "  Calling variants..."
    samtools mpileup -aa -A -d 0 -Q 0 -r "$VIRAL_CONTIG" \
        --reference "$REFERENCE" "$input_bam" | \
        ivar variants -p "$out_dir/variants" \
            -r "$REFERENCE" \
            -m "$IVAR_MIN_VARIANT_DEPTH" \
            -t "$IVAR_MIN_VARIANT_FREQ"

    # Step 5: QC stats
    echo "  Generating QC stats..."
    local viral_reads
    viral_reads=$(samtools view -c "$out_dir/viral_only.bam")
    local variant_count=0
    if [ -f "$out_dir/variants.tsv" ]; then
        variant_count=$(tail -n +2 "$out_dir/variants.tsv" | wc -l)
    fi

    {
        echo "Sample: $sample_name"
        echo "Viral reads: $viral_reads"
        echo "Variants: $variant_count"
        echo "Reference: $REFERENCE"
        echo "Contig: $VIRAL_CONTIG"
    } > "$out_dir/qc_stats.txt"

    echo "  ✓ Complete: $sample_name ($variant_count variants)"
    echo ""
}

# Process single sample or all
if [ -n "${1:-}" ]; then
    process_sample "$INPUT_DIR/${1}.bam"
else
    # All samples
    FORCE_RECALL="${FORCE_RECALL:-false}"
    if [ "$FORCE_RECALL" = "true" ]; then
        echo "FORCE_RECALL=true: Cleaning previous iVar results..."
        rm -rf "$IVAR_OUT"/*
    else
        echo "FORCE_RECALL=false: Keeping previous iVar results, skipping already processed samples."
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

    if [ "${#pids[@]}" -gt 0 ]; then
        wait "${pids[@]}"
    fi
fi

echo "═══════════════════════════════════════════════════════════"
echo "  iVar pipeline complete for: $DATASET"
echo "  Results: $IVAR_OUT"
echo "═══════════════════════════════════════════════════════════"
