#!/bin/bash
# =====================================================================
# extract_viral_bams.sh — Pre-extract viral contig from all BAMs to SSD
#
# This is the critical optimization: reading 3.4 GB BAMs from HDD via
# WSL's /mnt/d (9p filesystem) is very slow. We extract the viral-only
# BAM (~250 MB) once, then all downstream steps read from SSD.
#
# Usage:
#   DATASET=mouse_veev bash extract_viral_bams.sh
# =====================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

# Output directory on WSL ext4 (SSD-backed)
VIRAL_BAM_DIR="${VARIANT_ROOT}/input/${DATASET}/viral_bams"
mkdir -p "$VIRAL_BAM_DIR"

echo "═══════════════════════════════════════════════════════════"
echo "  Pre-extracting viral BAMs: $DATASET"
echo "  Contig: $VIRAL_CONTIG"
echo "  Source: $INPUT_DIR (HDD via /mnt/d)"
echo "  Target: $VIRAL_BAM_DIR (SSD)"
echo "═══════════════════════════════════════════════════════════"

# Function to extract one sample
extract_sample() {
    local bam_file="$1"
    local sample_name
    sample_name=$(basename "$bam_file" .bam)
    local out_bam="$VIRAL_BAM_DIR/${sample_name}.viral_only.bam"

    # Skip if already exists and is non-empty
    if [ -f "$out_bam" ] && [ "$(stat -c%s "$out_bam" 2>/dev/null || echo 0)" -gt 1000 ]; then
        echo "  ✓ $sample_name: already extracted, skipping"
        return 0
    fi

    # Check viral reads via idxstats (fast, reads index only)
    local viral_reads
    viral_reads=$(samtools idxstats "$bam_file" 2>/dev/null | grep "^${VIRAL_CONTIG}" | awk '{print $3}')
    viral_reads=${viral_reads:-0}

    if [ "$viral_reads" -eq 0 ]; then
        echo "  ⊘ $sample_name: 0 viral reads (Mock?), skipping"
        # Create empty marker file
        echo "0" > "$VIRAL_BAM_DIR/${sample_name}.no_virus"
        return 0
    fi

    echo "  → $sample_name: $viral_reads viral reads, extracting..."
    samtools view -b -@ 2 "$bam_file" "$VIRAL_CONTIG" | \
        samtools sort -@ 2 -o "$out_bam" -
    samtools index "$out_bam"
    echo "  ✓ $sample_name: done ($(stat -c%s "$out_bam" | numfmt --to=iec))"
}

export -f extract_sample
export VIRAL_BAM_DIR VIRAL_CONTIG INPUT_DIR

# Process all BAMs in parallel (8 jobs × 2 threads = 16 cores)
MAX_EXTRACT_JOBS=8
job_count=0
pids=()

for bam_file in "$INPUT_DIR"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in $INPUT_DIR"
        break
    fi

    extract_sample "$bam_file" &
    pids+=($!)
    job_count=$((job_count + 1))

    if [ "$job_count" -ge "$MAX_EXTRACT_JOBS" ]; then
        echo "  Waiting for $MAX_EXTRACT_JOBS jobs to finish..."
        wait "${pids[@]}"
        pids=()
        job_count=0
    fi
done

if [ "${#pids[@]}" -gt 0 ]; then
    wait "${pids[@]}"
fi

# Summary
echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  Extraction complete"
echo "═══════════════════════════════════════════════════════════"
extracted=$(ls "$VIRAL_BAM_DIR"/*.viral_only.bam 2>/dev/null | wc -l)
skipped=$(ls "$VIRAL_BAM_DIR"/*.no_virus 2>/dev/null | wc -l)
total_size=$(du -sh "$VIRAL_BAM_DIR" 2>/dev/null | awk '{print $1}')
echo "  Extracted: $extracted samples"
echo "  Skipped (no virus): $skipped samples"
echo "  Total size: $total_size"
echo ""
echo "  Viral BAMs ready at: $VIRAL_BAM_DIR"
echo "═══════════════════════════════════════════════════════════"
