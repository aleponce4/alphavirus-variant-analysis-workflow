#!/bin/bash
# =====================================================================
# calculate_coverage.sh — Assess viral coverage per sample
#
# Critical for RNA-seq data: coverage is uneven, so we need to
# determine appropriate depth thresholds before variant calling.
#
# Usage:
#   DATASET=mouse_veev bash calculate_coverage.sh
# =====================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

mkdir -p "$COVERAGE_OUT"

echo "═══════════════════════════════════════════════════════════"
echo "  Calculating viral coverage: $DATASET"
echo "  Viral contig: $VIRAL_CONTIG"
echo "═══════════════════════════════════════════════════════════"

# Summary file
SUMMARY="$COVERAGE_OUT/coverage_summary.tsv"
printf "sample\tmean_depth\tmin_depth\tmax_depth\tpercent_above_100x\tpercent_above_1000x\tpercent_above_5000x\n" > "$SUMMARY"

for bam_file in "$INPUT_DIR"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in $INPUT_DIR"
        break
    fi

    sample_name=$(basename "$bam_file" .bam)
    output_file="$COVERAGE_OUT/${sample_name}_coverage.txt"

    echo "Processing: $sample_name"

    # Per-position depth for viral contig only
    samtools depth -a -r "$VIRAL_CONTIG" "$bam_file" > "$output_file"

    # Calculate stats with awk
    awk -v sample="$sample_name" '
    BEGIN { min=999999; max=0; sum=0; n=0; above100=0; above1000=0; above5000=0 }
    {
        depth = $3
        sum += depth
        n++
        if (depth < min) min = depth
        if (depth > max) max = depth
        if (depth >= 100) above100++
        if (depth >= 1000) above1000++
        if (depth >= 5000) above5000++
    }
    END {
        if (n == 0) {
            printf "%s\t0\t0\t0\t0\t0\t0\n", sample
        } else {
            mean = sum / n
            printf "%s\t%.1f\t%d\t%d\t%.1f\t%.1f\t%.1f\n", \
                sample, mean, min, max, \
                (above100/n)*100, (above1000/n)*100, (above5000/n)*100
        }
    }' "$output_file" >> "$SUMMARY"
done

echo ""
echo "✓ Coverage summary saved to: $SUMMARY"
echo ""
echo "Quick stats:"
column -t "$SUMMARY" | head -20
echo ""
echo "Use this to tune IVAR_MIN_VARIANT_DEPTH and LOFREQ_MIN_VARIANT_DEPTH in config.sh"
