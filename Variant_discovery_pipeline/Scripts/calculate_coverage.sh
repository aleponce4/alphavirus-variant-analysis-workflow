#!/bin/bash

# Calculate per-position depth for all BAM files
# Output: One file per sample in Input/Coverage/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${REPO_DIR}/config.sh"

INPUT_DIR="${BAM_OUTPUT_DIR:-${REPO_DIR}/Input/BAMs}"
OUTPUT_DIR="${REPO_DIR}/Input/Coverage"
TARGET_CONTIG="${TARGET_CONTIG:-target_contig}"

mkdir -p "$OUTPUT_DIR"
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: BAM input directory not found: ${INPUT_DIR}"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools command not found"
    exit 1
fi

echo "Calculating coverage for all samples..."
found_bam=0

for bam_file in "$INPUT_DIR"/*.bam; do
    if [ ! -f "$bam_file" ]; then
        continue
    fi
    found_bam=1
    sample_name=$(basename "$bam_file" .bam)
    output_file="$OUTPUT_DIR/${sample_name}_coverage.txt"

    echo "Processing: $sample_name"
    if [ -n "$TARGET_CONTIG" ]; then
        samtools depth -a -r "$TARGET_CONTIG" "$bam_file" > "$output_file"
    else
        samtools depth -a "$bam_file" > "$output_file"
    fi
done

if [ "$found_bam" -eq 0 ]; then
    echo "ERROR: no BAM files found in ${INPUT_DIR}"
    exit 1
fi

echo "OK Done. Coverage files saved in $OUTPUT_DIR/"
