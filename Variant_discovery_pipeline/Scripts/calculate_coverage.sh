#!/bin/bash

# Calculate per-position depth for all BAM files
# Output: One file per sample in Input/Coverage/
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null || echo "Warning: No conda environment activated"

INPUT_DIR="Input/BAMs"
OUTPUT_DIR="Input/Coverage"
VIRAL_CONTIG="VEEV_INH"

mkdir -p "$OUTPUT_DIR"

echo "Calculating coverage for all samples..."

for bam_file in "$INPUT_DIR"/*.bam; do
    sample_name=$(basename "$bam_file" .bam)
    output_file="$OUTPUT_DIR/${sample_name}_coverage.txt"
    
    echo "Processing: $sample_name"
    samtools depth -a -r "$VIRAL_CONTIG" "$bam_file" > "$output_file"
done

echo "✓ Done. Coverage files saved in $OUTPUT_DIR/"