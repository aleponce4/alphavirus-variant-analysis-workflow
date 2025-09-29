#!/bin/bash

# How to run: ./Scripts/run_lofreq.sh  
# LoFreq Pipeline for multi-sample statistical variant calling
# Processes all BAM files in Input/BAMs/ directory
source ./config.sh

# Setup logging
LOG_FILE="LoFreq/run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p LoFreq
exec > >(tee -a "$LOG_FILE") 2>&1

eval "$(conda shell.bash hook)"
conda activate lofreq-env

echo "Starting LoFreq pipeline..."

# Clean previous results
echo "Cleaning previous LoFreq results..."
rm -rf LoFreq/*

# Find reference file
reference=$(find Input/Reference -name "*.fasta" -o -name "*.fa" | head -n 1)
if [ -z "$reference" ]; then
    echo "ERROR: No reference file found in Input/Reference/"
    exit 1
fi

# Check for BED file (optional)
bed_file=$(find Input -name "*.bed" 2>/dev/null | head -n 1)

# Index reference
echo "Indexing reference file..."
samtools faidx "$reference"

# Process each BAM file
for bam_file in Input/BAMs/*.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in Input/BAMs/"
        continue
    fi
    
    # Extract sample name  
    sample_name=$(basename "$bam_file" .bam)
    out_dir="LoFreq/${sample_name}"
    mkdir -p "$out_dir"
    
    echo "Processing sample: $sample_name"
    
    # Restrict to viral contig before viterbi to avoid spliced host reads
    echo "Extracting viral reads from contig: $VIRAL_CONTIG"
    samtools view -b "$bam_file" "$VIRAL_CONTIG" | samtools sort -o "$out_dir/viral_only.bam" -
    samtools index "$out_dir/viral_only.bam"
    input_bam="$out_dir/viral_only.bam"
    
    # Viterbi realignment (conditional)
    if [ "$LOFREQ_BAQ" -eq 1 ]; then
        echo "Viterbi realignment..."
        lofreq viterbi -f "$reference" "$input_bam" | \
            samtools sort -o "$out_dir/realigned.bam" -
        samtools index "$out_dir/realigned.bam"
        processing_bam="$out_dir/realigned.bam"
    else
        echo "Skipping Viterbi realignment..."
        processing_bam="$input_bam"
    fi
    
    # Add indel qualities (conditional)
    if [ "$LOFREQ_ENABLE_INDELQUAL" -eq 1 ]; then
        echo "Adding indel qualities..."
        lofreq indelqual --dindel -f "$reference" -o "$out_dir/indelqual.bam" "$processing_bam"
        samtools index "$out_dir/indelqual.bam"
        final_bam="$out_dir/indelqual.bam"
    else
        echo "Skipping indel quality adjustment..."
        final_bam="$processing_bam"
    fi
    
    # Call variants
    echo "Calling variants..."
    lofreq call-parallel --pp-threads $THREADS -f "$reference" \
        --min-cov "$LOFREQ_MIN_VARIANT_DEPTH" \
        --min-bq "$LOFREQ_MIN_BASE_QUALITY" \
        --min-alt-bq "$LOFREQ_MIN_BASE_QUALITY" \
        --min-mq "$LOFREQ_MIN_MAP_QUALITY" \
        --sig "$LOFREQ_CALL_SIG" \
        -o "$out_dir/variants.vcf" \
        "$final_bam"
    
    # Filter variants
    echo "Filtering variants..."
    lofreq filter -i "$out_dir/variants.vcf" -o "$out_dir/variants.filtered.vcf" \
        --snvqual-thresh 20 --indelqual-thresh 20
    
    # Compress and index
    bgzip -f "$out_dir/variants.filtered.vcf"
    tabix -f -p vcf "$out_dir/variants.filtered.vcf.gz"
    
    # Generate basic QC stats
    echo "Generating basic QC stats..."
    echo "Sample: $sample_name" > "$out_dir/qc_stats.txt"
    if [ -f "$out_dir/variants.vcf" ]; then
        variant_count=$(grep -c -v '^#' "$out_dir/variants.vcf")
        echo "Raw variants: $variant_count" >> "$out_dir/qc_stats.txt"
    fi
    if [ -f "$out_dir/variants.filtered.vcf.gz" ]; then
        filtered_count=$(zcat "$out_dir/variants.filtered.vcf.gz" | grep -c -v '^#')
        echo "Filtered variants: $filtered_count" >> "$out_dir/qc_stats.txt"
    fi
    
    echo "Completed processing: $sample_name"
    echo ""
done

echo "LoFreq pipeline completed. Processed all samples."
