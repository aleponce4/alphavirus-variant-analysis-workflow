#!/bin/bash

# How to run: ./Scripts/run_ivar.sh
# iVar Pipeline for multi-sample viral variant calling and consensus generation
# Processes all BAM files in Input/BAMs/ directory
source ./config.sh

# Setup logging
LOG_FILE="Ivar/run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p Ivar
exec > >(tee -a "$LOG_FILE") 2>&1

eval "$(conda shell.bash hook)"
conda activate ivar_env

echo "Starting iVar pipeline..."

# Clean previous results
echo "Cleaning previous iVar results..."
rm -rf Ivar/*

# Find reference file
reference=$(find Input/Reference -name "*.fasta" -o -name "*.fa" | head -n 1)
if [ -z "$reference" ]; then
    echo "ERROR: No reference file found in Input/Reference/"
    exit 1
fi

# Check for primer file (optional)
primer_bed=$(find Input/Primers -name "*.bed" 2>/dev/null | head -n 1)

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
    out_dir="Ivar/${sample_name}"
    mkdir -p "$out_dir"
    
    echo "Processing sample: $sample_name"
    
    # Restrict to viral contig first to avoid host reads
    echo "Extracting viral reads from contig: $VIRAL_CONTIG"
    samtools view -b "$bam_file" "$VIRAL_CONTIG" | samtools sort -o "$out_dir/viral_only.bam" -
    samtools index "$out_dir/viral_only.bam"
    viral_bam="$out_dir/viral_only.bam"
    
    # Conditional primer trimming
    if [ -n "$primer_bed" ] && [ -f "$primer_bed" ]; then
        echo "Trimming primers using: $(basename "$primer_bed")"
        ivar trim -i "$viral_bam" -b "$primer_bed" -p "$out_dir/trimmed" \
            -m $IVAR_PRIMER_TRIM_QUALITY -q $IVAR_PRIMER_BASE_QUALITY -s $IVAR_PRIMER_WINDOW_SIZE -e
        
        # Sort and index trimmed BAM
        samtools sort -o "$out_dir/trimmed.sorted.bam" "$out_dir/trimmed.bam"
        samtools index "$out_dir/trimmed.sorted.bam"
        input_bam="$out_dir/trimmed.sorted.bam"
    else
        echo "No primer file found. Skipping primer trimming step."
        # Use viral-only BAM file
        input_bam="$viral_bam"
    fi
    
    # Generate consensus
    echo "Creating consensus..."
    samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$input_bam" | \
        ivar consensus -p "$out_dir/consensus" -m $IVAR_MIN_CONSENSUS_COVERAGE -t 0.5 -q $IVAR_MIN_BASE_QUALITY
    
    # Call variants
    echo "Calling variants..."
    samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$input_bam" | \
        ivar variants -p "$out_dir/variants" -r "$reference" \
        -m $IVAR_MIN_VARIANT_DEPTH -t $IVAR_MIN_VARIANT_FREQ
    
    # Generate basic QC stats
    echo "Generating basic QC stats..."
    echo "Sample: $sample_name" > "$out_dir/qc_stats.txt"
    echo "Reference length: $(grep -v '^>' "$reference" | tr -d '\n' | wc -c)" >> "$out_dir/qc_stats.txt"
    if [ -f "$out_dir/consensus.fa" ]; then
        echo "Consensus length: $(grep -v '^>' "$out_dir/consensus.fa" | tr -d '\n' | wc -c)" >> "$out_dir/qc_stats.txt"
        echo "N positions: $(grep -v '^>' "$out_dir/consensus.fa" | tr -d '\n' | grep -o 'N' | wc -l)" >> "$out_dir/qc_stats.txt"
    fi
    if [ -f "$out_dir/variants.tsv" ]; then
        variant_count=$(tail -n +2 "$out_dir/variants.tsv" | wc -l)
        echo "Variants called: $variant_count" >> "$out_dir/qc_stats.txt"
    fi
    
    # Clean up intermediate files
    if [ -f "$out_dir/trimmed.bam" ]; then
        rm -f "$out_dir/trimmed.bam"
    fi
    
    echo "Completed processing: $sample_name"
    echo ""
done

echo "iVar pipeline completed. Processed all samples."
