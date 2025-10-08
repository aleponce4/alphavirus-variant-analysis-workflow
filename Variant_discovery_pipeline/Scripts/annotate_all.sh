#!/bin/bash

# Annotate all VCF files using bcftools csq
# CRITICAL: Uses --local-csq for iVar files to fix amino acid position errors

# Activate conda environment with bcftools
eval "$(conda shell.bash hook)"
conda activate annotation-env

# Configuration
REFERENCE="Input/Reference/inh.fasta"
VIRAL_REFERENCE="/tmp/viral_only.fasta"  # Viral-only reference for annotation
ANNOTATION="Input/Reference/VEEV_INH_fromGenbank.gff3" # GFF3 annotation file 
OUTPUT_DIR="Annotated_variants"
LOFREQ_DIR="$OUTPUT_DIR/LoFreq"
IVAR_DIR="$OUTPUT_DIR/Ivar"

# Create output directories
mkdir -p "$LOFREQ_DIR"
mkdir -p "$IVAR_DIR"

echo "Starting variant annotation pipeline..."
echo "Reference: $REFERENCE"
echo "Annotation: $ANNOTATION"
echo "Output directory: $OUTPUT_DIR"

# Check required files exist
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference file not found: $REFERENCE"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo "ERROR: Annotation file not found: $ANNOTATION"
    exit 1
fi

echo "Files validated. Starting annotation..."

# Function to annotate LoFreq files (bcftools needs a GFF even with --local-csq)
annotate_lofreq() {
  local input_file=$1
  local output_file=$2

  echo "  Annotating LoFreq: $(basename "$input_file") [local-csq mode]"
  if bcftools csq \
      -f "$VIRAL_REFERENCE" \
      -g "$ANNOTATION" \
      --local-csq \
      "$input_file" \
      -o "$output_file"; then
    echo "    ✓ Success"
  else
    echo "    ✗ Failed (exit code $?)"
    return 1
  fi
}

# Function to annotate iVar files (bcftools needs a GFF even with --local-csq)
annotate_ivar() {
  local input_file=$1
  local output_file=$2

  echo "  Annotating iVar: $(basename "$input_file") [local-csq mode]"
  if bcftools csq \
      -f "$VIRAL_REFERENCE" \
      -g "$ANNOTATION" \
      --local-csq \
      "$input_file" \
      -o "$output_file"; then
    echo "    ✓ Success"
  else
    echo "    ✗ Failed (exit code $?)"
    return 1
  fi
}


# Process LoFreq files
echo ""
echo "Processing LoFreq VCF files..."
lofreq_count=0
for sample_dir in LoFreq/*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        
        # Check for filtered VCF first
        if [ -f "$sample_dir/variants.filtered.vcf.gz" ]; then
            vcf_file="$sample_dir/variants.filtered.vcf.gz"
            output_file="$LOFREQ_DIR/${sample_name}_filtered.vcf"
            
            # Check if VCF has variants (skip if empty to avoid segfault)
            variant_count=$(zcat "$vcf_file" | grep -v "^#" | wc -l)
            if [ $variant_count -gt 0 ]; then
                echo "  Processing $sample_name ($variant_count variants)"
                if annotate_lofreq "$vcf_file" "$output_file"; then
                    ((lofreq_count++))
                else
                    echo "    Failed to annotate $sample_name - continuing..."
                fi
            else
                echo "  Skipping $sample_name (no variants - would cause segfault)"
            fi
            
        elif [ -f "$sample_dir/variants.vcf" ]; then
            vcf_file="$sample_dir/variants.vcf"
            output_file="$LOFREQ_DIR/${sample_name}.vcf"
            
            # Check if VCF has variants (skip if empty to avoid segfault)
            variant_count=$(grep -v "^#" "$vcf_file" | wc -l)
            if [ $variant_count -gt 0 ]; then
                echo "  Processing $sample_name ($variant_count variants)"
                if annotate_lofreq "$vcf_file" "$output_file"; then
                    ((lofreq_count++))
                else
                    echo "    Failed to annotate $sample_name - continuing..."
                fi
            else
                echo "  Skipping $sample_name (no variants - would cause segfault)"
            fi
        fi
    fi
done

if [ $lofreq_count -eq 0 ]; then
    echo "  No LoFreq VCF files found in LoFreq/*/"
else
    echo "  Processed $lofreq_count LoFreq files"
fi

# Process iVar files (convert TSV to VCF first)
echo ""
echo "Processing iVar TSV files..."
ivar_count=0
for sample_dir in Ivar/*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        
        if [ -f "$sample_dir/variants.tsv" ]; then
            # Check if TSV has variants first (skip if empty to avoid segfault)
            variant_count=$(tail -n +2 "$sample_dir/variants.tsv" | wc -l)
            
            if [ $variant_count -gt 0 ]; then
                # Convert TSV to VCF using Python script
                tsv_file="$sample_dir/variants.tsv"
                temp_vcf="$sample_dir/variants.vcf"
                
                echo "  Converting TSV to VCF: $sample_name ($variant_count variants)"
                python Scripts/ivar_variants_to_vcf.py "$tsv_file" "$temp_vcf" "$REFERENCE"
                
                if [ -f "$temp_vcf" ]; then
                    output_file="$IVAR_DIR/${sample_name}.vcf"
                    echo "  Annotating iVar: $sample_name"
                    if annotate_ivar "$temp_vcf" "$output_file"; then
                        ((ivar_count++))
                    else
                        echo "    Failed to annotate $sample_name - continuing..."
                    fi
                    
                    # Clean up temporary VCF
                    rm -f "$temp_vcf"
                fi
            else
                echo "  Skipping $sample_name (no variants - would cause segfault)"
            fi
        fi
    fi
done

if [ $ivar_count -eq 0 ]; then
    echo "  No iVar TSV files found in Ivar/*/"
else
    echo "  Processed $ivar_count iVar files"
fi

echo ""
echo "Annotation completed. Results in $OUTPUT_DIR/"
echo "LoFreq annotated files: $lofreq_count"
echo "iVar annotated files: $ivar_count"
echo "Total files processed: $((lofreq_count + ivar_count))"