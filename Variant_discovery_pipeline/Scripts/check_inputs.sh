#!/bin/bash

# How to run: ./Scripts/check_inputs.sh
# Input validation and filename sanitization script
# This script checks and sanitizes filenames with spaces or special characters

# Activate conda environment for samtools
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null || echo "Warning: No conda environment activated"

echo "Checking and sanitizing input files..."

# Function to sanitize filename
sanitize_filename() {
    local filepath="$1"
    local dir=$(dirname "$filepath")
    local filename=$(basename "$filepath")
    local extension="${filename##*.}"
    local basename="${filename%.*}"
    
    # Replace spaces and special characters with underscores
    local clean_basename=$(echo "$basename" | sed 's/[^A-Za-z0-9._-]/_/g')
    local clean_filename="${clean_basename}.${extension}"
    local new_filepath="${dir}/${clean_filename}"
    
    if [ "$filepath" != "$new_filepath" ]; then
        echo "  Renaming: $(basename "$filepath") → $(basename "$new_filepath")"
        mv "$filepath" "$new_filepath"
        return 0
    else
        return 1
    fi
}

# Check Input directory structure
echo ""
echo "Checking directory structure..."
for required_dir in "Input/BAMs" "Input/Reference"; do
    if [ ! -d "$required_dir" ]; then
        echo "WARNING: Required directory not found: $required_dir"
    else
        echo "✓ Found: $required_dir"
    fi
done

# Optional directories
for optional_dir in "Input/Primers" "Input"; do
    if [ -d "$optional_dir" ]; then
        echo "✓ Found: $optional_dir"
    fi
done

# Sanitize BAM files
echo ""
echo "Checking BAM files..."
bam_count=0
if [ -d "Input/BAMs" ]; then
    for bam_file in Input/BAMs/*.bam; do
        if [ -f "$bam_file" ]; then
            sanitize_filename "$bam_file"
            ((bam_count++))
        fi
    done
    
    # Also sanitize BAI files
    for bai_file in Input/BAMs/*.bai; do
        if [ -f "$bai_file" ]; then
            sanitize_filename "$bai_file"
        fi
    done
    
    echo "Found $bam_count BAM files"
else
    echo "ERROR: Input/BAMs directory not found"
    exit 1
fi

# Sanitize reference files
echo ""
echo "Checking reference files..."
ref_count=0
if [ -d "Input/Reference" ]; then
    for ref_file in Input/Reference/*.fasta Input/Reference/*.fa Input/Reference/*.fai; do
        if [ -f "$ref_file" ]; then
            sanitize_filename "$ref_file"
            ((ref_count++))
        fi
    done
    echo "Found $ref_count reference-related files"
else
    echo "WARNING: Input/Reference directory not found"
fi

# Check and sanitize primer files (optional)
echo ""
echo "Checking primer files..."
primer_count=0
if [ -d "Input/Primers" ]; then
    for primer_file in Input/Primers/*.bed; do
        if [ -f "$primer_file" ]; then
            sanitize_filename "$primer_file"
            ((primer_count++))
        fi
    done
    if [ $primer_count -gt 0 ]; then
        echo "Found $primer_count primer files"
    else
        echo "No primer files found (optional)"
    fi
else
    echo "No Primers directory found (optional)"
fi

# Check and sanitize annotation files (optional)
echo ""
echo "Checking annotation files..."
annot_count=0
for annot_file in Input/*.gff3 Input/*.gff Input/*.gtf Input/Reference/*.gff3 Input/Reference/*.gff Input/Reference/*.gtf; do
    if [ -f "$annot_file" ]; then
        sanitize_filename "$annot_file"
        ((annot_count++))
    fi
done
if [ $annot_count -gt 0 ]; then
    echo "Found $annot_count annotation files"
else
    echo "No annotation files found (optional)"
fi

# Index BAM files
echo ""
echo "Indexing BAM files..."
for bam_file in Input/BAMs/*.bam; do
    if [ -f "$bam_file" ]; then
        bai_file="${bam_file}.bai"
        if [ ! -f "$bai_file" ]; then
            echo "  Indexing: $(basename "$bam_file")"
            samtools index "$bam_file"
        else
            echo "  Already indexed: $(basename "$bam_file")"
        fi
    fi
done

# Validate BAM headers
echo ""
echo "Validating BAM headers..."
for bam_file in Input/BAMs/*.bam; do
    if [ -f "$bam_file" ]; then
        sample_name=$(basename "$bam_file" .bam)
        header_check=$(samtools view -H "$bam_file" | head -n 1)
        if [[ $header_check == @HD* ]]; then
            echo "  ✓ Valid header: $sample_name"
        else
            echo "  ✗ Invalid header: $sample_name"
        fi
    fi
done

echo ""
echo "Input validation completed!"
echo "Summary:"
echo "  BAM files: $bam_count"
echo "  Reference files: $ref_count"
echo "  Primer files: $primer_count"
echo "  Annotation files: $annot_count"