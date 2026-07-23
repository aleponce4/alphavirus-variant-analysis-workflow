#!/bin/bash
# =====================================================================
# annotate_all.sh — Annotate variants with bcftools csq
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Dataset-aware paths
#   - Requires viral_only.gff3 with CDS features (see PLAN.md §2.1)
#
# Usage:
#   DATASET=mouse_veev bash annotate_all.sh
# =====================================================================
set -uo pipefail
# Note: no 'set -e' — annotation failures on individual samples should not kill the whole run

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Activate conda
eval "$(conda shell.bash hook)"
conda activate annotation-env

mkdir -p "$ANNOTATED_OUT/LoFreq" "$ANNOTATED_OUT/Ivar"

echo "═══════════════════════════════════════════════════════════"
echo "  Annotating variants for: $DATASET"
echo "  Reference:  $REFERENCE"
echo "  Annotation: $ANNOTATION"
echo "═══════════════════════════════════════════════════════════"

# Check GFF3 exists
if [ ! -f "$ANNOTATION" ]; then
    echo "ERROR: Annotation GFF3 not found: $ANNOTATION"
    echo "  Generate CDS annotations first (see PLAN.md §2.1)"
    exit 1
fi

# Check for CDS features in GFF3
if ! grep -q "CDS" "$ANNOTATION"; then
    echo "ERROR: No CDS features found in $ANNOTATION"
    echo "  bcftools csq requires CDS features for amino acid prediction"
    exit 1
fi

# --- Annotate LoFreq VCFs ---
echo ""
echo "Processing LoFreq VCF files..."
lofreq_count=0

for sample_dir in "$LOFREQ_OUT"/*/; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi

    sample_name=$(basename "$sample_dir")

    if [ -f "$sample_dir/variants.filtered.vcf.gz" ]; then
        vcf_file="$sample_dir/variants.filtered.vcf.gz"
        output_file="$ANNOTATED_OUT/LoFreq/${sample_name}_filtered.vcf"

        FORCE_RECALL="${FORCE_RECALL:-false}"
        if [ "$FORCE_RECALL" = "false" ] && [ -f "$output_file" ]; then
            echo "  ✓ LoFreq: $sample_name already annotated, skipping"
            lofreq_count=$((lofreq_count + 1))
            continue
        fi

        echo "  Annotating LoFreq: $sample_name"
        if bcftools csq -f "$REFERENCE" -g "$ANNOTATION" --local-csq \
            "$vcf_file" -o "$output_file"; then
            echo "    ✓ Success"
            lofreq_count=$((lofreq_count + 1))
        else
            echo "    ✗ Failed (exit code $?)"
        fi
    fi
done

# --- Annotate iVar VCFs ---
echo ""
echo "Processing iVar variant files..."
ivar_count=0

for sample_dir in "$IVAR_OUT"/*/; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi

    sample_name=$(basename "$sample_dir")

    if [ -f "$sample_dir/variants.tsv" ]; then
        variant_count=$(tail -n +2 "$sample_dir/variants.tsv" | wc -l)

        if [ "$variant_count" -gt 0 ]; then
            tsv_file="$sample_dir/variants.tsv"
            temp_vcf="$sample_dir/variants.vcf"
            output_file="$ANNOTATED_OUT/Ivar/${sample_name}.vcf"

            FORCE_RECALL="${FORCE_RECALL:-false}"
            if [ "$FORCE_RECALL" = "false" ] && [ -f "$output_file" ]; then
                echo "  ✓ iVar: $sample_name already annotated, skipping"
                ivar_count=$((ivar_count + 1))
                continue
            fi

            echo "  Converting TSV→VCF: $sample_name ($variant_count variants)"
            python "$SCRIPT_DIR/ivar_variants_to_vcf.py" \
                "$tsv_file" "$temp_vcf" "$REFERENCE"

            if [ -f "$temp_vcf" ]; then
                echo "  Annotating iVar: $sample_name"
                if bcftools csq -f "$REFERENCE" -g "$ANNOTATION" --local-csq \
                    "$temp_vcf" -o "$output_file"; then
                    ivar_count=$((ivar_count + 1))
                else
                    echo "    Failed to annotate $sample_name"
                fi
                rm -f "$temp_vcf"
            fi
        else
            echo "  Skipping $sample_name (no variants)"
        fi
    fi
done

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  Annotation complete for: $DATASET"
echo "  LoFreq annotated: $lofreq_count"
echo "  iVar annotated:   $ivar_count"
echo "  Output: $ANNOTATED_OUT"
echo "═══════════════════════════════════════════════════════════"
