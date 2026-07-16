#!/bin/bash
# =====================================================================
# prepare_inputs.sh — Stage BAM files for a dataset
#
# Creates symlinks from nf-core star_salmon output to the variant
# pipeline's expected Input/BAMs/ directory.
#
# Usage:
#   DATASET=mouse_veev bash prepare_inputs.sh
#   DATASET=mouse_eeev bash prepare_inputs.sh
#   DATASET=rat_veev   bash prepare_inputs.sh
# =====================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Source config
export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

echo "═══════════════════════════════════════════════════════════"
echo "  Staging BAM files for: $DATASET"
echo "═══════════════════════════════════════════════════════════"

# Source BAM directory (nf-core star_salmon output, at project root level)
BAM_SOURCE="$PROJECT_ROOT/nfcore_results/${DATASET}/star_salmon"

echo "  Source: $BAM_SOURCE"
echo "  Target: $INPUT_DIR"
echo ""

mkdir -p "$INPUT_DIR"

# Check source exists
if [ ! -d "$BAM_SOURCE" ]; then
    echo "ERROR: BAM source directory not found: $BAM_SOURCE"
    echo "  Expected nf-core star_salmon output at: $DATASET/star_salmon/"
    echo "  BAMs may still be downloading."
    exit 1
fi

# Find all markdup.sorted.bam files
bam_count=0
missing_index=0

for bam_file in "$BAM_SOURCE"/*.markdup.sorted.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "No .markdup.sorted.bam files found in $BAM_SOURCE"
        break
    fi

    sample_name=$(basename "$bam_file" .markdup.sorted.bam)
    link_path="$INPUT_DIR/${sample_name}.bam"

    # Create symlink (relative for portability)
    ln -sf "$bam_file" "$link_path"
    bam_count=$((bam_count + 1))

    # Check for index
    bai_file="${bam_file}.bai"
    if [ ! -f "$bai_file" ]; then
        echo "  ⚠ Missing index: $(basename "$bai_file")"
        missing_index=$((missing_index + 1))
    else
        ln -sf "$bai_file" "$INPUT_DIR/${sample_name}.bam.bai"
    fi
done

echo ""
echo "  Staged BAM files: $bam_count"
echo "  Missing indices:  $missing_index"

if [ "$missing_index" -gt 0 ]; then
    echo ""
    echo "  ⚠ Some BAMs lack indices. Run:"
    echo "    for f in $INPUT_DIR/*.bam; do samtools index \"\$f\"; done"
fi

# Verify viral contig is present in at least one BAM
if [ "$bam_count" -gt 0 ]; then
    echo ""
    echo "  Checking viral contig '$VIRAL_CONTIG' in first BAM..."
    first_bam=$(ls "$INPUT_DIR"/*.bam 2>/dev/null | head -1)
    if [ -n "$first_bam" ]; then
        if samtools view -H "$first_bam" 2>/dev/null | grep -q "SN:${VIRAL_CONTIG}"; then
            echo "  ✓ Viral contig '$VIRAL_CONTIG' found in BAM header"
        else
            echo "  ✗ Viral contig '$VIRAL_CONTIG' NOT found in BAM header!"
            echo "    Available contigs (last 5):"
            samtools view -H "$first_bam" 2>/dev/null | grep "^@SQ" | tail -5
        fi
    fi
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  Input staging complete for: $DATASET"
echo "═══════════════════════════════════════════════════════════"
