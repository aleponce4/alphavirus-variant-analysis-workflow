#!/bin/bash
# =====================================================================
# prepare_references.sh — Prepare viral-only reference files
#
# Extracts viral-only FASTA from the full virus reference and creates
# the .fai index. Also checks for CDS annotations in GFF3 format.
#
# Usage:
#   bash variant_analysis/workflow/Scripts/prepare_references.sh
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
SRC_REFS="$PROJECT_ROOT/references"
DST_REFS="$PROJECT_ROOT/variant_analysis/references"

echo "═══════════════════════════════════════════════════════════"
echo "  Preparing viral-only reference files"
echo "═══════════════════════════════════════════════════════════"

# --- VEEV INH-9813 ---
VEEV_SRC="$SRC_REFS/VEEV_INH/virus.fa"
VEEV_DST="$DST_REFS/VEEV_INH"
VEEV_GTF_SRC="$SRC_REFS/VEEV_INH/virus.gtf"

echo ""
echo "── VEEV INH-9813 ──"
mkdir -p "$VEEV_DST"

# Copy viral-only FASTA (it's already viral-only)
cp -f "$VEEV_SRC" "$VEEV_DST/viral_only.fasta"
echo "  ✓ FASTA: viral_only.fasta"

# Index reference
if command -v samtools &>/dev/null; then
    samtools faidx "$VEEV_DST/viral_only.fasta"
    echo "  ✓ Index: viral_only.fasta.fai"
else
    echo "  ⚠ samtools not found — create .fai manually: samtools faidx viral_only.fasta"
fi

# Copy GTF (for SNPGenie)
if [ -f "$VEEV_DST/viral_only.gff3" ]; then
    echo "  Generating SNPGenie-compatible GTF from GFF3..."
    python3 "$SCRIPT_DIR/Helpers/convert_gff3_to_gtf.py" --gff3 "$VEEV_DST/viral_only.gff3" --gtf "$VEEV_DST/viral_only.gtf"
else
    cp -f "$VEEV_GTF_SRC" "$VEEV_DST/viral_only.gtf"
    echo "  ⚠ GFF3 missing, copied default GTF (might lack CDS for SNPGenie)"
fi

# Check for GFF3 with CDS
if [ -f "$VEEV_DST/viral_only.gff3" ]; then
    echo "  ✓ GFF3:  viral_only.gff3 (exists)"
else
    echo "  ✗ GFF3:  MISSING — viral_only.gff3 with CDS features required for bcftools csq"
    echo "          See PLAN.md §2.1 for CDS annotation generation options"
fi

# --- EEEV FL93 ---
EEEV_SRC="$SRC_REFS/EEEV/virus.fa"
EEEV_DST="$DST_REFS/EEEV_FL93"
EEEV_GTF_SRC="$SRC_REFS/EEEV/virus.gtf"

echo ""
echo "── EEEV FL93 ──"
mkdir -p "$EEEV_DST"

cp -f "$EEEV_SRC" "$EEEV_DST/viral_only.fasta"
echo "  ✓ FASTA: viral_only.fasta"

if command -v samtools &>/dev/null; then
    samtools faidx "$EEEV_DST/viral_only.fasta"
    echo "  ✓ Index: viral_only.fasta.fai"
fi

# Copy GTF (for SNPGenie)
if [ -f "$EEEV_DST/viral_only.gff3" ]; then
    echo "  Generating SNPGenie-compatible GTF from GFF3..."
    python3 "$SCRIPT_DIR/Helpers/convert_gff3_to_gtf.py" --gff3 "$EEEV_DST/viral_only.gff3" --gtf "$EEEV_DST/viral_only.gtf"
else
    cp -f "$EEEV_GTF_SRC" "$EEEV_DST/viral_only.gtf"
    echo "  ⚠ GFF3 missing, copied default GTF (might lack CDS for SNPGenie)"
fi

if [ -f "$EEEV_DST/viral_only.gff3" ]; then
    echo "  ✓ GFF3:  viral_only.gff3 (exists)"
else
    echo "  ✗ GFF3:  MISSING — viral_only.gff3 with CDS features required for bcftools csq"
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  Reference preparation complete"
echo "═══════════════════════════════════════════════════════════"
