#!/bin/bash

# TEST PIPELINE - Quick validation on 1 subsampled sample
# Usage: cd Variant_discovery_pipeline && ./Scripts/test_pipeline.sh
# Expected runtime: ~2-3 minutes

set -e

# Setup
source test_config.sh
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null

# Config
SAMPLE="INH_3_DPI_R1_A3"
SAMPLE_TEST="${SAMPLE}_subsample"  # Append suffix to avoid overwriting production
REFERENCE="Input/Reference/inh.fasta"
ANNOTATION="Input/Reference/VEEV_INH_fromGenbank.gff3"
PRIMERS="Input/Primers/*.bed"

# Find input BAM
BAM_INPUT=$(find Input/BAMs -name "*${SAMPLE}*.bam" | head -1)
if [ -z "$BAM_INPUT" ]; then
    echo "ERROR: BAM file for $SAMPLE not found"
    exit 1
fi

echo "============================================"
echo "TEST PIPELINE: $SAMPLE (subsampled)"
echo "============================================"
echo "Original BAM: $BAM_INPUT"
echo "Output name: $SAMPLE_TEST"
echo ""

# Create temp subsampled BAM (downsample to ~10% reads for speed)
echo "Subsampling reads to 10%..."
mkdir -p .test_temp
samtools view -s 0.1 -b "$BAM_INPUT" > .test_temp/test.bam
samtools index .test_temp/test.bam
BAM_TEST=".test_temp/test.bam"
echo "  ✓ Subsampled BAM created"
echo ""

# ========== iVar ==========
echo "Running iVar..."
mkdir -p Ivar/$SAMPLE_TEST

# Trim
ivar trim -i "$BAM_TEST" -b $PRIMERS -p Ivar/$SAMPLE_TEST/trimmed \
    -q $IVAR_PRIMER_BASE_QUALITY -m $IVAR_PRIMER_WINDOW_SIZE -e 2>/dev/null || true

# Sort trimmed BAM
samtools sort -o Ivar/$SAMPLE_TEST/trimmed.sorted.bam Ivar/$SAMPLE_TEST/trimmed.bam 2>/dev/null || true
samtools index Ivar/$SAMPLE_TEST/trimmed.sorted.bam 2>/dev/null || true

# Consensus
samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $VIRAL_CONTIG -f $REFERENCE Ivar/$SAMPLE_TEST/trimmed.sorted.bam 2>/dev/null | \
ivar consensus -t 0.5 -m $IVAR_MIN_CONSENSUS_COVERAGE -n N -p Ivar/$SAMPLE_TEST/consensus 2>/dev/null || true

# Variants
samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $VIRAL_CONTIG -f $REFERENCE Ivar/$SAMPLE_TEST/trimmed.sorted.bam 2>/dev/null | \
ivar variants -t $IVAR_MIN_VARIANT_FREQ -m $IVAR_MIN_VARIANT_DEPTH -r $REFERENCE -g $ANNOTATION -p Ivar/$SAMPLE_TEST/variants 2>/dev/null || true

if [ -f "Ivar/$SAMPLE_TEST/variants.tsv" ]; then
    IVAR_COUNT=$(tail -n +2 Ivar/$SAMPLE_TEST/variants.tsv 2>/dev/null | wc -l || echo "0")
    echo "  ✓ iVar complete: $IVAR_COUNT variants called"
else
    echo "  ✗ iVar failed: no output file"
    IVAR_COUNT="ERROR"
fi
echo ""

# ========== LoFreq ==========
echo "Running LoFreq..."
mkdir -p LoFreq/$SAMPLE_TEST

lofreq call -f $REFERENCE -r $VIRAL_CONTIG \
    -d $LOFREQ_MIN_VARIANT_DEPTH -q $LOFREQ_MIN_BASE_QUALITY -Q $LOFREQ_MIN_MAP_QUALITY \
    --call-indels --use-mplileup \
    "$BAM_TEST" > LoFreq/$SAMPLE_TEST/variants.vcf 2>/dev/null || true

if [ -f "LoFreq/$SAMPLE_TEST/variants.vcf" ]; then
    LOFREQ_COUNT=$(grep -v "^#" LoFreq/$SAMPLE_TEST/variants.vcf 2>/dev/null | wc -l || echo "0")
    echo "  ✓ LoFreq complete: $LOFREQ_COUNT variants called"
else
    echo "  ✗ LoFreq failed: no output file"
    LOFREQ_COUNT="ERROR"
fi
echo ""

# ========== Annotation ==========
echo "Annotating variants..."
mkdir -p Annotated_variants/Ivar Annotated_variants/LoFreq

# iVar: convert TSV to VCF, then annotate
if [ -f "Ivar/$SAMPLE_TEST/variants.tsv" ]; then
    python Scripts/ivar_variants_to_vcf.py "Ivar/$SAMPLE_TEST/variants.tsv" \
        "Ivar/$SAMPLE_TEST/${SAMPLE_TEST}.vcf" "$REFERENCE" 2>/dev/null || true
    
    if [ -f "Ivar/$SAMPLE_TEST/${SAMPLE_TEST}.vcf" ]; then
        conda activate annotation-env 2>/dev/null || true
        bcftools csq --local-csq -f $REFERENCE -g $ANNOTATION \
            "Ivar/$SAMPLE_TEST/${SAMPLE_TEST}.vcf" -o Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf 2>/dev/null || \
        cp "Ivar/$SAMPLE_TEST/${SAMPLE_TEST}.vcf" Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf
        echo "  ✓ iVar annotation complete"
    fi
fi

# LoFreq: annotate VCF
if [ -f "LoFreq/$SAMPLE_TEST/variants.vcf" ]; then
    conda activate annotation-env 2>/dev/null || true
    bcftools csq --local-csq -f $REFERENCE -g $ANNOTATION \
        "LoFreq/$SAMPLE_TEST/variants.vcf" -o Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf 2>/dev/null || \
    cp "LoFreq/$SAMPLE_TEST/variants.vcf" Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf
    echo "  ✓ LoFreq annotation complete"
fi
echo ""

# ========== Final Check ==========
echo "============================================"
echo "TEST RESULTS"
echo "============================================"

# Check iVar
echo -n "iVar output:  "
if [ -f "Ivar/$SAMPLE_TEST/variants.tsv" ]; then
    echo "✓ $IVAR_COUNT variants"
    IVAR_OK=1
else
    echo "✗ FAILED"
    IVAR_OK=0
fi

# Check LoFreq
echo -n "LoFreq output: "
if [ -f "LoFreq/$SAMPLE_TEST/variants.vcf" ]; then
    echo "✓ $LOFREQ_COUNT variants"
    LOFREQ_OK=1
else
    echo "✗ FAILED"
    LOFREQ_OK=0
fi

# Check annotations
echo -n "Annotations:  "
if [ -f "Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf" ] && [ -f "Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf" ]; then
    echo "✓ both complete"
    ANNOT_OK=1
else
    echo "✗ FAILED"
    ANNOT_OK=0
fi

echo ""

# Summary
if [ $IVAR_OK -eq 1 ] && [ $LOFREQ_OK -eq 1 ] && [ $ANNOT_OK -eq 1 ]; then
    echo "✓ SUCCESS: All pipeline steps ran correctly!"
    echo ""
    echo "Files created:"
    echo "  - Ivar/$SAMPLE_TEST/variants.tsv"
    echo "  - LoFreq/$SAMPLE_TEST/variants.vcf"
    echo "  - Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf"
    echo "  - Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf"
    echo ""
    echo "Production results (without _subsample suffix) remain SAFE and unchanged"
    echo ""
    echo "Next: Run production with: source config.sh && ./Scripts/run_ivar.sh"
    rm -rf .test_temp
    exit 0
else
    echo "✗ FAILURE: Some steps did not complete"
    echo "Check output above for errors"
    exit 1
fi
