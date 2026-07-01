#!/bin/bash

# Run single-sample legacy pipeline path.
# Usage: ./Scripts/run_single_sample.sh <sample_name>

set -euo pipefail

source config.sh
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null

if [ -z "${1:-}" ]; then
    echo "Usage: ./Scripts/run_single_sample.sh <sample_name>"
    echo "Example: ./Scripts/run_single_sample.sh INH_3_DPI_R1_A3"
    exit 1
fi

SAMPLE="$1"
BAM="Input/BAMs/${SAMPLE}.bam"
REFERENCE="${REFERENCE_FASTA}"
ANNOTATION="${ANNOTATION_GFF}"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

echo "=================================================="
echo "Running single-sample pipeline on: $SAMPLE"
echo "=================================================="
echo ""

mkdir -p Ivar/$SAMPLE

# iVar
echo "iVar: primer trimming..."
if ! ivar trim -i "$BAM" -b Input/Primers/*.bed -p Ivar/$SAMPLE/trimmed \
    -q $IVAR_PRIMER_BASE_QUALITY -m $IVAR_PRIMER_WINDOW_SIZE -e; then
    echo "ERROR: iVar trim failed"
    exit 1
fi

if ! samtools sort -o Ivar/$SAMPLE/trimmed.sorted.bam Ivar/$SAMPLE/trimmed.bam; then
    echo "ERROR: sorting trimmed BAM failed"
    exit 1
fi
if ! samtools index Ivar/$SAMPLE/trimmed.sorted.bam; then
    echo "ERROR: indexing trimmed BAM failed"
    exit 1
fi

echo "iVar: calling consensus..."
if ! samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $TARGET_CONTIG -f $REFERENCE Ivar/$SAMPLE/trimmed.sorted.bam \
    | ivar consensus -t 0.5 -m $IVAR_MIN_CONSENSUS_COVERAGE -n N -p Ivar/$SAMPLE/consensus; then
    echo "ERROR: iVar consensus failed"
    exit 1
fi

echo "iVar: calling variants..."
if ! samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $TARGET_CONTIG -f $REFERENCE Ivar/$SAMPLE/trimmed.sorted.bam \
    | ivar variants -t $IVAR_MIN_VARIANT_FREQ -m $IVAR_MIN_VARIANT_DEPTH -r $REFERENCE -g $ANNOTATION -p Ivar/$SAMPLE/variants; then
    echo "ERROR: iVar variants failed"
    exit 1
fi

if [ -f "Ivar/$SAMPLE/variants.tsv" ]; then
    COUNT=$(tail -n +2 Ivar/$SAMPLE/variants.tsv 2>/dev/null | wc -l)
    echo "  iVar: $COUNT variants"
fi
echo ""

# LoFreq
echo "LoFreq: calling variants..."
mkdir -p LoFreq/$SAMPLE
if ! lofreq call -f "$REFERENCE" -r "$TARGET_CONTIG" \
    -d $LOFREQ_MIN_VARIANT_DEPTH -q $LOFREQ_MIN_BASE_QUALITY -Q $LOFREQ_MIN_MAP_QUALITY \
    --call-indels --use-mpileup \
    "$BAM" > LoFreq/$SAMPLE/variants.vcf; then
    echo "ERROR: LoFreq call failed"
    exit 1
fi

if [ -f "LoFreq/$SAMPLE/variants.vcf" ]; then
    COUNT=$(grep -v "^#" LoFreq/$SAMPLE/variants.vcf 2>/dev/null | wc -l)
    echo "  LoFreq: $COUNT variants"
fi
echo ""

# Annotation
echo "Annotation: iVar..."
mkdir -p Annotated_variants/Ivar Annotated_variants/LoFreq
if [ -f "Ivar/$SAMPLE/variants.tsv" ]; then
    python Scripts/ivar_variants_to_vcf.py "Ivar/$SAMPLE/variants.tsv" \
        "Ivar/$SAMPLE/${SAMPLE}.vcf" "$REFERENCE"
fi

if [ -f "Ivar/$SAMPLE/${SAMPLE}.vcf" ]; then
    if conda activate annotation-env 2>/dev/null; then
        if ! bcftools csq --local-csq -f "$REFERENCE" -g "$ANNOTATION" \
            "Ivar/$SAMPLE/${SAMPLE}.vcf" -o Annotated_variants/Ivar/${SAMPLE}_annotated.vcf; then
            cp "Ivar/$SAMPLE/${SAMPLE}.vcf" Annotated_variants/Ivar/${SAMPLE}_annotated.vcf
        fi
    else
        echo "WARNING: annotation-env unavailable; copying iVar VCF."
        cp "Ivar/$SAMPLE/${SAMPLE}.vcf" Annotated_variants/Ivar/${SAMPLE}_annotated.vcf
    fi
else
    echo "WARNING: iVar VCF conversion did not run"
fi

echo "Annotation: LoFreq..."
if [ -f "LoFreq/$SAMPLE/variants.vcf" ]; then
    if conda activate annotation-env 2>/dev/null; then
        if ! bcftools csq --local-csq -f "$REFERENCE" -g "$ANNOTATION" \
            "LoFreq/$SAMPLE/variants.vcf" -o Annotated_variants/LoFreq/${SAMPLE}_annotated.vcf; then
            cp "LoFreq/$SAMPLE/variants.vcf" Annotated_variants/LoFreq/${SAMPLE}_annotated.vcf
        fi
    else
        echo "WARNING: annotation-env unavailable; copying LoFreq VCF."
        cp "LoFreq/$SAMPLE/variants.vcf" Annotated_variants/LoFreq/${SAMPLE}_annotated.vcf
    fi
else
    echo "WARNING: LoFreq output missing"
fi

echo ""
echo "=================================================="
echo "Single-sample pipeline complete: $SAMPLE"
echo "=================================================="

