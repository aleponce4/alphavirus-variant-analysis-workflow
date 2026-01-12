#!/bin/bash

# Run full pipeline on a single sample
# Usage: ./Scripts/run_single_sample.sh INH_3_DPI_R1_A3

source config.sh
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null

if [ -z "$1" ]; then
    echo "Usage: ./Scripts/run_single_sample.sh <sample_name>"
    echo "Example: ./Scripts/run_single_sample.sh INH_3_DPI_R1_A3"
    exit 1
fi

SAMPLE="$1"
BAM="Input/BAMs/${SAMPLE}.bam"
REFERENCE="Input/Reference/inh.fasta"
ANNOTATION="Input/Reference/VEEV_INH_fromGenbank.gff3"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

echo "=================================================="
echo "Running full pipeline on: $SAMPLE"
echo "=================================================="
echo ""

# ========== iVar ==========
echo "iVar: primer trimming..."
mkdir -p Ivar/$SAMPLE
ivar trim -i "$BAM" -b Input/Primers/*.bed -p Ivar/$SAMPLE/trimmed \
    -q $IVAR_PRIMER_BASE_QUALITY -m $IVAR_PRIMER_WINDOW_SIZE -e 2>&1 | grep -v "samtools:" || true

samtools sort -o Ivar/$SAMPLE/trimmed.sorted.bam Ivar/$SAMPLE/trimmed.bam 2>/dev/null
samtools index Ivar/$SAMPLE/trimmed.sorted.bam 2>/dev/null

echo "iVar: calling consensus..."
samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $VIRAL_CONTIG -f $REFERENCE Ivar/$SAMPLE/trimmed.sorted.bam 2>/dev/null | \
ivar consensus -t 0.5 -m $IVAR_MIN_CONSENSUS_COVERAGE -n N -p Ivar/$SAMPLE/consensus 2>/dev/null

echo "iVar: calling variants..."
samtools mpileup -aa -A -d 600000 -B -Q $IVAR_MIN_BASE_QUALITY -r $VIRAL_CONTIG -f $REFERENCE Ivar/$SAMPLE/trimmed.sorted.bam 2>/dev/null | \
ivar variants -t $IVAR_MIN_VARIANT_FREQ -m $IVAR_MIN_VARIANT_DEPTH -r $REFERENCE -g $ANNOTATION -p Ivar/$SAMPLE/variants 2>/dev/null

if [ -f "Ivar/$SAMPLE/variants.tsv" ]; then
    COUNT=$(tail -n +2 Ivar/$SAMPLE/variants.tsv 2>/dev/null | wc -l)
    echo "  ✓ iVar: $COUNT variants"
fi
echo ""

# ========== LoFreq ==========
echo "LoFreq: calling variants..."
mkdir -p LoFreq/$SAMPLE
lofreq call -f $REFERENCE -r $VIRAL_CONTIG \
    -d $LOFREQ_MIN_VARIANT_DEPTH -q $LOFREQ_MIN_BASE_QUALITY -Q $LOFREQ_MIN_MAP_QUALITY \
    --call-indels --use-mplileup \
    "$BAM" > LoFreq/$SAMPLE/variants.vcf 2>/dev/null

if [ -f "LoFreq/$SAMPLE/variants.vcf" ]; then
    COUNT=$(grep -v "^#" LoFreq/$SAMPLE/variants.vcf 2>/dev/null | wc -l)
    echo "  ✓ LoFreq: $COUNT variants"
fi
echo ""

# ========== Annotation ==========
echo "Annotation: iVar..."
mkdir -p Annotated_variants/Ivar Annotated_variants/LoFreq

python Scripts/ivar_variants_to_vcf.py "Ivar/$SAMPLE/variants.tsv" \
    "Ivar/$SAMPLE/${SAMPLE}.vcf" "$REFERENCE" 2>/dev/null || true

conda activate annotation-env 2>/dev/null || true
if [ -f "Ivar/$SAMPLE/${SAMPLE}.vcf" ]; then
    bcftools csq --local-csq -f $REFERENCE -g $ANNOTATION \
        "Ivar/$SAMPLE/${SAMPLE}.vcf" -o Annotated_variants/Ivar/${SAMPLE}_annotated.vcf 2>/dev/null || \
    cp "Ivar/$SAMPLE/${SAMPLE}.vcf" Annotated_variants/Ivar/${SAMPLE}_annotated.vcf
fi

echo "Annotation: LoFreq..."
if [ -f "LoFreq/$SAMPLE/variants.vcf" ]; then
    bcftools csq --local-csq -f $REFERENCE -g $ANNOTATION \
        "LoFreq/$SAMPLE/variants.vcf" -o Annotated_variants/LoFreq/${SAMPLE}_annotated.vcf 2>/dev/null || \
    cp "LoFreq/$SAMPLE/variants.vcf" Annotated_variants/LoFreq/${SAMPLE}_annotated.vcf
fi

echo ""
echo "=================================================="
echo "✓ Complete: $SAMPLE"
echo "=================================================="
