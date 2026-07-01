#!/bin/bash

# TEST PIPELINE - Quick validation on 1 subsampled sample
# Usage: cd Variant_discovery_pipeline && ./Scripts/test_pipeline.sh

set -e

source test_config.sh
source config.sh
eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null || conda activate lofreq-env 2>/dev/null

SAMPLE="INH_3_DPI_R1_A3"
SAMPLE_TEST="${SAMPLE}_subsample"
REFERENCE="${REFERENCE_FASTA}"
ANNOTATION="${ANNOTATION_GFF}"
PRIMERS="Input/Primers/*.bed"

BAM_INPUT=$(find Input/BAMs -name "*${SAMPLE}*.bam" | head -1)
if [ -z "${BAM_INPUT}" ]; then
    echo "ERROR: BAM file for ${SAMPLE} not found."
    exit 1
fi

echo "============================================"
echo "TEST PIPELINE: ${SAMPLE} (subsampled)"
echo "============================================"
echo "Original BAM: ${BAM_INPUT}"
echo "Output name: ${SAMPLE_TEST}"
echo ""

echo "Subsampling reads to 10%..."
mkdir -p .test_temp
samtools view -s 0.1 -b "${BAM_INPUT}" > .test_temp/test.bam
samtools index .test_temp/test.bam
BAM_TEST=".test_temp/test.bam"
echo "  OK Subsampled BAM created"
echo ""

mkdir -p Ivar/"${SAMPLE_TEST}"

ivar trim -i "${BAM_TEST}" -b ${PRIMERS} -p Ivar/"${SAMPLE_TEST}"/trimmed \
    -q "${IVAR_PRIMER_BASE_QUALITY}" -m "${IVAR_PRIMER_WINDOW_SIZE}" -e
samtools sort -o Ivar/"${SAMPLE_TEST}"/trimmed.sorted.bam Ivar/"${SAMPLE_TEST}"/trimmed.bam
samtools index Ivar/"${SAMPLE_TEST}"/trimmed.sorted.bam
samtools mpileup -aa -A -d 600000 -B -Q "${IVAR_MIN_BASE_QUALITY}" -r "${TARGET_CONTIG}" -f "${REFERENCE}" \
    Ivar/"${SAMPLE_TEST}"/trimmed.sorted.bam | \
    ivar consensus -t 0.5 -m "${IVAR_MIN_CONSENSUS_COVERAGE}" -n N -p Ivar/"${SAMPLE_TEST}"/consensus
samtools mpileup -aa -A -d 600000 -B -Q "${IVAR_MIN_BASE_QUALITY}" -r "${TARGET_CONTIG}" -f "${REFERENCE}" \
    Ivar/"${SAMPLE_TEST}"/trimmed.sorted.bam | \
    ivar variants -t "${IVAR_MIN_VARIANT_FREQ}" -m "${IVAR_MIN_VARIANT_DEPTH}" -r "${REFERENCE}" -g "${ANNOTATION}" \
    -p Ivar/"${SAMPLE_TEST}"/variants

if [ -f "Ivar/${SAMPLE_TEST}/variants.tsv" ]; then
    IVAR_COUNT=$(tail -n +2 "Ivar/${SAMPLE_TEST}/variants.tsv" 2>/dev/null | wc -l || true)
    echo "iVar output: ${IVAR_COUNT} variants"
else
    echo "iVar output missing"
    IVAR_COUNT="ERROR"
fi
echo ""

mkdir -p LoFreq/"${SAMPLE_TEST}"
lofreq call -f "${REFERENCE}" -r "${TARGET_CONTIG}" \
    -d "${LOFREQ_MIN_VARIANT_DEPTH}" -q "${LOFREQ_MIN_BASE_QUALITY}" -Q "${LOFREQ_MIN_MAP_QUALITY}" \
    --call-indels --use-mpileup \
    "${BAM_TEST}" > LoFreq/"${SAMPLE_TEST}"/variants.vcf

if [ -f "LoFreq/${SAMPLE_TEST}/variants.vcf" ]; then
    LOFREQ_COUNT=$(grep -v "^#" "LoFreq/${SAMPLE_TEST}/variants.vcf" 2>/dev/null | wc -l || true)
    echo "LoFreq output: ${LOFREQ_COUNT} variants"
else
    echo "LoFreq output missing"
    LOFREQ_COUNT="ERROR"
fi
echo ""

mkdir -p Annotated_variants/Ivar Annotated_variants/LoFreq
if [ -f "Ivar/${SAMPLE_TEST}/variants.tsv" ]; then
    python Scripts/ivar_variants_to_vcf.py "Ivar/${SAMPLE_TEST}/variants.tsv" \
        "Ivar/${SAMPLE_TEST}/${SAMPLE_TEST}.vcf" "${REFERENCE}"
fi

ivar_annotation_ok=1
if [ -f "Ivar/${SAMPLE_TEST}/${SAMPLE_TEST}.vcf" ]; then
    if ! conda activate annotation-env 2>/dev/null; then
        ivar_annotation_ok=0
    elif ! bcftools csq --local-csq -f "${REFERENCE}" -g "${ANNOTATION}" \
        "Ivar/${SAMPLE_TEST}/${SAMPLE_TEST}.vcf" -o Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf; then
        ivar_annotation_ok=0
    fi
else
    ivar_annotation_ok=0
fi

lffreq_annotation_ok=1
if [ -f "LoFreq/${SAMPLE_TEST}/variants.vcf" ]; then
    if ! conda activate annotation-env 2>/dev/null; then
        lffreq_annotation_ok=0
    elif ! bcftools csq --local-csq -f "${REFERENCE}" -g "${ANNOTATION}" \
        "LoFreq/${SAMPLE_TEST}/variants.vcf" -o Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf; then
        lffreq_annotation_ok=0
    fi
else
    lffreq_annotation_ok=0
fi

echo "============================================"
echo "TEST RESULTS"
echo "============================================"

if [ -f "Ivar/${SAMPLE_TEST}/variants.tsv" ]; then
    echo "iVar: OK (${IVAR_COUNT} variants)"
    IVAR_OK=1
else
    echo "iVar: FAILED"
    IVAR_OK=0
fi

if [ -f "LoFreq/${SAMPLE_TEST}/variants.vcf" ]; then
    echo "LoFreq: OK (${LOFREQ_COUNT} variants)"
    LOFREQ_OK=1
else
    echo "LoFreq: FAILED"
    LOFREQ_OK=0
fi

if [ "${ivar_annotation_ok}" -eq 1 ] && [ "${lffreq_annotation_ok}" -eq 1 ]; then
    echo "Annotations: OK"
    ANNOT_OK=1
else
    echo "Annotations: FAILED"
    ANNOT_OK=0
fi

echo ""

if [ "${IVAR_OK}" -eq 1 ] && [ "${LOFREQ_OK}" -eq 1 ] && [ "${ANNOT_OK}" -eq 1 ]; then
    echo "SUCCESS: pipeline test completed"
    echo ""
    echo "Files created:"
    echo "  - Ivar/${SAMPLE_TEST}/variants.tsv"
    echo "  - LoFreq/${SAMPLE_TEST}/variants.vcf"
    echo "  - Annotated_variants/Ivar/${SAMPLE_TEST}_annotated.vcf"
    echo "  - Annotated_variants/LoFreq/${SAMPLE_TEST}_annotated.vcf"
    echo ""
    rm -rf .test_temp
    exit 0
fi

echo "FAILURE: some steps did not complete"
exit 1

