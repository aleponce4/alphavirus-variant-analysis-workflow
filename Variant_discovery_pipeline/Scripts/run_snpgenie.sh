#!/bin/bash
# =====================================================================
# run_snpgenie.sh — Run SNPGenie selection analysis (dataset-aware)
#
# Adapted from: aleponce4/alphavirus-variant-analysis-workflow
# Changes:
#   - Uses dataset-aware config.sh paths and variables
#   - Auto-clones SNPGenie perl repository if not present
#   - Saves all intermediate and final outputs inside results/$DATASET/SNPGenie
#   - Fixed exit-on-error bug (set -uo pipefail)
#   - Uses manifest-aware helper scripts to dynamically support DPI timepoints
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATASET="${DATASET:-mouse_veev}"
source "$SCRIPT_DIR/../config.sh"

# Paths
SNPGENIE_DIR="${RESULTS_DIR}/SNPGenie"
VCF_DIR="${SNPGENIE_DIR}/vcf_by_sample"
MANIFEST_DIR="${SNPGENIE_DIR}/manifest"
MANIFEST_FILE="${MANIFEST_DIR}/samples.tsv"
RUNS_DIR="${SNPGENIE_DIR}/runs"
OUTPUT_DIR="${SNPGENIE_DIR}/output"
LOGS_DIR="${SNPGENIE_DIR}/logs"
SUMMARY_DIR="${SNPGENIE_DIR}/summary"
ANALYSIS_DIR="${SNPGENIE_DIR}/analysis/delta_selection"
REPORT_DIR="${SNPGENIE_DIR}/analysis/summary"

# SNPGenie perl repo location
REPO_DIR="${VARIANT_ROOT}/workflow/original_repo/SNPGenie/repo"
SNPGENIE_SCRIPT="${REPO_DIR}/snpgenie.pl"

echo "═══════════════════════════════════════════════════════════"
# Check/Clone SNPGenie
if [ ! -f "$SNPGENIE_SCRIPT" ]; then
    echo "  SNPGenie script not found. Cloning repository to:"
    echo "  $REPO_DIR"
    mkdir -p "$(dirname "$REPO_DIR")"
    git clone https://github.com/chasewnelson/SNPGenie.git "$REPO_DIR"
fi
echo "  Using SNPGenie script: $SNPGENIE_SCRIPT"
echo "═══════════════════════════════════════════════════════════"

# Create directories
mkdir -p "${VCF_DIR}" "${MANIFEST_DIR}" "${RUNS_DIR}" "${OUTPUT_DIR}" "${LOGS_DIR}" "${SUMMARY_DIR}" "${ANALYSIS_DIR}" "${REPORT_DIR}"

# 1. Stage VCF inputs from LoFreq outputs
echo "Staging VCF files for $DATASET..."
printf "sample\tsource_filtered_vcf\tstaged_vcf\n" > "${MANIFEST_FILE}"

sample_count=0
for sample_dir in "${LOFREQ_OUT}"/*/; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    sample=$(basename "$sample_dir")
    filtered_vcf_gz="${sample_dir}/variants.filtered.vcf.gz"
    staged_vcf="${VCF_DIR}/${sample}.vcf"

    if [ -f "$filtered_vcf_gz" ]; then
        gzip -dc "$filtered_vcf_gz" > "$staged_vcf"
        
        # Ensure it's a valid VCF
        if grep -q '^#CHROM' "$staged_vcf"; then
            printf "%s\t%s\t%s\n" "${sample}" "${filtered_vcf_gz}" "${staged_vcf}" >> "${MANIFEST_FILE}"
            sample_count=$((sample_count + 1))
        else
            echo "  ⚠ Warning: Invalid VCF header in $staged_vcf, skipping"
            rm -f "$staged_vcf"
        fi
    fi
done

echo "Staged $sample_count samples."

if [ "$sample_count" -eq 0 ]; then
    echo "No valid VCF files found to analyze. Exiting."
    exit 0
fi

# References
FASTA_FILE="${REFERENCE}"
GTF_FILE="${REFERENCE%.fasta}.gtf"

if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: SNPGenie-compatible GTF file not found: $GTF_FILE"
    echo "Please run Scripts/prepare_references.sh first."
    exit 1
fi

DATE_TAG="$(date +%Y%m%d_%H%M%S)"

# 2. Run SNPGenie Perl Script
run_one_threshold() {
    local minfreq="$1"
    local label="$2"
    local threshold_root="${RUNS_DIR}/${label}"
    local out_root="${OUTPUT_DIR}/${label}"
    local status_tsv="${LOGS_DIR}/run_status_${label}_${DATE_TAG}.tsv"

    echo "Running SNPGenie for threshold $label (minfreq=$minfreq)..."
    mkdir -p "${threshold_root}" "${out_root}"
    printf "sample\tstatus\tlog\n" > "${status_tsv}"

    tail -n +2 "${MANIFEST_FILE}" | while IFS=$'\t' read -r sample _ staged_vcf; do
        sample_work="${threshold_root}/${sample}"
        sample_out="${out_root}/${sample}"
        sample_log="${LOGS_DIR}/${label}_${sample}_${DATE_TAG}.log"

        mkdir -p "${sample_work}" "${sample_out}"
        cp -f "${staged_vcf}" "${sample_work}/variants.vcf"
        cp -f "${FASTA_FILE}" "${sample_work}/viral_only.fasta"
        cp -f "${GTF_FILE}" "${sample_work}/viral_only.gtf"

        (
            cd "${sample_work}"
            perl "${SNPGENIE_SCRIPT}" \
                --vcfformat=2 \
                --minfreq="${minfreq}" \
                --snpreport="variants.vcf" \
                --fastafile="viral_only.fasta" \
                --gtffile="viral_only.gtf"
        ) > "${sample_log}" 2>&1
        code=$?

        if [ ${code} -eq 0 ] && [ -d "${sample_work}/SNPGenie_Results" ]; then
            rm -rf "${sample_out}"
            mv "${sample_work}/SNPGenie_Results" "${sample_out}"
            printf "%s\tSUCCESS\t%s\n" "${sample}" "${sample_log}" >> "${status_tsv}"
        else
            printf "%s\tFAILED\t%s\n" "${sample}" "${sample_log}" >> "${status_tsv}"
            echo "  ✗ Failed for sample $sample (check log: $sample_log)"
        fi
    done
    echo "Completed threshold $label."
}

run_one_threshold "0.001" "minfreq_0p001"
run_one_threshold "0.01" "minfreq_0p01"

# 3. Merging and Downstream statistics
echo ""
echo "Merging and summarizing results..."
eval "$(conda shell.bash hook)"
conda activate annotation-env

# Summarize (using our custom helper)
python3 "${SCRIPT_DIR}/Helpers/summarize_snpgenie_outputs.py" \
    --output-dir "${OUTPUT_DIR}" \
    --summary-dir "${SUMMARY_DIR}"

# Analyze delta selection (Kruskal-Wallis)
python3 "${SCRIPT_DIR}/Helpers/analyze_delta_selection.py" \
    --input "${SUMMARY_DIR}/product_results_all_samples.tsv" \
    --outdir "${ANALYSIS_DIR}" \
    --manifest "${VARIANT_ROOT}/config/samples_manifest.tsv" \
    --write-detailed

# Analyze delta limma (R script)
echo "Running limma differential selection analysis in R..."
Rscript "${SCRIPT_DIR}/Helpers/analyze_delta_limma.R" \
    --delta-file="${ANALYSIS_DIR}/delta_per_sample.tsv" \
    --outdir="${ANALYSIS_DIR}" \
    --manifest="${VARIANT_ROOT}/config/samples_manifest.tsv" \
    --write-detailed

# Build compact selection tables
python3 "${SCRIPT_DIR}/Helpers/build_compact_selection_tables.py" \
    --base "${ANALYSIS_DIR}" \
    --product-summary "${SUMMARY_DIR}/product_results_all_samples.tsv" \
    --report-dir "${REPORT_DIR}" \
    --manifest "${VARIANT_ROOT}/config/samples_manifest.tsv" \
    --threshold "minfreq_0p01"

echo "═══════════════════════════════════════════════════════════"
echo "  SNPGenie selection analysis complete"
echo "  Results: $SNPGENIE_DIR"
echo "═══════════════════════════════════════════════════════════"
