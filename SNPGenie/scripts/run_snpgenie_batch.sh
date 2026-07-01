#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR="${ROOT_DIR}/input"
REF_DIR="${INPUT_DIR}/reference"
VCF_DIR="${INPUT_DIR}/vcf_by_sample"
MANIFEST_FILE="${INPUT_DIR}/manifest/samples.tsv"
RUNS_DIR="${ROOT_DIR}/runs"
OUTPUT_DIR="${ROOT_DIR}/output"
LOGS_DIR="${ROOT_DIR}/logs"
REPO_DIR="${ROOT_DIR}/repo"

SNPGENIE_SCRIPT="${REPO_DIR}/snpgenie.pl"

if [[ ! -f "${SNPGENIE_SCRIPT}" ]]; then
  echo "Cannot find SNPGenie script at ${SNPGENIE_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${MANIFEST_FILE}" ]]; then
  echo "Manifest missing at ${MANIFEST_FILE}. Run scripts/prepare_inputs.sh first." >&2
  exit 1
fi

mkdir -p "${RUNS_DIR}" "${OUTPUT_DIR}" "${LOGS_DIR}"

FASTA_FILE="${REF_DIR}/target_only.fasta"
GTF_FILE="${REF_DIR}/target_reference.gtf"
DATE_TAG="$(date +%Y%m%d_%H%M%S)"

run_one_threshold() {
  local minfreq="$1"
  local label="$2"
  local threshold_root="${RUNS_DIR}/${label}"
  local out_root="${OUTPUT_DIR}/${label}"
  local status_tsv="${LOGS_DIR}/run_status_${label}_${DATE_TAG}.tsv"

  mkdir -p "${threshold_root}" "${out_root}"
  printf "sample\tstatus\tlog\n" > "${status_tsv}"

  tail -n +2 "${MANIFEST_FILE}" | while IFS=$'\t' read -r sample _ staged_vcf; do
    sample_work="${threshold_root}/${sample}"
    sample_out="${out_root}/${sample}"
    sample_log="${LOGS_DIR}/${label}_${sample}_${DATE_TAG}.log"

    mkdir -p "${sample_work}" "${sample_out}"
    cp -f "${staged_vcf}" "${sample_work}/variants.vcf"
    cp -f "${FASTA_FILE}" "${sample_work}/target_only.fasta"
    cp -f "${GTF_FILE}" "${sample_work}/target_reference.gtf"

    set +e
    (
      cd "${sample_work}"
      perl "${SNPGENIE_SCRIPT}" \
        --vcfformat=2 \
        --minfreq="${minfreq}" \
        --snpreport="variants.vcf" \
        --fastafile="target_only.fasta" \
        --gtffile="target_reference.gtf"
    ) > "${sample_log}" 2>&1
    code=$?
    set -e

    if [[ ${code} -eq 0 && -d "${sample_work}/SNPGenie_Results" ]]; then
      rm -rf "${sample_out}"
      mv "${sample_work}/SNPGenie_Results" "${sample_out}"
      printf "%s\tSUCCESS\t%s\n" "${sample}" "${sample_log}" >> "${status_tsv}"
    else
      printf "%s\tFAILED\t%s\n" "${sample}" "${sample_log}" >> "${status_tsv}"
    fi
  done

  echo "Completed threshold ${label}. Status: ${status_tsv}"
}

run_one_threshold "0.001" "minfreq_0p001"
run_one_threshold "0.01" "minfreq_0p01"


