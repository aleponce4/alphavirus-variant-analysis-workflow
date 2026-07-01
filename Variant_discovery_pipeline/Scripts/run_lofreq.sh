#!/bin/bash

# Multi-sample LoFreq variant calling with optional Stage-1 manifest support.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

source ./config.sh

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

usage() {
    echo "Usage: ./run_lofreq.sh [--manifest /path/to/manifest.tsv] [--status /path/to/status.tsv] [--reference /path/to/ref.fa] [--threads N]"
}

resolve_path() {
    local path="$1"
    if [ -z "${path}" ]; then
        echo ""
    elif [[ "${path}" == /* ]]; then
        echo "${path}"
    else
        echo "${REPO_DIR}/${path}"
    fi
}

resolve_reference() {
    local ref="$1"
    if [ -n "${ref}" ] && [ -f "${ref}" ]; then
        echo "${ref}"
        return 0
    fi

    while IFS= read -r -d '' candidate; do
        echo "${candidate}"
        return 0
    done < <(find "${REPO_DIR}/Input/Reference" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
}

log_status() {
    local sample_id="$1"
    local status="$2"
    local start_ts="$3"
    local end_ts="$4"
    local error="$5"
    printf "%s\tlofreq\t%s\t%s\t%s\t%s\n" "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}" >> "${STATUS_PATH}"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --manifest)
                MANIFEST_PATH="$2"
                shift 2
                ;;
            --status)
                STATUS_PATH="$2"
                shift 2
                ;;
            --reference)
                REFERENCE_OVERRIDE="$2"
                shift 2
                ;;
            --threads)
                THREADS_OVERRIDE="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                echo "ERROR: unknown option $1"
                usage
                exit 1
                ;;
        esac
    done
}

MANIFEST_PATH="${STAGE1_SAMPLE_MANIFEST:-}"
STATUS_PATH="${STAGE1_MODULE_STATUS:-LoFreq/run_status.tsv}"
REFERENCE_OVERRIDE="${REFERENCE_FASTA:-}"
THREADS_OVERRIDE="${THREADS:-4}"
parse_args "$@"
REFERENCE_FASTA="${REFERENCE_OVERRIDE}"
THREADS="${THREADS_OVERRIDE}"
export REFERENCE_FASTA
export THREADS

BAM_DIR="$(resolve_path "${BAM_OUTPUT_DIR:-Input/BAMs}")"
REFERENCE_DEFAULT=""
if [ -n "${REFERENCE_FASTA}" ]; then
    REFERENCE_DEFAULT="$(resolve_path "${REFERENCE_FASTA}")"
fi
if [ -z "${REFERENCE_DEFAULT}" ] || [ ! -f "${REFERENCE_DEFAULT}" ]; then
    while IFS= read -r -d '' candidate; do
        REFERENCE_DEFAULT="${candidate}"
        break
    done < <(find "${REPO_DIR}/Input/Reference" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
fi

mkdir -p LoFreq
mkdir -p "$(dirname "${STATUS_PATH}")"
printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${STATUS_PATH}"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not available"
    exit 1
fi

eval "$(conda shell.bash hook)"
conda activate lofreq-env 2>/dev/null \
    || conda activate ivar_env 2>/dev/null \
    || { echo "ERROR: could not activate lofreq-env or ivar_env"; exit 1; }

if ! command -v lofreq >/dev/null 2>&1; then
    echo "ERROR: lofreq command not found"
    exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools command not found"
    exit 1
fi
if ! command -v bgzip >/dev/null 2>&1; then
    echo "ERROR: bgzip command not found"
    exit 1
fi
if ! command -v tabix >/dev/null 2>&1; then
    echo "ERROR: tabix command not found"
    exit 1
fi
if ! command -v zcat >/dev/null 2>&1; then
    echo "ERROR: zcat command not found"
    exit 1
fi

declare -A SAMPLE_STATUS
declare -A SAMPLE_BAM
declare -A SAMPLE_REFERENCE
declare -A SAMPLE_CONTIG
declare -a SAMPLE_IDS

if [ -n "${MANIFEST_PATH}" ] && [ -f "${MANIFEST_PATH}" ]; then
    while IFS='' read -r line; do
        sample_id="$(printf '%s' "${line}" | cut -f1)"
        if [ "${sample_id}" = "sample_id" ] || [ -z "${sample_id}" ]; then
            continue
        fi
        if [ -n "${SAMPLE_BAM["${sample_id}"]+x}" ]; then
            continue
        fi
        sample_status="$(printf '%s' "${line}" | cut -f8)"
        sample_bam="$(printf '%s' "${line}" | cut -f3)"
        sample_reference_col="$(printf '%s' "${line}" | cut -f6)"
        sample_contig="$(printf '%s' "${line}" | cut -f7)"
        sample_status="${sample_status//$'\r'/}"
        sample_bam="${sample_bam//$'\r'/}"
        sample_reference_col="${sample_reference_col//$'\r'/}"
        sample_contig="${sample_contig//$'\r'/}"
        if [ -z "${sample_status}" ]; then
            sample_status="MANIFEST_MISSING_STATUS"
        fi
        if [ -z "${sample_bam}" ]; then
            sample_bam=""
        fi

        SAMPLE_STATUS["${sample_id}"]="${sample_status}"
        SAMPLE_BAM["${sample_id}"]="${sample_bam}"
        SAMPLE_REFERENCE["${sample_id}"]="$(resolve_path "${sample_reference_col}")"
        SAMPLE_CONTIG["${sample_id}"]="${sample_contig}"
        SAMPLE_IDS+=("${sample_id}")
    done < <(tail -n +2 "${MANIFEST_PATH}")
else
    while IFS= read -r -d '' bam_file; do
        sample_id="$(basename "${bam_file}" .bam)"
        SAMPLE_STATUS["${sample_id}"]="OK"
        SAMPLE_BAM["${sample_id}"]="${bam_file}"
        SAMPLE_IDS+=("${sample_id}")
        SAMPLE_REFERENCE["${sample_id}"]="${REFERENCE_DEFAULT}"
        SAMPLE_CONTIG["${sample_id}"]="${TARGET_CONTIG:-}"
    done < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*.bam" -print0)
fi

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "ERROR: No input samples to run LoFreq on."
    exit 1
fi

for sample_id in "${SAMPLE_IDS[@]}"; do
    rm -rf "LoFreq/${sample_id}"
done

failed_count=0

process_sample() {
    local sample_id="$1"
    local bam_file="$2"
    local input_status="$3"
    local sample_reference="$4"
    local sample_contig="$5"
    local start_ts
    local end_ts
    local status="FAILED"
    local error="unknown"
    local out_dir="LoFreq/${sample_id}"
    local target_bam="${out_dir}/target_only.bam"
    local proc_bam="${out_dir}/target_only.bam"
    local realigned_bam="${out_dir}/realigned.bam"
    local indel_bam="${out_dir}/indelqual.bam"

    start_ts="$(ts_now)"
    end_ts="${start_ts}"

    mkdir -p "${out_dir}"

    if [ "${input_status}" != "OK" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "${input_status}"
        return 0
    fi

    if [ -z "${sample_reference}" ] || [ ! -f "${sample_reference}" ]; then
        sample_reference="$(resolve_reference "${REFERENCE_DEFAULT}")"
        if [ -z "${sample_reference}" ] || [ ! -f "${sample_reference}" ]; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "missing_reference"
            return 1
        fi
    fi

    if [ -z "${sample_contig}" ]; then
        sample_contig="${TARGET_CONTIG:-}"
    fi

    if [ -z "${bam_file}" ] || [ ! -f "${bam_file}" ]; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "missing_bam"
        return 1
    fi

    if ! samtools faidx "${sample_reference}"; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "reference_index_failed"
        return 1
    fi

    if [ -n "${sample_contig}" ]; then
        if ! samtools view -b "${bam_file}" "${sample_contig}" | samtools sort -o "${target_bam}" -; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_failed"
            return 1
        fi
    else
        if ! samtools view -b "${bam_file}" | samtools sort -o "${target_bam}" -; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_failed"
            return 1
        fi
    fi
    if ! samtools index "${target_bam}"; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_index_failed"
        return 1
    fi

    if [ "${LOFREQ_BAQ:-1}" = "1" ]; then
        if ! lofreq viterbi -f "${sample_reference}" "${target_bam}" | samtools sort -o "${realigned_bam}" -; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "viterbi_failed"
            return 1
        fi
        if ! samtools index "${realigned_bam}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "realigned_bam_index_failed"
            return 1
        fi
        proc_bam="${realigned_bam}"
    fi

    if [ "${LOFREQ_ENABLE_INDELQUAL:-1}" = "1" ]; then
        if ! lofreq indelqual --dindel -f "${sample_reference}" -o "${indel_bam}" "${proc_bam}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "indelqual_failed"
            return 1
        fi
        if ! samtools index "${indel_bam}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "indel_bam_index_failed"
            return 1
        fi
        proc_bam="${indel_bam}"
    fi

    if ! lofreq call-parallel --pp-threads "${THREADS:-4}" -f "${sample_reference}" \
        --min-cov "${LOFREQ_MIN_VARIANT_DEPTH}" \
        --min-bq "${LOFREQ_MIN_BASE_QUALITY}" \
        --min-alt-bq "${LOFREQ_MIN_BASE_QUALITY}" \
        --min-mq "${LOFREQ_MIN_MAP_QUALITY}" \
        --sig "${LOFREQ_CALL_SIG}" \
        -o "${out_dir}/variants.vcf" \
        "${proc_bam}"; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "lofreq_call_failed"
        return 1
    fi

    if ! lofreq filter -i "${out_dir}/variants.vcf" -o "${out_dir}/variants.filtered.vcf" \
        --snvqual-thresh 20 --indelqual-thresh 20; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "lofreq_filter_failed"
        return 1
    fi

    if [ -f "${out_dir}/variants.filtered.vcf" ]; then
        if ! bgzip -f "${out_dir}/variants.filtered.vcf"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "bgzip_failed"
            return 1
        fi
        if ! tabix -f -p vcf "${out_dir}/variants.filtered.vcf.gz"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "tabix_failed"
            return 1
        fi
    fi

    if ! samtools view -H "${proc_bam}" >/dev/null 2>&1; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "invalid_processing_bam"
        return 1
    fi

    {
        echo "Sample: ${sample_id}" > "${out_dir}/qc_stats.txt"
        echo "Reference: ${sample_reference}" >> "${out_dir}/qc_stats.txt"
        echo "Raw variant lines: $(grep -vc '^#' "${out_dir}/variants.vcf" 2>/dev/null || echo 0)" >> "${out_dir}/qc_stats.txt"
        if [ -f "${out_dir}/variants.filtered.vcf.gz" ]; then
            echo "Filtered variant lines: $(zcat "${out_dir}/variants.filtered.vcf.gz" | grep -vc '^#' 2>/dev/null || echo 0)" >> "${out_dir}/qc_stats.txt"
        fi
    }

    if [ -f "${out_dir}/variants.vcf" ] || [ -f "${out_dir}/variants.filtered.vcf.gz" ]; then
        status="OK"
        error=""
    else
        status="FAILED"
        error="missing_lofreq_outputs"
    fi

    end_ts="$(ts_now)"
    log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}"

    if [ "${status}" != "OK" ]; then
        return 1
    fi
    return 0
}

for sample_id in "${SAMPLE_IDS[@]}"; do
    status_in="${SAMPLE_STATUS["${sample_id}"]:-OK}"
    bam_file="${SAMPLE_BAM["${sample_id}"]:-}"
    if [ -z "${bam_file}" ]; then
        bam_file="${BAM_DIR}/${sample_id}.bam"
    fi
    sample_reference="${SAMPLE_REFERENCE["${sample_id}"]:-${REFERENCE_DEFAULT}}"
    sample_contig="${SAMPLE_CONTIG["${sample_id}"]:-${TARGET_CONTIG:-}}"

    if ! process_sample "${sample_id}" "${bam_file}" "${status_in}" "${sample_reference}" "${sample_contig}"; then
        failed_count=$((failed_count + 1))
    fi
done

if [ "${failed_count}" -gt 0 ]; then
    echo "LoFreq completed with ${failed_count} failed sample(s)."
    exit 1
fi

echo "LoFreq pipeline completed."

