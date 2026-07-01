#!/bin/bash

# Multi-sample iVar variant calling with optional Stage-1 manifest support.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

log_status() {
    local sample_id="$1"
    local status="$2"
    local start_ts="$3"
    local end_ts="$4"
    local error="$5"
    printf "%s\tivar\t%s\t%s\t%s\t%s\n" "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}" >> "${STATUS_PATH}"
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

usage() {
    echo "Usage: ./run_ivar.sh [--manifest /path/to/manifest.tsv] [--status /path/to/status.tsv] [--reference /path/to/ref.fa] [--threads N]"
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

source ./config.sh
MANIFEST_PATH="${STAGE1_SAMPLE_MANIFEST:-}"
STATUS_PATH="${STAGE1_MODULE_STATUS:-Ivar/run_status.tsv}"
REFERENCE_OVERRIDE="${REFERENCE_FASTA:-}"
THREADS_OVERRIDE="${THREADS:-4}"
parse_args "$@"
REFERENCE_FASTA="${REFERENCE_OVERRIDE}"
THREADS="${THREADS_OVERRIDE}"
export REFERENCE_FASTA
export THREADS

REFERENCE_DEFAULT="$(resolve_path "${REFERENCE_FASTA}")"
if [ -z "${REFERENCE_DEFAULT}" ] || [ ! -f "${REFERENCE_DEFAULT}" ]; then
    while IFS= read -r -d '' candidate; do
        REFERENCE_DEFAULT="${candidate}"
        break
    done < <(find "${REPO_DIR}/Input/Reference" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not available"
    exit 1
fi

mkdir -p Ivar
mkdir -p "$(dirname "${STATUS_PATH}")"
printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${STATUS_PATH}"

eval "$(conda shell.bash hook)"
conda activate ivar_env 2>/dev/null \
    || { echo "ERROR: could not activate ivar_env"; exit 1; }

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools command not found"
    exit 1
fi
if ! command -v ivar >/dev/null 2>&1; then
    echo "ERROR: ivar command not found"
    exit 1
fi

BAM_DIR="$(resolve_path "${BAM_OUTPUT_DIR:-Input/BAMs}")"
PRIMER_DIR="$(resolve_path "Input/Primers")"

declare -A SAMPLE_STATUS
declare -A SAMPLE_BAM
declare -a SAMPLE_IDS
declare -A SAMPLE_REFERENCE
declare -A SAMPLE_CONTIG

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
        SAMPLE_REFERENCE["${sample_id}"]="${REFERENCE_FASTA}"
        SAMPLE_CONTIG["${sample_id}"]="${TARGET_CONTIG:-}"
    done < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*.bam" -print0)
fi

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "ERROR: No input samples to run iVar on."
    exit 1
fi

for sample_id in "${SAMPLE_IDS[@]}"; do
    sample_dir="Ivar/${sample_id}"
    mkdir -p "${sample_dir}"
done

primer_bed=$(find "${PRIMER_DIR}" -maxdepth 1 -type f -name "*.bed" 2>/dev/null | head -n 1 || true)
thread_count="${THREADS:-4}"
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
    local sample_out="Ivar/${sample_id}"
    local target_bam="${sample_out}/target_only.bam"
    local working_bam="${target_bam}"
    local qc_file="${sample_out}/qc_stats.txt"
    local variants_file="${sample_out}/variants.tsv"
    local consensus_file="${sample_out}/consensus.fa"

    start_ts="$(ts_now)"
    end_ts="${start_ts}"

    if [ "${input_status}" != "OK" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "${input_status}"
        return 0
    fi

    if [ -z "${bam_file}" ] || [ ! -f "${bam_file}" ]; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "missing_bam"
        return 1
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

    if ! samtools faidx "${sample_reference}"; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "reference_index_failed"
        return 1
    fi

    if [ -n "${sample_contig}" ]; then
        if ! samtools view -@ "${thread_count}" -b "${bam_file}" "${sample_contig}" | \
            samtools sort -@ "${thread_count}" -o "${target_bam}" - ; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_failed"
            return 1
        fi
    else
        if ! samtools view -@ "${thread_count}" -b "${bam_file}" | \
            samtools sort -@ "${thread_count}" -o "${target_bam}" - ; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_failed"
            return 1
        fi
    fi

    if ! samtools index "${target_bam}"; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "target_bam_index_failed"
        return 1
    fi

    if [ -n "${primer_bed}" ] && [ -f "${primer_bed}" ]; then
        if ! ivar trim -i "${target_bam}" -b "${primer_bed}" -p "${sample_out}/trimmed" \
            -m "${IVAR_PRIMER_TRIM_QUALITY}" -q "${IVAR_PRIMER_BASE_QUALITY}" -s "${IVAR_PRIMER_WINDOW_SIZE}" -e; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "ivar_trim_failed"
            return 1
        fi
        if ! samtools sort -@ "${thread_count}" -o "${sample_out}/trimmed.sorted.bam" "${sample_out}/trimmed.bam"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "trim_sort_failed"
            return 1
        fi
        if ! samtools index "${sample_out}/trimmed.sorted.bam"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "trim_index_failed"
            return 1
        fi
        working_bam="${sample_out}/trimmed.sorted.bam"
    fi

    if [ -n "${sample_contig}" ]; then
        if ! samtools mpileup -aa -A -d 0 -Q 0 -r "${sample_contig}" --reference "${sample_reference}" "${working_bam}" | \
            ivar consensus -p "${sample_out}/consensus" -m "${IVAR_MIN_CONSENSUS_COVERAGE}" -t 0.5 -q "${IVAR_MIN_BASE_QUALITY}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "ivar_consensus_failed"
            return 1
        fi
    else
        if ! samtools mpileup -aa -A -d 0 -Q 0 --reference "${sample_reference}" "${working_bam}" | \
            ivar consensus -p "${sample_out}/consensus" -m "${IVAR_MIN_CONSENSUS_COVERAGE}" -t 0.5 -q "${IVAR_MIN_BASE_QUALITY}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "ivar_consensus_failed"
            return 1
        fi
    fi

    if [ -n "${sample_contig}" ]; then
        if ! samtools mpileup -aa -A -d 0 -Q 0 -r "${sample_contig}" --reference "${sample_reference}" "${working_bam}" | \
            ivar variants -p "${sample_out}/variants" \
            -r "${sample_reference}" -m "${IVAR_MIN_VARIANT_DEPTH}" -t "${IVAR_MIN_VARIANT_FREQ}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "ivar_variants_failed"
            return 1
        fi
    else
        if ! samtools mpileup -aa -A -d 0 -Q 0 --reference "${sample_reference}" "${working_bam}" | \
            ivar variants -p "${sample_out}/variants" \
            -r "${sample_reference}" -m "${IVAR_MIN_VARIANT_DEPTH}" -t "${IVAR_MIN_VARIANT_FREQ}"; then
            log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "ivar_variants_failed"
            return 1
        fi
    fi

    {
        echo "Sample: ${sample_id}"
        echo "Reference length: $(grep -v '^>' "${sample_reference}" | tr -d '\n' | wc -c)"
        echo "Consensus length: $(grep -v '^>' "${consensus_file}" | tr -d '\n' | wc -c 2>/dev/null || echo 0)"
        echo "N positions: $(grep -v '^>' "${consensus_file}" 2>/dev/null | tr -d '\n' | grep -o 'N' | wc -l 2>/dev/null || echo 0)"
        echo "Variants called: $(tail -n +2 "${variants_file}" 2>/dev/null | wc -l || echo 0)"
    } > "${qc_file}"

    if [ -f "${qc_file}" ] && [ -f "${variants_file}" ]; then
        status="OK"
        error=""
    else
        status="FAILED"
        error="missing_ivar_outputs"
    fi

    end_ts="$(ts_now)"
    log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}"
    if [ "${status}" != "OK" ]; then
        return 1
    fi
    return 0
}

for sample_id in "${SAMPLE_IDS[@]}"; do
    input_status="${SAMPLE_STATUS["${sample_id}"]:-OK}"
    bam_file="${SAMPLE_BAM["${sample_id}"]:-}"
    if [ -z "${bam_file}" ]; then
        bam_file="${BAM_DIR}/${sample_id}.bam"
    fi
    sample_reference="${SAMPLE_REFERENCE["${sample_id}"]:-${REFERENCE_DEFAULT}}"
    sample_contig="${SAMPLE_CONTIG["${sample_id}"]:-${TARGET_CONTIG:-}}"
    if ! process_sample "${sample_id}" "${bam_file}" "${input_status}" "${sample_reference}" "${sample_contig}"; then
        failed_count=$((failed_count + 1))
    fi
done

if [ "${failed_count}" -gt 0 ]; then
    echo "iVar completed with ${failed_count} failed sample(s)."
    exit 1
fi

echo "iVar pipeline completed."

