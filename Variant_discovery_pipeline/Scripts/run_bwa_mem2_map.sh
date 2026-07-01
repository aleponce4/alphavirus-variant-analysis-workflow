#!/bin/bash

# Mapping script used by Stage 1 (FASTQ mode)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${REPO_ROOT}/config.sh"

MANIFEST_PATH="${STAGE1_SAMPLE_MANIFEST:-}"
STATUS_PATH="${STAGE1_MODULE_STATUS:-${REPO_ROOT}/work/stage1/mapping_status.tsv}"

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

usage() {
    echo "Usage: ./Scripts/run_bwa_mem2_map.sh [--manifest /path/to/sample_manifest.tsv] [--status /path/to/status.tsv]"
}

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
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

resolve_path() {
    local path="$1"
    if [ -z "${path}" ]; then
        echo ""
    elif [[ "${path}" == /* ]]; then
        echo "${path}"
    else
        echo "${REPO_ROOT}/${path}"
    fi
}

log_status() {
    local sample_id="$1"
    local status="$2"
    local start_ts="$3"
    local end_ts="$4"
    local error="$5"
    printf "%s\tmapping\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}" >> "${STATUS_PATH}"
}

collect_fastq_candidates() {
    local -a patterns=(
        "*_R1*.fastq.gz"
        "*_R1*.fq.gz"
        "*_R1*.fastq"
        "*_R1*.fq"
        "*_1*.fastq.gz"
        "*_1*.fq.gz"
        "*_1*.fastq"
        "*_1*.fq"
    )
    local pattern

    for pattern in "${patterns[@]}"; do
        while IFS= read -r -d '' file; do
            printf "%s\0" "${file}"
        done < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0)
    done
}

normalize_fastq_stem() {
    local filename="$1"
    local base="${filename}"
    base="${base%.gz}"
    base="${base%.fastq}"
    base="${base%.fq}"
    printf "%s" "${base}"
}

find_matching_fastq() {
    local stem="$1"
    local pattern
    local match_file=""

    for pattern in "*${stem}*.fastq.gz" "*${stem}*.fq.gz" "*${stem}*.fastq" "*${stem}*.fq"; do
        while IFS= read -r -d '' file; do
            match_file="${file}"
            break
        done < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0)
        if [ -n "${match_file}" ]; then
            printf "%s" "${match_file}"
            return 0
        fi
    done
    return 1
}

infer_fastq_sample() {
    local base="$1"
    local sample_id=""
    local pair_style="single"
    local pair_stem=""

    if [[ "${base}" == *_R1* ]]; then
        sample_id="${base%_R1*}"
        pair_stem="${base/_R1/_R2}"
        pair_style="paired"
    elif [[ "${base}" == *_1* ]]; then
        sample_id="${base%_1*}"
        pair_stem="${base/_1/_2}"
        pair_style="paired"
    else
        sample_id="${base}"
    fi

    printf "%s\t%s\t%s\n" "${sample_id}" "${pair_style}" "${pair_stem}"
}

if ! command -v bwa-mem2 >/dev/null 2>&1; then
    echo "ERROR: bwa-mem2 command not found."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools command not found."
    exit 1
fi

if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "${BWAMEM2_ENV_NAME:-bwa_mem2_env}" 2>/dev/null \
      || conda activate ivar_env 2>/dev/null \
      || conda activate lofreq-env 2>/dev/null \
      || { echo "ERROR: could not activate bwa-mem2 environment"; exit 1; }
fi

REFERENCE="${REFERENCE_FASTA:-}"
if [ -n "${REFERENCE}" ]; then
    REFERENCE="$(resolve_path "${REFERENCE}")"
fi
if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
    while IFS= read -r -d '' ref; do
        REFERENCE="${ref}"
        break
    done < <(find "${REPO_ROOT}/Input/Reference" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
fi
if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: reference FASTA not found."
    exit 1
fi

FASTQ_DIR="$(resolve_path "${INPUT_FASTQ_DIR:-Input/FASTQ}")"
BAM_DIR="$(resolve_path "${BAM_OUTPUT_DIR:-Input/BAMs}")"
mkdir -p "${BAM_DIR}"
mkdir -p "$(dirname "${STATUS_PATH}")"

if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

if [ ! -f "${REFERENCE}.amb" ] || [ ! -f "${REFERENCE}.ann" ] || [ ! -f "${REFERENCE}.pac" ] || [ ! -f "${REFERENCE}.bwt.2bit.64" ] || [ ! -f "${REFERENCE}.sa" ]; then
    bwa-mem2 index "${REFERENCE}"
fi

FORCE_REMAP="${FORCE_REMAP:-0}"
SKIP_EXISTING_BAM_MAP="${SKIP_EXISTING_BAM_MAP:-1}"

THREADS="${BWAMEM2_THREADS:-${THREADS:-4}}"
SORT_THREADS=$((THREADS > 1 ? THREADS / 2 : 1))
if [ "${SORT_THREADS}" -lt 1 ]; then
    SORT_THREADS=1
fi

read -r -a EXTRA_ARGS <<< "${BWAMEM2_EXTRA_ARGS:-}"

declare -A SAMPLE_STATUS
declare -A SAMPLE_R1
declare -A SAMPLE_R2
declare -A SAMPLE_REFERENCE

if [ -n "${MANIFEST_PATH}" ] && [ -f "${MANIFEST_PATH}" ]; then
    while IFS='' read -r line; do
        sample_id="$(printf '%s' "${line}" | cut -f1)"
        if [ -z "${sample_id}" ] || [ "${sample_id}" = "sample_id" ]; then
            continue
        fi
        sample_status="$(printf '%s' "${line}" | cut -f8)"
        input_bam="$(printf '%s' "${line}" | cut -f3)"
        input_r1="$(printf '%s' "${line}" | cut -f4)"
        input_r2="$(printf '%s' "${line}" | cut -f5)"
        reference="$(printf '%s' "${line}" | cut -f6)"
        sample_status="${sample_status//$'\r'/}"
        input_bam="${input_bam//$'\r'/}"
        input_r1="${input_r1//$'\r'/}"
        input_r2="${input_r2//$'\r'/}"
        reference="${reference//$'\r'/}"
        sample_status="${sample_status//$'\r'/}"
        if [ -z "${SAMPLE_R1["${sample_id}"]+x}" ]; then
            SAMPLE_STATUS["${sample_id}"]="${sample_status}"
            SAMPLE_R1["${sample_id}"]="${input_r1}"
            SAMPLE_R2["${sample_id}"]="${input_r2}"
            SAMPLE_REFERENCE["${sample_id}"]="$(resolve_path "${reference}")"
        fi
    done < <(tail -n +2 "${MANIFEST_PATH}")
else
    while IFS= read -r -d '' r1_file; do
        r1_base="$(normalize_fastq_stem "$(basename "${r1_file}")")"
        sample_id=""
        pair_style="single"
        pair_stem=""
        while IFS=$'\t' read -r sid pstyle pstem; do
            sample_id="${sid}"
            pair_style="${pstyle}"
            pair_stem="${pstem}"
        done < <(infer_fastq_sample "${r1_base}")

        if [ -z "${sample_id}" ]; then
            continue
        fi
        if [ -n "${SAMPLE_R1["${sample_id}"]+x}" ]; then
            continue
        fi

        if [ "${pair_style}" = "paired" ] && [ -n "${pair_stem}" ]; then
            if ! r2_file="$(find_matching_fastq "${pair_stem}")"; then
                r2_file=""
            fi
        else
            r2_file=""
        fi

        SAMPLE_STATUS["${sample_id}"]="OK"
        SAMPLE_R1["${sample_id}"]="${r1_file}"
        SAMPLE_R2["${sample_id}"]="${r2_file:-}"
        SAMPLE_REFERENCE["${sample_id}"]="${REFERENCE}"
    done < <(collect_fastq_candidates)
fi

if [ "${#SAMPLE_STATUS[@]}" -eq 0 ]; then
    echo "ERROR: No FASTQ samples discovered in ${FASTQ_DIR}"
    exit 1
fi

mapfile -t sample_ids < <(printf '%s\n' "${!SAMPLE_STATUS[@]}" | sort)

printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${STATUS_PATH}"

mapped_count=0
failed_count=0
skipped_count=0

for sample_id in "${sample_ids[@]}"; do
    start_ts="$(ts_now)"
    end_ts="${start_ts}"
    sample_reference="${SAMPLE_REFERENCE["${sample_id}"]:-${REFERENCE}}"
    if [ -z "${sample_reference}" ] || [ ! -f "${sample_reference}" ]; then
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "missing_reference"
        failed_count=$((failed_count + 1))
        continue
    fi

    status_in="${SAMPLE_STATUS["${sample_id}"]}"
    if [ "${status_in}" != "OK" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "${status_in}"
        skipped_count=$((skipped_count + 1))
        continue
    fi

    r1_file="${SAMPLE_R1["${sample_id}"]:-}"
    r2_file="${SAMPLE_R2["${sample_id}"]:-}"
    if [ -z "${r1_file}" ] || [ ! -f "${r1_file}" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "missing_r1"
        skipped_count=$((skipped_count + 1))
        continue
    fi

    bam_out="${BAM_DIR}/${sample_id}.bam"
    if [ -f "${bam_out}" ] && [ "${FORCE_REMAP}" != "1" ] && [ "${SKIP_EXISTING_BAM_MAP}" = "1" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "existing_bam"
        skipped_count=$((skipped_count + 1))
        continue
    fi

    if [ -f "${bam_out}" ] && [ "${FORCE_REMAP}" = "1" ]; then
        rm -f "${bam_out}" "${bam_out}.bai" "${bam_out%.bam}.tmp" 2>/dev/null || true
    fi

    if [ -n "${r2_file}" ] && [ -f "${r2_file}" ]; then
        echo "Mapping paired sample: ${sample_id}"
        if bwa-mem2 mem -t "${THREADS}" "${EXTRA_ARGS[@]}" \
            -R $'@RG\tID:'"${sample_id}"$'\tSM:'"${sample_id}"$'\tPL:ILLUMINA' \
            "${sample_reference}" "${r1_file}" "${r2_file}" \
            | samtools sort -@ "${SORT_THREADS}" -o "${bam_out}" - ; then
            if samtools index "${bam_out}"; then
                status="OK"
                mapped_count=$((mapped_count + 1))
                end_ts="$(ts_now)"
                log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" ""
            else
                status="FAILED"
                failed_count=$((failed_count + 1))
                end_ts="$(ts_now)"
                log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "samtools_index_failed"
            fi
        else
            status="FAILED"
            failed_count=$((failed_count + 1))
            end_ts="$(ts_now)"
            log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "bwa_mem2_failed"
        fi
    else
        echo "Mapping single sample: ${sample_id}"
        if bwa-mem2 mem -t "${THREADS}" "${EXTRA_ARGS[@]}" \
            -R $'@RG\tID:'"${sample_id}"$'\tSM:'"${sample_id}"$'\tPL:ILLUMINA' \
            "${sample_reference}" "${r1_file}" \
            | samtools sort -@ "${SORT_THREADS}" -o "${bam_out}" - ; then
            if samtools index "${bam_out}"; then
                status="OK"
                mapped_count=$((mapped_count + 1))
                end_ts="$(ts_now)"
                log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" ""
            else
                status="FAILED"
                failed_count=$((failed_count + 1))
                end_ts="$(ts_now)"
                log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "samtools_index_failed"
            fi
        else
            status="FAILED"
            failed_count=$((failed_count + 1))
            end_ts="$(ts_now)"
            log_status "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "bwa_mem2_failed"
        fi
    fi
done

echo "Mapping complete."
echo "  total_samples: ${#sample_ids[@]}"
echo "  mapped: ${mapped_count}"
echo "  skipped: ${skipped_count}"
echo "  failed: ${failed_count}"

if [ "${failed_count}" -gt 0 ]; then
    exit 1
fi
