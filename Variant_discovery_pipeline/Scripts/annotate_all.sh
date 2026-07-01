#!/bin/bash

# Annotate iVar and LoFreq outputs with bcftools csq.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

source ./config.sh

usage() {
    echo "Usage: ./annotate_all.sh [--manifest /path/to/manifest.tsv] [--status /path/to/status.tsv] [--reference /path/to/ref.fa] [--annotation /path/to/annotation.gff3]"
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
            --annotation)
                ANNOTATION_OVERRIDE="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                # preserve legacy positional reference argument
                if [ -z "${REFERENCE_OVERRIDE:+x}" ]; then
                    REFERENCE_OVERRIDE="$1"
                    shift
                    continue
                fi
                echo "ERROR: unknown option $1"
                usage
                exit 1
                ;;
        esac
    done
}

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

log_status() {
    local sample_id="$1"
    local status="$2"
    local start_ts="$3"
    local end_ts="$4"
    local error="$5"
    printf "%s\tannotation\t%s\t%s\t%s\t%s\n" "${sample_id}" "${status}" "${start_ts}" "${end_ts}" "${error}" >> "${STATUS_PATH}"
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

MANIFEST_PATH="${STAGE1_SAMPLE_MANIFEST:-}"
STATUS_PATH="${STAGE1_MODULE_STATUS:-Annotated_variants/run_status.tsv}"
REFERENCE_OVERRIDE="${REFERENCE_FASTA:-}"
ANNOTATION_OVERRIDE="${ANNOTATION_GFF:-}"
parse_args "$@"

mkdir -p Annotated_variants
mkdir -p "$(dirname "${STATUS_PATH}")"
printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${STATUS_PATH}"

if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate annotation-env 2>/dev/null \
      || { echo "ERROR: could not activate annotation-env"; exit 1; }
else
    echo "ERROR: conda not available"
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools not found."
    exit 1
fi

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python not found."
    exit 1
fi
if ! command -v zcat >/dev/null 2>&1; then
    echo "ERROR: zcat not found."
    exit 1
fi

if [ -z "${REFERENCE_OVERRIDE}" ]; then
    while IFS= read -r -d '' ref; do
        REFERENCE_OVERRIDE="${ref}"
        break
    done < <(find "${REPO_DIR}/Input/Reference" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
fi

REFERENCE="$(resolve_path "${REFERENCE_OVERRIDE}")"
ANNOTATION="$(resolve_path "${ANNOTATION_OVERRIDE}")"

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: reference not found: ${REFERENCE}"
    exit 1
fi
if [ ! -f "${ANNOTATION}" ]; then
    echo "ERROR: annotation file not found: ${ANNOTATION}"
    exit 1
fi

has_non_header_records() {
    local file="$1"
    if [ ! -f "${file}" ]; then
        return 1
    fi
    if [[ "${file}" == *.gz ]]; then
        if zcat "${file}" | grep -qv '^#'; then
            return 0
        fi
    else
        if grep -qv '^#' "${file}"; then
            return 0
        fi
    fi
    return 1
}

annotate_with_csq() {
    local input_file="$1"
    local output_file="$2"
    local ref_fasta="$3"
    local gff="$4"

    bcftools csq --local-csq -f "${ref_fasta}" -g "${gff}" "${input_file}" -o "${output_file}"
}

declare -A SAMPLE_STATUS
declare -A SAMPLE_REFERENCE
declare -A SAMPLE_I_STAT
declare -A SAMPLE_L_STAT
declare -A SAMPLE_INPUT_STATUS
declare -a SAMPLE_IDS

if [ -n "${MANIFEST_PATH}" ] && [ -f "${MANIFEST_PATH}" ]; then
    while IFS='' read -r line; do
        sample_id="$(printf '%s' "${line}" | cut -f1)"
        if [ "${sample_id}" = "sample_id" ] || [ -z "${sample_id}" ]; then
            continue
        fi
        if [ -n "${SAMPLE_STATUS["${sample_id}"]+x}" ]; then
            continue
        fi
        sample_status="$(printf '%s' "${line}" | cut -f8)"
        sample_reference="$(printf '%s' "${line}" | cut -f6)"
        sample_status="${sample_status//$'\r'/}"
        sample_reference="${sample_reference//$'\r'/}"
        if [ -z "${sample_status}" ]; then
            sample_status="MANIFEST_MISSING_STATUS"
        fi
        SAMPLE_INPUT_STATUS["${sample_id}"]="${sample_status}"
        SAMPLE_REFERENCE["${sample_id}"]="$(resolve_path "${sample_reference}")"
        SAMPLE_IDS+=("${sample_id}")
        SAMPLE_STATUS["${sample_id}"]=1
    done < <(tail -n +2 "${MANIFEST_PATH}")
else
    for sample_dir in Ivar/*/ LoFreq/*/; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        sample_id="$(basename "${sample_dir%/}")"
        if [ -n "${SAMPLE_STATUS["${sample_id}"]+x}" ]; then
            continue
        fi
        SAMPLE_IDS+=("${sample_id}")
        SAMPLE_INPUT_STATUS["${sample_id}"]="OK"
        SAMPLE_STATUS["${sample_id}"]=1
        SAMPLE_REFERENCE["${sample_id}"]="${REFERENCE}"
    done
fi

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "ERROR: no samples available for annotation"
    exit 1
fi

for sample_id in "${SAMPLE_IDS[@]}"; do
    sample_ref="${SAMPLE_REFERENCE["${sample_id}"]:-${REFERENCE}}"
    if [ -z "${sample_ref}" ] || [ ! -f "${sample_ref}" ]; then
        sample_ref="${REFERENCE}"
    fi
    SAMPLE_REFERENCE["${sample_id}"]="${sample_ref}"

    lofreq_filtered="${REPO_DIR}/LoFreq/${sample_id}/variants.filtered.vcf.gz"
    lofreq_raw="${REPO_DIR}/LoFreq/${sample_id}/variants.vcf"
    ivar_tsv="${REPO_DIR}/Ivar/${sample_id}/variants.tsv"

    if has_non_header_records "${lofreq_filtered}"; then
        SAMPLE_L_STAT["${sample_id}"]="filtered"
    elif [ -f "${lofreq_raw}" ] && has_non_header_records "${lofreq_raw}"; then
        SAMPLE_L_STAT["${sample_id}"]="raw"
    else
        SAMPLE_L_STAT["${sample_id}"]="none"
    fi

    if [ -f "${ivar_tsv}" ] && has_non_header_records "${ivar_tsv}"; then
        SAMPLE_I_STAT["${sample_id}"]="has_variants"
    else
        SAMPLE_I_STAT["${sample_id}"]="none"
    fi
done

lofreq_count=0
ivar_count=0
annotated_ok_count=0
annotated_skip_count=0
annotated_fail_count=0

LOFREQ_OUT_DIR="Annotated_variants/LoFreq"
IVAR_OUT_DIR="Annotated_variants/Ivar"
mkdir -p "${LOFREQ_OUT_DIR}" "${IVAR_OUT_DIR}"

for sample_id in "${SAMPLE_IDS[@]}"; do
    start_ts="$(ts_now)"
    end_ts="${start_ts}"
    input_status="${SAMPLE_INPUT_STATUS["${sample_id}"]:-OK}"

    if [ "${input_status}" != "OK" ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "${input_status}"
        annotated_skip_count=$((annotated_skip_count + 1))
        continue
    fi

    sample_ref="${SAMPLE_REFERENCE["${sample_id}"]}"
    sample_has_output=0
    sample_error=""

    if [ "${SAMPLE_I_STAT["${sample_id}"]}" = "has_variants" ]; then
        tsv_file="${REPO_DIR}/Ivar/${sample_id}/variants.tsv"
        temp_vcf="${REPO_DIR}/Ivar/${sample_id}/_tmp_variants_${sample_id}.vcf"
        if python Scripts/ivar_variants_to_vcf.py "${tsv_file}" "${temp_vcf}" "${sample_ref}" ; then
            if annotate_with_csq "${temp_vcf}" "${IVAR_OUT_DIR}/${sample_id}.vcf" "${sample_ref}" "${ANNOTATION}"; then
                ivar_count=$((ivar_count + 1))
                sample_has_output=1
            else
                sample_error="ivar_annotation_failed"
            fi
        else
            sample_error="ivar_vcf_conversion_failed"
        fi
        rm -f "${temp_vcf}"
    fi

    lofreq_source="none"
    if [ "${SAMPLE_L_STAT["${sample_id}"]}" = "filtered" ]; then
        lofreq_source="${REPO_DIR}/LoFreq/${sample_id}/variants.filtered.vcf.gz"
        lofreq_target="${LOFREQ_OUT_DIR}/${sample_id}_filtered.vcf"
    elif [ "${SAMPLE_L_STAT["${sample_id}"]}" = "raw" ]; then
        lofreq_source="${REPO_DIR}/LoFreq/${sample_id}/variants.vcf"
        lofreq_target="${LOFREQ_OUT_DIR}/${sample_id}.vcf"
    fi

    if [ "${lofreq_source}" != "none" ]; then
        if annotate_with_csq "${lofreq_source}" "${lofreq_target}" "${sample_ref}" "${ANNOTATION}"; then
            lofreq_count=$((lofreq_count + 1))
            sample_has_output=1
        else
            sample_error="lofreq_annotation_failed"
        fi
    fi

    if [ "${sample_has_output}" -eq 0 ]; then
        log_status "${sample_id}" "SKIPPED" "${start_ts}" "${end_ts}" "no_variants"
        annotated_skip_count=$((annotated_skip_count + 1))
        continue
    fi

    if [ -n "${sample_error}" ]; then
        end_ts="$(ts_now)"
        log_status "${sample_id}" "FAILED" "${start_ts}" "${end_ts}" "${sample_error}"
        annotated_fail_count=$((annotated_fail_count + 1))
        continue
    fi

    end_ts="$(ts_now)"
    log_status "${sample_id}" "OK" "${start_ts}" "${end_ts}" ""
    annotated_ok_count=$((annotated_ok_count + 1))
done

echo "Annotation completed."
echo "LoFreq annotated files: ${lofreq_count}"
echo "iVar annotated files: ${ivar_count}"
echo "OK samples: ${annotated_ok_count}"
echo "Skipped samples: ${annotated_skip_count}"
echo "Failed samples: ${annotated_fail_count}"

if [ "${annotated_fail_count}" -gt 0 ]; then
    exit 1
fi
