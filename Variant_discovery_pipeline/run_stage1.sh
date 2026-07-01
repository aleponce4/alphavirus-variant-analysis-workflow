#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

source ./config.sh

usage() {
    echo "Usage: ./run_stage1.sh [--mode {bam|fastq}] [--reference /path/to/ref.fa] [--manifest /path/to/manifest.tsv] [--threads N] [--target-contig CONTIG]"
}

MODE="${PIPELINE_INPUT_MODE:-bam}"
REFERENCE_OVERRIDE=""
MANIFEST_PATH="${STAGE1_MANIFEST}"
THREADS_OVERRIDE="${THREADS:-4}"
TARGET_CONTIG_OVERRIDE=""
INPUT_FASTQ_DIR_OVERRIDE="${INPUT_FASTQ_DIR:-}"
BAM_OUTPUT_DIR_OVERRIDE="${BAM_OUTPUT_DIR:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="$2"
            shift 2
            ;;
        --reference)
            REFERENCE_OVERRIDE="$2"
            shift 2
            ;;
        --manifest)
            MANIFEST_PATH="$2"
            shift 2
            ;;
        --threads)
            THREADS_OVERRIDE="$2"
            shift 2
            ;;
        --target-contig|--Target-contig)
            TARGET_CONTIG_OVERRIDE="$2"
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

if [[ "${MODE}" != "bam" && "${MODE}" != "fastq" ]]; then
    echo "ERROR: invalid mode '${MODE}'"
    exit 1
fi

PIPELINE_INPUT_MODE="${MODE}"
REFERENCE_FASTA="${REFERENCE_OVERRIDE:-${REFERENCE_FASTA}}"
THREADS="${THREADS_OVERRIDE}"
if [ -n "${INPUT_FASTQ_DIR_OVERRIDE}" ]; then
    INPUT_FASTQ_DIR="${INPUT_FASTQ_DIR_OVERRIDE}"
fi
if [ -n "${BAM_OUTPUT_DIR_OVERRIDE}" ]; then
    BAM_OUTPUT_DIR="${BAM_OUTPUT_DIR_OVERRIDE}"
fi
if [ -n "${TARGET_CONTIG_OVERRIDE}" ]; then
    TARGET_CONTIG="${TARGET_CONTIG_OVERRIDE}"
fi

export PIPELINE_INPUT_MODE
export REFERENCE_FASTA
export THREADS
export INPUT_FASTQ_DIR
export BAM_OUTPUT_DIR
export TARGET_CONTIG

abs_path() {
    local path="$1"
    if [[ "${path}" == /* ]]; then
        echo "${path}"
    else
        echo "${SCRIPT_DIR}/${path}"
    fi
}

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

run_stage1_dir="$(abs_path "${STAGE1_DIR}")"
mkdir -p "${run_stage1_dir}"
RUN_STATUS_FILE="${run_stage1_dir}/run_status.tsv"
RUN_REPORT_FILE="${run_stage1_dir}/stage1_report.tsv"
MANIFEST_PATH="$(abs_path "${MANIFEST_PATH}")"
export STAGE1_MANIFEST="${MANIFEST_PATH}"

printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${RUN_STATUS_FILE}"
printf "metric\tvalue\n" > "${RUN_REPORT_FILE}"

log_metric() {
    local key="$1"
    local value="$2"
    printf "%s\t%s\n" "${key}" "${value}" >> "${RUN_REPORT_FILE}"
}

log_stage_row() {
    local sample_id="$1"
    local stage="$2"
    local status="$3"
    local start_ts="$4"
    local end_ts="$5"
    local error="$6"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${sample_id}" "${stage}" "${status}" "${start_ts}" "${end_ts}" "${error}" >> "${RUN_STATUS_FILE}"
}

read_manifest() {
    declare -g -a SAMPLE_IDS=()
    declare -g -A SAMPLE_INPUT_STATUS=()
    declare -g -A SAMPLE_REFERENCE=()
    declare -g -A SAMPLE_CONTIG=()

    while IFS='' read -r line; do
        sample_id="$(printf '%s' "${line}" | cut -f1)"
        sample_id="${sample_id//$'\r'/}"
        if [ "${sample_id}" = "sample_id" ] || [ -z "${sample_id}" ]; then
            continue
        fi

        sample_status="$(printf '%s' "${line}" | cut -f8)"
        sample_reference="$(printf '%s' "${line}" | cut -f6)"
        sample_contig="$(printf '%s' "${line}" | cut -f7)"
        sample_status="${sample_status//$'\r'/}"
        sample_reference="${sample_reference//$'\r'/}"
        sample_contig="${sample_contig//$'\r'/}"

        if [ -z "${SAMPLE_INPUT_STATUS["${sample_id}"]+x}" ]; then
            SAMPLE_IDS+=("${sample_id}")
            SAMPLE_INPUT_STATUS["${sample_id}"]="${sample_status}"
            SAMPLE_REFERENCE["${sample_id}"]="${sample_reference}"
            SAMPLE_CONTIG["${sample_id}"]="${sample_contig}"
        else
            if [ "${SAMPLE_INPUT_STATUS["${sample_id}"]}" = "OK" ] && [ "${sample_status}" != "OK" ]; then
                SAMPLE_INPUT_STATUS["${sample_id}"]="${sample_status}"
            fi
        fi
    done < <(tail -n +2 "${MANIFEST_PATH}")
}

parse_module_status() {
    local file="$1"
    local -n out_status="$2"
    local -n out_error="$3"
    local -n out_start="$4"
    local -n out_end="$5"
    local line
    local line_sample_id
    local line_status
    local line_start_ts
    local line_end_ts
    local line_error
    if [ ! -f "${file}" ]; then
        return
    fi
    while IFS='' read -r line; do
        if [ -z "${line}" ]; then
            continue
        fi

        line_sample_id="$(printf '%s' "${line}" | cut -f1)"
        line_status="$(printf '%s' "${line}" | cut -f3)"
        line_start_ts="$(printf '%s' "${line}" | cut -f4)"
        line_end_ts="$(printf '%s' "${line}" | cut -f5)"
        line_error="$(printf '%s' "${line}" | cut -f6)"

        line_sample_id="${line_sample_id//$'\r'/}"
        line_status="${line_status//$'\r'/}"
        line_start_ts="${line_start_ts//$'\r'/}"
        line_end_ts="${line_end_ts//$'\r'/}"
        line_error="${line_error//$'\r'/}"

        if [ "${line_sample_id}" = "sample_id" ] || [ -z "${line_sample_id}" ]; then
            continue
        fi

        out_status["${line_sample_id}"]="${line_status}"
        out_error["${line_sample_id}"]="${line_error}"
        out_start["${line_sample_id}"]="${line_start_ts}"
        out_end["${line_sample_id}"]="${line_end_ts}"
    done < "${file}"
}

run_stage() {
    local stage_name="$1"
    local script_path="$2"
    local status_file="$3"
    local -n out_failed="$4"
    local -n out_ok="$5"
    local -n out_skipped="$6"
    shift 6
    local -a module_args=("$@")
    local start_ts
    local end_ts
    local rc=0

    rm -f "${status_file}"
    mkdir -p "$(dirname "${status_file}")"

    export STAGE1_SAMPLE_MANIFEST="${MANIFEST_PATH}"
    export STAGE1_MODULE_STATUS="${status_file}"
    export FORCE_REMAP="${FORCE_REMAP:-0}"
    export SKIP_EXISTING_BAM_MAP="${SKIP_EXISTING_BAM_MAP:-1}"
    export BWAMEM2_EXTRA_ARGS="${BWAMEM2_EXTRA_ARGS:-}"
    export BWAMEM2_ENV_NAME="${BWAMEM2_ENV_NAME:-bwa_mem2_env}"
    export BWAMEM2_THREADS="${THREADS}"

    start_ts="$(ts_now)"
    if [ "${stage_name}" = "mapping" ] && [ "${PIPELINE_INPUT_MODE}" = "bam" ]; then
        # Explicit BAM mode behavior: mapping should be skipped for every sample.
        printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${status_file}"
        for sample_id in "${SAMPLE_IDS[@]}"; do
            input_status="${SAMPLE_INPUT_STATUS["${sample_id}"]-mode_bam}"
            if [ "${input_status}" = "OK" ] || [ -z "${input_status}" ]; then
                sample_reason="mode_bam"
            else
                sample_reason="${input_status}"
            fi

            if [ -z "${sample_reason}" ]; then
                sample_reason="mode_bam"
            fi
            printf "%s\tmapping\tSKIPPED\t%s\t%s\t%s\n" "${sample_id}" "${start_ts}" "${start_ts}" "${sample_reason}" >> "${status_file}"
        done
    else
        if bash "${script_path}" "${module_args[@]}"; then
            rc=0
        else
            rc=$?
        fi
    fi
    end_ts="$(ts_now)"

    if [ "${stage_name}" = "mapping" ] && [ "${PIPELINE_INPUT_MODE}" = "bam" ]; then
        : 
    elif [ ! -f "${status_file}" ] || [ ! -s "${status_file}" ]; then
        # ensure each sample gets a row even if the script crashes before writing a status file
        printf "sample_id\tstage\tstatus\tstart_ts\tend_ts\terror\n" > "${status_file}"
        for sample_id in "${SAMPLE_IDS[@]}"; do
            printf "%s\t${stage_name}\tFAILED\t%s\t%s\tmissing_module_status\n" "${sample_id}" "${start_ts}" "${end_ts}" >> "${status_file}"
        done
    fi

    declare -A MODULE_STATUS
    declare -A MODULE_ERROR
    declare -A MODULE_START
    declare -A MODULE_END
    parse_module_status "${status_file}" MODULE_STATUS MODULE_ERROR MODULE_START MODULE_END

    out_failed=0
    out_ok=0
    out_skipped=0

    local sample_id
    local sample_status
    local sample_error
    local sample_start
    local sample_end

    for sample_id in "${SAMPLE_IDS[@]}"; do
        sample_status="FAILED"
        sample_error="module_missing_status"
        sample_start="${start_ts}"
        sample_end="${end_ts}"
        if [ "${SAMPLE_INPUT_STATUS["${sample_id}"]}" != "OK" ]; then
            sample_status="SKIPPED"
            sample_error="${SAMPLE_INPUT_STATUS["${sample_id}"]}"
        elif [ -n "${MODULE_STATUS["${sample_id}"]+x}" ]; then
            sample_status="${MODULE_STATUS["${sample_id}"]}"
            sample_error="${MODULE_ERROR["${sample_id}"]-}"
            if [ -n "${MODULE_START["${sample_id}"]+x}" ]; then
                sample_start="${MODULE_START["${sample_id}"]}"
            fi
            if [ -n "${MODULE_END["${sample_id}"]+x}" ]; then
                sample_end="${MODULE_END["${sample_id}"]}"
            fi
            if [ "${sample_status}" = "OK" ]; then
                out_ok=$((out_ok + 1))
            elif [ "${sample_status}" = "SKIPPED" ]; then
                out_skipped=$((out_skipped + 1))
            elif [ "${sample_status}" = "FAILED" ]; then
                out_failed=$((out_failed + 1))
            else
                sample_status="FAILED"
                sample_error="${sample_error:-invalid_stage_status}"
                out_failed=$((out_failed + 1))
            fi
        else
            # If module returns no status for a sample, treat as skip if there was no module error.
            if [ "${rc}" -ne 0 ]; then
                sample_status="FAILED"
                sample_error="module_failed_exit_${rc}"
                out_failed=$((out_failed + 1))
            else
                sample_status="SKIPPED"
                sample_error="module_status_missing"
                out_skipped=$((out_skipped + 1))
            fi
        fi
        if [ "${sample_start}" = "" ]; then
            sample_start="${start_ts}"
        fi
        if [ "${sample_end}" = "" ]; then
            sample_end="${end_ts}"
        fi
        log_stage_row "${sample_id}" "${stage_name}" "${sample_status}" "${sample_start}" "${sample_end}" "${sample_error}"
    done

    log_metric "${stage_name}_ok" "${out_ok}"
    log_metric "${stage_name}_skipped" "${out_skipped}"
    log_metric "${stage_name}_failed" "${out_failed}"
    log_metric "${stage_name}_exit" "${rc}"

    return "${rc}"
}

run_start="$(ts_now)"
log_metric "started_at" "${run_start}"
log_metric "mode" "${MODE}"
log_metric "threads" "${THREADS}"
log_metric "TARGET_CONTIG" "${TARGET_CONTIG}"
log_metric "reference" "${REFERENCE_FASTA}"
log_metric "manifest" "${MANIFEST_PATH}"
log_metric "stage1_dir" "${run_stage1_dir}"

echo "Stage-1 start: ${run_start}"
echo "Mode: ${MODE}"
echo "Manifest: ${MANIFEST_PATH}"

bash "${SCRIPT_DIR}/Scripts/check_inputs.sh" \
    --mode "${MODE}" \
    --reference "${REFERENCE_FASTA}" \
    --manifest "${MANIFEST_PATH}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "ERROR: manifest not generated: ${MANIFEST_PATH}"
    exit 1
fi

read_manifest

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "ERROR: manifest is empty: ${MANIFEST_PATH}"
    exit 1
fi

manifest_ok=0
manifest_not_ok=0
for sample_id in "${SAMPLE_IDS[@]}"; do
    if [ "${SAMPLE_INPUT_STATUS["${sample_id}"]}" = "OK" ]; then
        manifest_ok=$((manifest_ok + 1))
    else
        manifest_not_ok=$((manifest_not_ok + 1))
    fi
done
log_metric "manifest_ok" "${manifest_ok}"
log_metric "manifest_not_ok" "${manifest_not_ok}"

mapping_status_file="${run_stage1_dir}/mapping_status.tsv"
run_stage "mapping" "${SCRIPT_DIR}/Scripts/run_bwa_mem2_map.sh" "${mapping_status_file}" MAPPING_FAILED MAPPING_OK MAPPING_SKIPPED --manifest "${MANIFEST_PATH}" --status "${mapping_status_file}" || true

ivar_status_file="${run_stage1_dir}/ivar_status.tsv"
run_stage "ivar" "${SCRIPT_DIR}/Scripts/run_ivar.sh" "${ivar_status_file}" IVAR_FAILED IVAR_OK IVAR_SKIPPED --manifest "${MANIFEST_PATH}" --status "${ivar_status_file}" || true

lofreq_status_file="${run_stage1_dir}/lofreq_status.tsv"
run_stage "lofreq" "${SCRIPT_DIR}/Scripts/run_lofreq.sh" "${lofreq_status_file}" LOFREQ_FAILED LOFREQ_OK LOFREQ_SKIPPED --manifest "${MANIFEST_PATH}" --status "${lofreq_status_file}" || true

annotation_status_file="${run_stage1_dir}/annotation_status.tsv"
run_stage "annotation" "${SCRIPT_DIR}/Scripts/annotate_all.sh" "${annotation_status_file}" ANNOTATION_FAILED ANNOTATION_OK ANNOTATION_SKIPPED --manifest "${MANIFEST_PATH}" --status "${annotation_status_file}" || true

run_end="$(ts_now)"
log_metric "completed_at" "${run_end}"

total_failed=$((MAPPING_FAILED + IVAR_FAILED + LOFREQ_FAILED + ANNOTATION_FAILED))
log_metric "total_failed_samples" "${total_failed}"

if [ "${total_failed}" -gt 0 ]; then
    log_metric "overall_status" "FAILED"
    echo "FAILED: ${total_failed} sample-level failures."
    echo "Run status: ${RUN_STATUS_FILE}"
    echo "Run report: ${RUN_REPORT_FILE}"
    exit 1
fi

log_metric "overall_status" "SUCCESS"
echo "SUCCESS"
echo "Run status: ${RUN_STATUS_FILE}"
echo "Run report: ${RUN_REPORT_FILE}"

