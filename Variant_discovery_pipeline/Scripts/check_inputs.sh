#!/bin/bash

# Input validation + sample manifest generation for Stage 1.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${REPO_ROOT}/config.sh"

DEFAULT_MANIFEST_PATH="${STAGE1_MANIFEST:-${REPO_ROOT}/work/stage1/manifest/sample_inputs.tsv}"
MODE="${PIPELINE_INPUT_MODE:-bam}"
REFERENCE_PATH="${REFERENCE_FASTA:-}"
MANIFEST_PATH="${DEFAULT_MANIFEST_PATH}"

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

ts_now() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

usage() {
    echo "Usage: ./Scripts/check_inputs.sh [--mode {bam|fastq}] [--reference /path/to/ref.fa] [--manifest /path/to/manifest.tsv]"
    echo "  --mode       Input mode used for this run (bam or fastq)."
    echo "  --reference  Override REFERENCE_FASTA."
    echo "  --manifest   Output manifest path."
}

append_manifest_row() {
    local sample_id="$1"
    local mode="$2"
    local input_bam="$3"
    local input_r1="$4"
    local input_r2="$5"
    local reference="$6"
    local TARGET_CONTIG="$7"
    local status="$8"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" "${mode}" "${input_bam}" "${input_r1}" "${input_r2}" \
        "${reference}" "${TARGET_CONTIG}" "${status}"
}

sanitize_filename() {
    local filepath="$1"
    local filename
    local directory
    local sanitized
    local stem
    local ext
    local target
    local source_abs
    local target_abs
    local suffix=1

    if [ ! -f "${filepath}" ]; then
        echo "ERROR: input file missing: ${filepath}" >&2
        return 1
    fi

    filename="$(basename "${filepath}")"
    directory="$(dirname "${filepath}")"
    sanitized="$(printf "%s" "${filename}" | sed -E 's/[^A-Za-z0-9._-]/_/g')"

    if [ "${filename}" = "${sanitized}" ]; then
        echo "${filepath}"
        return 0
    fi

    target="${directory}/${sanitized}"

    # Ensure deterministic suffix if the sanitized target already exists.
    if [ -e "${target}" ]; then
        if command -v readlink >/dev/null 2>&1; then
            source_abs="$(readlink -f "${filepath}")"
            target_abs="$(readlink -f "${target}" 2>/dev/null || true)"
        elif command -v realpath >/dev/null 2>&1; then
            source_abs="$(realpath "${filepath}")"
            target_abs="$(realpath "${target}" 2>/dev/null || true)"
        fi

        if [ -n "${source_abs}" ] && [ -n "${target_abs}" ] && [ "${source_abs}" = "${target_abs}" ]; then
            echo "${target}"
            return 0
        fi

        stem="${sanitized}"
        ext=""
        if [[ "${sanitized}" == *.* ]]; then
            stem="${sanitized%.*}"
            ext=".${sanitized##*.}"
        fi
        while :; do
            target="${directory}/${stem}_dup${suffix}${ext}"
            if [ ! -e "${target}" ]; then
                break
            fi
            suffix=$((suffix + 1))
        done
    fi

    if ! command -v ln >/dev/null 2>&1; then
        echo "ERROR: cannot create sanitized path for '${filepath}' (ln not available)" >&2
        return 1
    fi
    if ! ln -s "${filepath}" "${target}" 2>/dev/null; then
        echo "ERROR: could not create sanitized path for '${filepath}' -> '${target}'" >&2
        return 1
    fi

    echo "INFO: sanitized filename '${filename}' -> '$(basename "${target}")'" >&2
    echo "${target}"
}

normalize_fastq_stem() {
    local filename="$1"
    local base="${filename}"
    base="${base%.gz}"
    base="${base%.fastq}"
    base="${base%.fq}"
    printf "%s" "${base}"
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
        done < <(find "${INPUT_FASTQ_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0)
    done
}

find_matching_fastq() {
    local stem="$1"
    local pattern
    local match_file=""

    for pattern in "*${stem}*.fastq.gz" "*${stem}*.fq.gz" "*${stem}*.fastq" "*${stem}*.fq"; do
        while IFS= read -r -d '' file; do
            match_file="${file}"
            break
        done < <(find "${INPUT_FASTQ_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0)

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

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --mode)
                MODE="$2"
                shift 2
                ;;
            --reference)
                REFERENCE_PATH="$2"
                shift 2
                ;;
            --manifest)
                MANIFEST_PATH="$2"
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
}

parse_args "$@"

if [[ "${MODE}" != "fastq" && "${MODE}" != "bam" ]]; then
    echo "ERROR: Invalid mode: ${MODE}"
    echo "Set --mode bam or --mode fastq"
    exit 1
fi

mkdir -p "$(dirname "${MANIFEST_PATH}")"

INPUT_FASTQ_DIR="$(resolve_path "${INPUT_FASTQ_DIR:-Input/FASTQ}")"
INPUT_BAM_DIR="$(resolve_path "${BAM_OUTPUT_DIR:-Input/BAMs}")"
REFERENCE_DIR="${REPO_ROOT}/Input/Reference"

if [ -n "${REFERENCE_PATH}" ]; then
    REFERENCE_PATH="$(resolve_path "${REFERENCE_PATH}")"
fi

if [ -z "${REFERENCE_PATH}" ] || [ ! -f "${REFERENCE_PATH}" ]; then
    while IFS= read -r -d '' ref; do
        REFERENCE_PATH="${ref}"
        break
    done < <(find "${REFERENCE_DIR}" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print0)
fi

if [ -z "${REFERENCE_PATH}" ] || [ ! -f "${REFERENCE_PATH}" ]; then
    echo "ERROR: No FASTA reference found in ${REFERENCE_DIR}."
    echo "Set REFERENCE_FASTA in config.sh or pass --reference."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools command not found."
    exit 1
fi

if [ "${MODE}" == "fastq" ]; then
    required_dirs=("${INPUT_FASTQ_DIR}")
else
    required_dirs=("${INPUT_BAM_DIR}")
fi

for required_dir in "${required_dirs[@]}"; do
    if [ ! -d "${required_dir}" ]; then
        echo "ERROR: Required directory not found: ${required_dir}"
        exit 1
    fi
done
if [ -n "${REFERENCE_PATH}" ] && [ ! -f "${REFERENCE_PATH}" ]; then
    echo "ERROR: Reference file not found: ${REFERENCE_PATH}"
    exit 1
fi
if [ -z "${REFERENCE_PATH}" ] && [ ! -d "${REFERENCE_DIR}" ]; then
    echo "ERROR: Missing reference directory: ${REFERENCE_DIR}"
    exit 1
fi

echo "Pipeline mode: ${MODE}"
echo "Reference: ${REFERENCE_PATH}"
echo "Manifest output: ${MANIFEST_PATH}"

# Header is always written.
printf "sample_id\tmode\tinput_bam\tinput_r1\tinput_r2\treference\ttarget_contig\tstatus\n" > "${MANIFEST_PATH}"

declare -A seen_samples
sample_count=0
ok_count=0
declare -A status_counts
status_counts["OK"]=0
status_counts["DUPLICATE_SKIPPED"]=0
status_counts["MISSING_INPUT"]=0
status_counts["INVALID_BAM_HEADER"]=0
status_counts["INDEX_FAILED"]=0

if [ "${MODE}" = "bam" ]; then
    while IFS= read -r -d '' bam_file_raw; do
        if ! bam_file="$(sanitize_filename "${bam_file_raw}")"; then
            status="INVALID_INPUT_PATH"
            sample_id="$(basename "${bam_file_raw}" .bam)"
            if [ -z "${sample_id}" ]; then
                sample_id="$(basename "${bam_file_raw}")"
            fi
            append_manifest_row "${sample_id}" "${MODE}" "${bam_file_raw}" "" "" "${REFERENCE_PATH}" "${TARGET_CONTIG:-}" "${status}" >> "${MANIFEST_PATH}"
            ((sample_count += 1))
            ((status_counts["${status}"] += 1))
            continue
        fi
        sample_id="$(basename "${bam_file}" .bam)"
        if [ -z "${sample_id}" ]; then
            continue
        fi

        status="OK"
        if [ -n "${seen_samples["${sample_id}"]+x}" ]; then
            status="DUPLICATE_SKIPPED"
        else
            seen_samples["${sample_id}"]=1
            bai_file="${bam_file}.bai"
            if [ ! -f "${bai_file}" ] && [ ! -f "${bam_file}.bam.bai" ]; then
                if ! samtools index "${bam_file}" 2>/dev/null; then
                    status="INDEX_FAILED"
                fi
            fi
            if [ "${status}" = "OK" ] && ! samtools view -H "${bam_file}" | head -n 1 | grep -q '^@HD'; then
                status="INVALID_BAM_HEADER"
            fi
        fi

        append_manifest_row "${sample_id}" "${MODE}" "${bam_file}" "" "" "${REFERENCE_PATH}" "${TARGET_CONTIG:-}" "${status}" >> "${MANIFEST_PATH}"
        ((sample_count += 1))
        ((status_counts["${status}"] += 1))
        if [ "${status}" = "OK" ]; then
            ok_count=$((ok_count + 1))
        fi
    done < <(find "${INPUT_BAM_DIR}" -maxdepth 1 -type f -name "*.bam" -print0)
else
    while IFS= read -r -d '' r1_file_raw; do
        if ! r1_file="$(sanitize_filename "${r1_file_raw}")"; then
            sample_id="$(basename "${r1_file_raw}")"
            status="INVALID_INPUT_PATH"
            append_manifest_row "${sample_id}" "${MODE}" "" "${r1_file_raw}" "" "${REFERENCE_PATH}" "${TARGET_CONTIG:-}" "${status}" >> "${MANIFEST_PATH}"
            ((sample_count += 1))
            ((status_counts["${status}"] += 1))
            continue
        fi
        r1_base="$(normalize_fastq_stem "$(basename "${r1_file}")")"
        sample_id=""
        pair_style="single"
        pair_stem=""

        while IFS=$'\t' read -r sample_id_tmp pair_style_tmp pair_stem_tmp; do
            sample_id="${sample_id_tmp}"
            pair_style="${pair_style_tmp}"
            pair_stem="${pair_stem_tmp}"
        done < <(infer_fastq_sample "${r1_base}")

        if [ -z "${sample_id}" ]; then
            echo "WARN: Could not parse sample id from ${r1_file}; skipping."
            continue
        fi

        if [ -n "${seen_samples["${sample_id}"]+x}" ]; then
            status="DUPLICATE_SKIPPED"
        else
            seen_samples["${sample_id}"]=1
            r2_file=""
            if [ "${pair_style}" = "paired" ] && [ -n "${pair_stem}" ]; then
                if ! r2_file="$(find_matching_fastq "${pair_stem}")"; then
                    r2_file=""
                fi
            else
                r2_file=""
            fi

            status="OK"
            if [ "${pair_style}" = "paired" ] && [ -z "${r2_file}" ]; then
                status="MISSING_INPUT"
            fi
        fi

        append_manifest_row "${sample_id}" "${MODE}" "" "${r1_file}" "${r2_file:-}" "${REFERENCE_PATH}" "${TARGET_CONTIG:-}" "${status}" >> "${MANIFEST_PATH}"
        ((sample_count += 1))
        ((status_counts["${status}"] += 1))
        if [ "${status}" = "OK" ]; then
            ok_count=$((ok_count + 1))
        fi
    done < <(collect_fastq_candidates)
fi

if [ "${sample_count}" -eq 0 ]; then
    if [ "${MODE}" = "bam" ]; then
        echo "ERROR: No BAM files found in ${INPUT_BAM_DIR}"
    else
        echo "ERROR: No FASTQ files found in ${INPUT_FASTQ_DIR}"
    fi
    exit 1
fi

echo "Input validation completed: $(ts_now)"
echo "Manifest generated: ${MANIFEST_PATH}"
echo "Summary:"
echo "  mode: ${MODE}"
echo "  total_samples: ${sample_count}"
echo "  ok_samples: ${ok_count}"
echo "  missing_or_error: $((sample_count - ok_count))"

