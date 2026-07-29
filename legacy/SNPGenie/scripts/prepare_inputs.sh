#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"

REF_SRC_DIR="${PROJECT_ROOT}/Variant_discovery_pipeline/Input/Reference"
LOFREQ_DIR="${PROJECT_ROOT}/Variant_discovery_pipeline/LoFreq"

INPUT_DIR="${ROOT_DIR}/input"
REF_DIR="${INPUT_DIR}/reference"
VCF_DIR="${INPUT_DIR}/vcf_by_sample"
MANIFEST_DIR="${INPUT_DIR}/manifest"

mkdir -p "${REF_DIR}" "${VCF_DIR}" "${MANIFEST_DIR}"

FASTA_SRC="${REF_SRC_DIR}/viral_only.fasta"
GFF3_SRC="${REF_SRC_DIR}/VEEV_INH_fromGenbank.gff3"
GTF_OUT="${REF_DIR}/VEEV_INH_fromGenbank.gtf"

cp -f "${FASTA_SRC}" "${REF_DIR}/viral_only.fasta"
cp -f "${GFF3_SRC}" "${REF_DIR}/VEEV_INH_fromGenbank.gff3"

python3 "${ROOT_DIR}/scripts/convert_gff3_to_gtf.py" \
  --gff3 "${REF_DIR}/VEEV_INH_fromGenbank.gff3" \
  --gtf "${GTF_OUT}"

MANIFEST_FILE="${MANIFEST_DIR}/samples.tsv"
printf "sample\tsource_filtered_vcf\tstaged_vcf\n" > "${MANIFEST_FILE}"

sample_count=0
for sample_dir in "${LOFREQ_DIR}"/INH_*; do
  [[ -d "${sample_dir}" ]] || continue
  sample="$(basename "${sample_dir}")"
  filtered_vcf_gz="${sample_dir}/variants.filtered.vcf.gz"
  staged_vcf="${VCF_DIR}/${sample}.vcf"

  if [[ ! -f "${filtered_vcf_gz}" ]]; then
    echo "Missing filtered VCF: ${filtered_vcf_gz}" >&2
    exit 1
  fi

  gzip -dc "${filtered_vcf_gz}" > "${staged_vcf}"

  if ! grep -q '^#CHROM' "${staged_vcf}"; then
    echo "Invalid VCF header in ${staged_vcf}" >&2
    exit 1
  fi

  printf "%s\t%s\t%s\n" "${sample}" "${filtered_vcf_gz}" "${staged_vcf}" >> "${MANIFEST_FILE}"
  ((sample_count += 1))
done

echo "Prepared ${sample_count} samples."
echo "Manifest: ${MANIFEST_FILE}"
