#!/usr/bin/env bash
set -euo pipefail

echo "Running directory structure lint checks..."

MANDATORY_FILES=(
    "AGENT.md"
    "CHANGELOG.md"
    "main.nf"
    "nextflow.config"
    "nextflow_schema.json"
    "nf-test.config"
    ".gitignore"
    ".editorconfig"
    "assets/samplesheet.test.csv"
    "assets/snpgenie/PINNED_COMMIT"
    "assets/snpgenie/checksums.sha256"
    "conf/base.config"
    "conf/modules.config"
    "conf/containers.config"
    "conf/test.config"
    "conf/slurm.config"
    "conf/awsbatch.config"
    "bin/ivar_variants_to_vcf.py"
    "bin/summarize_coverage.py"
    "bin/summarize_snpgenie.py"
    "bin/analyze_delta_selection.py"
    "bin/build_selection_tables.py"
    "bin/analyze_delta_limma.R"
    "bin/generate_run_summary.py"
    "bin/plot_variants.py"
    "bin/plot_coverage.py"
    "bin/plot_haplotypes.py"
    "bin/export_excel_variants.py"
    "bin/parse_gb_to_gff3.py"
    "bin/convert_gff3_to_gtf.py"
    "modules/local/samtools/faidx/main.nf"
    "modules/local/extract_viral_bam/main.nf"
    "modules/local/ivar/variants/main.nf"
    "modules/local/ivar/consensus/main.nf"
    "modules/local/lofreq/call/main.nf"
    "modules/local/lofreq/filter/main.nf"
    "modules/local/bcftools/csq/main.nf"
    "modules/local/coverage/depth/main.nf"
    "modules/local/coverage/summarize/main.nf"
    "modules/local/snpgenie/run/main.nf"
    "modules/local/cliquesnv/main.nf"
    "modules/local/viloca/main.nf"
    "modules/local/report/run_summary/main.nf"
    "modules/local/report/variant_plots/main.nf"
    "modules/local/report/coverage_plots/main.nf"
    "modules/local/report/haplotype_plots/main.nf"
    "modules/local/report/excel_export/main.nf"
    "subworkflows/input_check/main.nf"
    "subworkflows/variant_calling/main.nf"
    "subworkflows/annotation/main.nf"
    "subworkflows/coverage_qc/main.nf"
    "subworkflows/selection/main.nf"
    "subworkflows/haplotype/main.nf"
    "subworkflows/reporting/main.nf"
    "tests/data/README.md"
    "tests/data/viral_ref.test.fasta"
    "tests/data/viral_ref.test.gff3"
    "tests/data/sampleA.test.bam"
    "tests/data/sampleA.test.bam.bai"
    "tests/data/sampleB.test.bam"
    "tests/data/sampleB.test.bam.bai"
)

ERRORS=0

for file in "${MANDATORY_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "MISSING MANDATORY FILE: $file"
        ERRORS=$((ERRORS + 1))
    fi
done

if [ "$ERRORS" -gt 0 ]; then
    echo "STRUCTURE LINT FAILED ($ERRORS missing files)"
    exit 1
fi

echo "STRUCTURE OK"
