# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-07-29

### Added
- Initial v1.0.0 release of the `alphavirus-variant-analysis` Nextflow DSL2 workflow.
- **Scaffolding & Infrastructure**: Built complete directory layout, `nextflow.config`, parameter defaults, container pins in `conf/containers.config`, base process rules, and GitHub Actions CI workflow.
- **Core Variant Calling Modules**:
  - `SAMTOOLS_FAIDX`: Reference FASTA indexing.
  - `EXTRACT_VIRAL_BAM`: Contig extraction and BAM sorting/indexing.
  - `IVAR_VARIANTS` & `IVAR_CONSENSUS`: iVar variant calling and consensus generation.
  - `LOFREQ_CALL` & `LOFREQ_FILTER`: LoFreq variant calling, filtering, and QC statistics computation.
- **Downstream Analytics Modules**:
  - `BCFTOOLS_CSQ`: In-silico variant consequence annotation.
  - `COVERAGE_DEPTH` & `COVERAGE_SUMMARIZE`: Per-base depth extraction and windowed coverage summarization.
  - `SNPGENIE_RUN`: Vendored Perl-based SNPGenie evolutionary selection analysis with pinned commit validation (`71584c6c9a30b2c159844f210d70eb89df0f4e19`).
  - `CLIQUESNV` & `VILOCA`: Haplotype reconstruction subworkflows.
  - `REPORTING`: Consolidated multi-sample execution summary, sliding-window variant frequency plots, coverage profiles, haplotype visualizations, and multi-tab Excel workbook generation.
- **CLI Utilities (`bin/`)**: 15 standalone, containerized Python and R CLI scripts equipped with strict flag handling (`--help`, `--version`, exit 2 on invalid flags).
- **Testing & Verification**: Built end-to-end `nf-test` suite covering module unit tests, subworkflow assembly, and full pipeline integration tests on synthetic VEEV fixtures.
- **Execution Profiles**: Configured self-contained runtime profiles for `test` (Docker), `slurm` (SLURM + Singularity), and `awsbatch` (AWS Batch + Docker).
