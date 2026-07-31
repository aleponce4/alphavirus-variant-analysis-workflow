# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-07-29

### Fixed
- **Scientific Audit Adjustments (RA-1 to RA-12)**:
  - **iVar Frequency Floor (RA-1)**: Lowered `ivar_min_freq` default from `0.01` (1%) to `0.001` (0.1%) to align reporting with 0.1%–1% iSNV discovery targets.
  - **Coverage Inclusion Floors (RA-2)**: Lowered inclusion depth floors `ivar_min_depth` and `lofreq_min_depth` from `1000` to `10` to eliminate artificial dead zones at amplicon shoulders and prevent downward bias in SNPGenie diversity denominators.
  - **Dead Parameter Cleanup (RA-3)**: Removed dead `lofreq_min_freq` parameter from `nextflow.config` and schema.
  - **LoFreq Mapping Quality Floor (RA-4)**: Relaxed `lofreq_min_mq` from `60` to `20` to avoid stripping informative reads from naturally variable/repeated viral regions.
  - **Haplotype Parameter Explicit Geometry (RA-5 & RA-6)**: Added explicit parameters `viloca_window` (`150`), `viloca_shift` (`50`), and `cliquesnv_min_freq` (`0.001`), wiring them into process invocations and recording them in `methods_key_parameters.tsv`.
  - **Reporting Container & Hard-Fail Statistical Guard (RA-7)**: Updated `container_python_reporting` to a multi-package data science biocontainer with pandas and scipy, and updated `analyze_delta_selection.py` to hard-fail (`exit 1`) if scipy is unavailable.
  - **iVar VCF Output Integrity (RA-8)**: Corrected fabricated fields in `ivar_variants_to_vcf.py` by emitting haploid GT (`1`) and converting binomial p-values to Phred-scaled QUAL (`-10*log10(pval)`).
  - **Proportion to Percent Reporting Contract (RA-9)**: Standardized raw machine files to proportions (0–1) and human-facing summary tables/plots to percentages (`100 × proportion`) with explicit `%` labels.
  - **Frequency Tiering Policy (RA-10)**: Documented 3-tier iSNV classification policy (`≥1%` high-confidence candidate, `0.1–1%` low-frequency candidate, `<0.1%` exploratory) in documentation and executive reporting.
  - **Explicit SNPGenie Frequency Floor (RA-11)**: Passed `--minfreq=0` explicitly to `snpgenie.pl` in `SNPGENIE_RUN`.
  - **Indel Scope Documentation (RA-12)**: Documented SNV-only indel scope under default `lofreq_enable_indelqual = false`.
  - **VILOCA Module Fix & Haplotype Reporting Integration**: Updated VILOCA module to execute real `viloca run` using pinned `quay.io/biocontainers/viloca:1.1.1--py310h563914a_0` container; split HAPLOTYPE subworkflow emits into dedicated CliqueSNV and VILOCA channels; added `build_haplotype_tables.py`, `summarize_linked_mutations.py`, and `plot_haplotypes.py` (composition stacked bar plot + Minimum Spanning Network with treatment pie nodes); integrated reporting into `REPORTING` subworkflow and `assets/executive_report.qmd`.

## [1.2.0] - 2026-07-29

### Added
- **Executive HTML Reporting**: Added self-contained HTML executive summary report module (`EXECUTIVE_REPORT`) rendering key metrics, depth profiles, variant density, selection statistics, and software provenance in a single report artifact.
- **Public & Hermetic Test Datasets**: Standardized `-profile test` for zero-friction demo execution, while retaining offline hermetic `nf-test` suite.

## [1.1.0] - 2026-07-29

### Added
- **Production Hardening & Schema Validation**: Integrated `plugin/nf-schema@2.1.1` for JSON Schema Draft 2020-12 input validation and fast parameter checking.
- **FASTQ Ingress & Read Preprocessing**: Added `READ_PREPROCESSING` subworkflow (`FASTQC`, `BWA_INDEX`, `BWA_MEM`, `SAMTOOLS_STATS`, `IVAR_TRIM`) for raw FASTQ input handling.
- **Subworkflow Architecture & Modularization**: Relocated `selection` and `haplotype` subworkflows to `subworkflows/local/` and wired explicit treatment group channels.
- **Software Provenance & MultiQC**: Added `DUMP_SOFTWARE_VERSIONS` module and `MULTIQC` terminal dashboard process.

## [1.0.0] - 2026-07-29

### Added
- Initial release of the `alphavirus-variant-analysis` Nextflow DSL2 workflow.
