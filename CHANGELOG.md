# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
