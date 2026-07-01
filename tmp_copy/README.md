
# Alphavirus Intra-host Variant Analysis Workflow

Variant calling workflow used for analysis of intra-host viral diversity in deeply sequenced VEEV datasets. This repository provides a standardized execution framework for iVar and LoFreq variant calling with downstream functional annotation and optional haplotype reconstruction.

The workflow was developed to support consistent processing of viral sequencing datasets and facilitate comparison of variant calls across samples and timepoints.

## Overview

The workflow performs the following steps:

1. Input validation and preprocessing checks
2. Variant calling with iVar
3. Variant calling with LoFreq
4. Functional annotation of variants
5. Optional haplotype reconstruction analyses
6. Preparation of downstream analysis inputs

The implementation focuses on reproducible execution and consistent filtering across datasets rather than development of new variant calling methods.

## Quick start

```bash
cd Variant_discovery_pipeline

# Stage-1 orchestrator (preferred)
source config.sh
./run_stage1.sh --mode bam --manifest work/stage1/manifest/sample_inputs.tsv

# FASTQ mode
./run_stage1.sh --mode fastq --reference Input/Reference/inh.fasta --threads 20 --viral-contig VEEV_INH

# Legacy module flow (compatible)
./Scripts/check_inputs.sh
./Scripts/run_ivar.sh
./Scripts/run_lofreq.sh
./Scripts/annotate_all.sh
```

### Stage-1 contract artifacts

- `work/stage1/manifest/sample_inputs.tsv` (input contract)
- `work/stage1/run_status.tsv` (per-sample, per-stage status)
- `work/stage1/stage1_report.tsv` (run summary metrics)

`sample_id, stage, status, start_ts, end_ts, error` are the run-status columns.
