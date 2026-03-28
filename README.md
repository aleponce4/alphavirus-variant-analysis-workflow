
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

# Validate inputs
./Scripts/check_inputs.sh

# Test execution on subsampled data
./Scripts/test_pipeline.sh

# Run full workflow
source config.sh
./Scripts/run_ivar.sh
./Scripts/run_lofreq.sh
./Scripts/annotate_all.sh
