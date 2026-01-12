# VEEV Variant Discovery Pipeline

Comprehensive viral variant discovery pipeline optimized for **detecting rare, low-frequency variants** in deeply sequenced samples. Compares iVar (frequency-based) and LoFreq (statistical) variant calling methods with functional annotation.

## Purpose

This pipeline is specifically designed to:
- **Detect ultra-rare variants** at frequencies as low as **0.1%** (1 in 1000 reads)
- **Leverage deep coverage** (typically 5000×+) to identify variants with high confidence
- **Compare two statistical approaches** (iVar vs LoFreq) to validate variant calls
- **Annotate functional consequences** to assess clinical/biological impact

This is particularly useful for viral samples with deep sequencing where rare variants may represent:
- Intra-host viral evolution
- Treatment-resistant mutations
- Minor viral populations in mixed infections

## Quick Start

```bash
cd Variant_discovery_pipeline

# Validate inputs
./Scripts/check_inputs.sh

# Test on subsampled sample (2-3 min)
./Scripts/test_pipeline.sh

# Run full pipeline
source config.sh
./Scripts/run_ivar.sh      # iVar variant calling
./Scripts/run_lofreq.sh    # LoFreq variant calling
./Scripts/annotate_all.sh  # Annotate all variants
```

## Documentation

See [Variant_discovery_pipeline/README.md](Variant_discovery_pipeline/README.md) for full documentation.

## Key Features

- **Two variant callers**: iVar (frequency-based) + LoFreq (statistical) for robust variant detection
- **Rare variant optimized**: Detects variants at 0.1% frequency in deeply covered samples (5000×+)
- **Quality filtering**: Stringent thresholds (Q30 base, Q60 mapping) ensure high-confidence calls
- **Functional annotation**: Predict amino acid consequences using GFF3 reference annotation
- **Test mode**: Validate pipeline on subsampled data before production run
- **Per-sample runner**: `./Scripts/run_single_sample.sh <sample>` for quick re-processing

## Output Structure

```
Ivar/                    → iVar variant calls (TSV)
LoFreq/                  → LoFreq variant calls (VCF)
Annotated_variants/      → Functional annotations
├── Ivar/                → Annotated iVar VCFs
└── LoFreq/              → Annotated LoFreq VCFs
```

## Configuration

Edit `config.sh` to adjust:
- **Coverage thresholds**: iVar 5000×, LoFreq 5000× (tuned for rare variant detection)
- **Frequency threshold**: 0.1% (detects 1 variant per 1000 reads)
- **Quality scores**: Base quality 30, map quality 60 (ensures high confidence)
- **Threads**: Default 20 for parallel processing

These defaults are optimized for deep sequencing of viral genomes. Adjust coverage thresholds lower only if using shallower sequencing.

## Quick Commands

| Task | Command |
|------|---------|
| Test pipeline | `./Scripts/test_pipeline.sh` |
| Run one sample | `./Scripts/run_single_sample.sh <sample>` |
| Run all samples (iVar) | `source config.sh && ./Scripts/run_ivar.sh` |
| Run all samples (LoFreq) | `source config.sh && ./Scripts/run_lofreq.sh` |
| Annotate results | `./Scripts/annotate_all.sh` |

## Requirements

- samtools, bcftools
- iVar, LoFreq
- Python 3, pandas, numpy
- Reference genome (FASTA) + GFF3 annotation
