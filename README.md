# VEEV Variant Discovery Pipeline

Comprehensive viral variant discovery pipeline comparing iVar (frequency-based) and LoFreq (statistical) variant calling methods with functional annotation.

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

- **Two variant callers**: iVar (frequency-based) + LoFreq (statistical)
- **Quality filtering**: Base quality, coverage, mapping quality thresholds
- **Functional annotation**: Predict amino acid consequences using GFF3
- **Test mode**: Validate pipeline on subsampled data before production run
- **Per-sample runner**: `./Scripts/run_single_sample.sh <sample>`

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
- Coverage thresholds (iVar: 5000×, LoFreq: 5000×)
- Frequency threshold (0.1%)
- Quality scores (base quality 30, map quality 60)
- Threads (default 20)

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
