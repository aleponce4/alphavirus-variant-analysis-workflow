# VEEV Variant Discovery Pipeline

Multi-sample viral variant calling pipeline using iVar and LoFreq with functional annotation.

## Quick Start

```bash
# 1. Validate inputs
./Scripts/check_inputs.sh

# 2. Run variant calling
./Scripts/run_ivar.sh      # iVar pipeline
./Scripts/run_lofreq.sh    # LoFreq pipeline

# 3. Annotate variants
./Scripts/annotate_all.sh
```

## Requirements

- Conda environments: `ivar_env`, `lofreq-env`, `annotation-env`
- Input BAMs in `Input/BAMs/`
- Reference FASTA in `Input/Reference/`
- Parameters in `config.sh`

## Features

- Viral read filtering (only viral contig)
- Multi-sample processing
- High-sensitivity variant calling (0.1% frequency, 5000× coverage)
- Automatic cleanup and logging

## Output

- `Ivar/` - iVar consensus sequences and variants  
- `LoFreq/` - LoFreq variant calls
- `Annotated_variants/` - Functional annotations

- bcftools (with csq plugin)
- Python 3
- Reference genome: Input/inh.fasta
- Annotation file: Input/INH.gff3

## Important Notes

- The `--local-csq` flag is CRITICAL for iVar processing - do not remove it
- iVar files must be converted from TSV to VCF format before annotation
- The pipeline expects specific directory structure as shown above