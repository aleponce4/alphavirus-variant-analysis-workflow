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

### GFF3 Annotation Format

`bcftools csq` requires an Ensembl-style GFF3 file with the following features:

- `gene` features with `ID=gene:NAME` and `biotype=protein_coding`
- `mRNA`/`transcript` features with `ID=transcript:NAME`, `Parent=gene:NAME`, and `biotype=protein_coding`
- `CDS` features with `Parent=transcript:NAME` (note the `transcript:` prefix)
- `exon` features with `Parent=transcript:NAME`

Example:
```
##gff-version 3
KP282671.1\tGenBank\tgene\t45\t7526\t.\t+\t.\tID=gene:nsP1234;Name=nsP1234;biotype=protein_coding
KP282671.1\tGenBank\tmRNA\t45\t7526\t.\t+\t.\tID=transcript:nsP1234;Parent=gene:nsP1234;Name=nsP1234;biotype=protein_coding
KP282671.1\tGenBank\tCDS\t45\t7526\t.\t+\t0\tID=cds:nsP1234;Parent=transcript:nsP1234;Name=nsP1234;biotype=protein_coding
KP282671.1\tGenBank\texon\t45\t7526\t.\t+\t.\tID=exon:nsP1234;Parent=transcript:nsP1234
```

GFF3 files can be generated from NCBI GenBank using Biopython. See the
`convert_genbank_to_gff3` approach in the variant_analysis workspace.

### RNA-seq Data Considerations

When using this pipeline with RNA-seq data (e.g., STAR-aligned BAMs):

- Set `LOFREQ_BAQ=0` in `config.sh` — `lofreq viterbi` crashes on spliced
  alignments (CIGAR `N` operations from STAR)
- Set `LOFREQ_ENABLE_INDELQUAL=0` in `config.sh` — `lofreq indelqual` also
  crashes on CIGAR `N` operations
- No primer trimming is needed (RNA-seq, not amplicon-based)
- Coverage may be uneven; tune `LOFREQ_MIN_VARIANT_DEPTH` and
  `IVAR_MIN_VARIANT_DEPTH` based on `calculate_coverage.sh` results