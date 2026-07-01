# TARGET Variant Discovery Pipeline

Multi-sample Target variant calling pipeline using iVar and LoFreq with functional annotation.

## Quick Start

```bash
# Stage-1 orchestrator (preferred)
source config.sh
./run_stage1.sh --mode bam --manifest work/stage1/manifest/sample_inputs.tsv

# FASTQ mode example
./run_stage1.sh --mode fastq --reference Input/Reference/inh.fasta \
  --threads 20 --target-contig target_contig

# Legacy compatibility path (kept for manual execution)
./Scripts/check_inputs.sh
./Scripts/run_ivar.sh
./Scripts/run_lofreq.sh
./Scripts/annotate_all.sh
```

## Requirements

- Conda environments: `ivar_env`, `lofreq-env`, `annotation-env`, `bwa_mem2_env`
- Input BAMs in `Input/BAMs/`
- Reference FASTA in `Input/Reference/` (default `Input/Reference/inh.fasta`)
- FASTQ mode inputs in `Input/FASTQ/` (optional)
- Parameters in `config.sh`

## Stage-1 Contract and Outputs

Preferred entrypoint is `run_stage1.sh`.

- `work/stage1/manifest/sample_inputs.tsv`
  - `sample_id`, `mode`, `input_bam`, `input_r1`, `input_r2`, `reference`, `TARGET_CONTIG`, `status`
- `work/stage1/run_status.tsv`
  - Per-sample, per-stage status rows:
    - `sample_id`, `stage`, `status`, `start_ts`, `end_ts`, `error`
- `work/stage1/stage1_report.tsv`
  - Run-level metrics and counters:
    - `started_at`, `mode`, `threads`, `manifest_ok`, `manifest_not_ok`, per-stage ok/skipped/failed counts, `overall_status`

Status meanings:
- `OK` -- stage completed for sample
- `SKIPPED` -- intentionally skipped (e.g. missing input, `no_variants`, existing BAM reuse, invalid manifest status)
- `FAILED` -- hard failure

## Output Layout (Preserved)

- `Ivar/` -- iVar consensus and variant outputs
- `LoFreq/` -- LoFreq variant outputs (including `variants.filtered.vcf.gz`)
- `Annotated_variants/` -- annotation outputs

`run_pipeline.sh` is a compatibility wrapper that delegates to `run_stage1.sh`.

## Important Notes

- `run_stage1.sh` performs validation via `check_inputs.sh`, writes manifest and status artifacts, and does not use workflow engines in Stage 1.
- The `--local-csq` flag is required for annotation.
- iVar TSV-to-VCF conversion happens in annotation stage.
- Empty variant outputs are supported and should not be treated as hard failures; annotation skips those samples (`no_variants`) by design.

## Optional mapping environment

```bash
conda env create -f envs/bwa_mem2_env.yml
conda activate bwa_mem2_env
```




