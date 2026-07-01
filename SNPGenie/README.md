# SNPGenie analysis for TARGET INH-9813 dataset

This folder contains a reproducible SNPGenie workflow for 12 LoFreq-called samples
from `Variant_discovery_pipeline/LoFreq`.

## Inputs used

- Reference FASTA: `Variant_discovery_pipeline/Input/Reference/target_only.fasta`
- Annotation GFF3: `Variant_discovery_pipeline/Input/Reference/target_reference.gff3`
- Variant source: `Variant_discovery_pipeline/LoFreq/*/variants.filtered.vcf.gz`

The GFF3 is converted to a SNPGenie-compatible GTF in `input/reference/`.

## Folder layout

- `repo/`: cloned SNPGenie repository
- `scripts/`: helper scripts for prep, run, and summary
- `input/`: staged FASTA/GFF3/GTF and per-sample VCFs
- `runs/`: temporary per-sample working directories
- `output/`: final SNPGenie outputs by threshold and sample
- `logs/`: per-sample SNPGenie logs and run status TSV files
- `summary/`: merged tables across all samples
- `analysis/delta_selection/`: downstream delta/Kruskal and limma statistics
- `analysis/summary/`: compact downstream reporting tables (clean final deliverables)

## Run

From workspace root (`044_NH`):

```bash
bash SNPGenie/scripts/prepare_inputs.sh
bash SNPGenie/scripts/run_snpgenie_batch.sh
python3 SNPGenie/scripts/summarize_snpgenie_outputs.py
python3 SNPGenie/scripts/analyze_delta_selection.py
Rscript SNPGenie/scripts/analyze_delta_limma.R
python3 SNPGenie/scripts/build_compact_selection_tables.py \
  --base SNPGenie/analysis/delta_selection/detailed \
  --report-dir SNPGenie/analysis/summary \
  --product-summary SNPGenie/summary/product_results_all_samples.tsv \
  --threshold minfreq_0p01
```

To generate compact outputs in a chosen report folder (and choose threshold for the
limma key table):

```bash
python3 SNPGenie/scripts/build_compact_selection_tables.py \
  --base SNPGenie/analysis/delta_selection/detailed \
  --report-dir SNPGenie/analysis/summary \
  --product-summary SNPGenie/summary/product_results_all_samples.tsv \
  --threshold minfreq_0p01
```

If needed, install dependency for statistics:

```bash
pip install scipy
```

## Thresholds implemented

- `minfreq_0p001` (0.1%)
- `minfreq_0p01` (1%)

## Notes

- Workflow uses `--vcfformat=2` for LoFreq pooled-depth VCFs (`INFO` includes `DP` and `AF`).
- One sample VCF can be sparse/empty at high stringency; this pipeline still stages it from
  filtered output and logs success/failure explicitly.
- This workflow runs SNPGenie in within-pool mode (`snpgenie.pl`) for each sample independently.
- `analyze_delta_selection.py` outputs:
  - `analysis/delta_selection/delta_per_sample.tsv` (piN, piS, Δ by sample/gene)
  - `analysis/delta_selection/delta_kruskal_by_gene.tsv` (Kruskal–Wallis by gene across dpi + BH-FDR)
  - add `--write-detailed` to also emit per-threshold files
- `analyze_delta_limma.R` outputs limma empirical Bayes model results for `Δ ~ dpi`:
  - overall moderated F-test by gene
  - planned contrasts (`dpi3-dpi1`, `dpi5-dpi3`, `dpi5-dpi1`) with moderated t statistics
  - BH-FDR adjusted p-values across genes
  - add `--write-detailed` to also emit per-threshold/per-contrast files
  - `build_compact_selection_tables.py` writes compact and key reporting tables:
    - `analysis/summary/selection_gene_summary.tsv`
    - `analysis/summary/selection_gene_contrasts.tsv`
    - `analysis/summary/selection_gene_key_table.tsv`
      - one row per gene
      - mean `piN`, `piS`, and `Δ = piN - piS` by dpi
      - limma BH-FDR (overall and `dpi5_vs_dpi1`)
      - dpi5 direction label based on mean Δ sign
    - `analysis/summary/methods_key_parameters.tsv`
      - ready-to-use key parameters for methods text (threshold, sample counts, model/FDR settings)

 Detailed and intermediate outputs can be kept in:
 - `analysis/delta_selection/detailed/`

