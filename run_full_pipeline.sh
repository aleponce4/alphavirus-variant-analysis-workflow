#!/bin/bash
# =====================================================================
# run_full_pipeline.sh — Run all phases sequentially (modular & dataset-aware)
#
# Allows running the core variant calling, extra selection analyses (SNPGenie),
# and/or haplotype reconstruction (CliqueSNV & VILOCA).
#
# Safe to run using existing variant results because variant calling scripts
# default to skipping completed samples (FORCE_RECALL=false).
#
# Usage (from WSL terminal):
#   # Run only selection analysis using existing variant results:
#   RUN_PHASE_EXTRACT=false RUN_PHASE_CALL=false RUN_PHASE_ANNOTATE=false \
#   RUN_PHASE_COVERAGE=false RUN_PHASE_SNPGENIE=true bash run_full_pipeline.sh
#
#   # Run normal variant calling + annotation (default behavior):
#   bash run_full_pipeline.sh
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export VARIANT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
export DATASET="${DATASET:-mouse_veev}"

# ── Phase Control Settings (Default: Run normal variant calling only) ────
export RUN_PHASE_EXTRACT="${RUN_PHASE_EXTRACT:-true}"
export RUN_PHASE_CALL="${RUN_PHASE_CALL:-true}"
export RUN_PHASE_ANNOTATE="${RUN_PHASE_ANNOTATE:-true}"
export RUN_PHASE_COVERAGE="${RUN_PHASE_COVERAGE:-true}"
export RUN_PHASE_SNPGENIE="${RUN_PHASE_SNPGENIE:-false}"
export RUN_PHASE_HAPLOTYPE="${RUN_PHASE_HAPLOTYPE:-false}"

# Global recall flag (default: false, skip completed steps/samples)
export FORCE_RECALL="${FORCE_RECALL:-false}"

# Master log file
MASTER_LOG="$VARIANT_ROOT/results/${DATASET}/full_pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$VARIANT_ROOT/results/${DATASET}"

# Print configuration summary
{
    echo "========================================================"
    echo "  VARIANT ANALYSIS PIPELINE: $DATASET"
    echo "  Started: $(date)"
    echo "  Log:     $MASTER_LOG"
    echo "========================================================"
    echo "  PHASE OPTIONS:"
    echo "    1. Extract Viral BAMs:  $RUN_PHASE_EXTRACT"
    echo "    2. LoFreq Calling:      $RUN_PHASE_CALL"
    echo "    3. Annotation:          $RUN_PHASE_ANNOTATE"
    echo "    4. Coverage QC:         $RUN_PHASE_COVERAGE"
    echo "    5. SNPGenie (extra):    $RUN_PHASE_SNPGENIE"
    echo "    6. Haplotypes (extra):   $RUN_PHASE_HAPLOTYPE"
    echo "  SETTINGS:"
    echo "    Force recall (re-run):  $FORCE_RECALL"
    echo "========================================================"
} | tee "$MASTER_LOG"

source ~/miniconda3/etc/profile.d/conda.sh

# ── Phase 1: Extract viral BAMs to SSD ──────────────────────
if [ "$RUN_PHASE_EXTRACT" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 1: Extracting viral BAMs to SSD" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    conda activate lofreq-env
    bash "$SCRIPT_DIR/Scripts/extract_viral_bams.sh" 2>&1 | tee -a "$MASTER_LOG"
    echo "    Phase 1 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 1 (Extract viral BAMs): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Phase 2: LoFreq variant calling ──────────────────────────
if [ "$RUN_PHASE_CALL" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 2: LoFreq variant calling (all samples)" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    conda activate lofreq-env
    bash "$SCRIPT_DIR/Scripts/run_lofreq.sh" 2>&1 | tee -a "$MASTER_LOG"
    echo "    Phase 2 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 2 (LoFreq calling): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Phase 3: Annotation ─────────────────────────────────────
if [ "$RUN_PHASE_ANNOTATE" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 3: Annotating variants with bcftools csq" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    conda activate annotation-env
    bash "$SCRIPT_DIR/Scripts/annotate_all.sh" 2>&1 | tee -a "$MASTER_LOG"
    echo "    Phase 3 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 3 (Annotation): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Phase 4: Coverage summary ──────────────────────────────
if [ "$RUN_PHASE_COVERAGE" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 4: Coverage summary" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    conda activate lofreq-env
    bash "$SCRIPT_DIR/Scripts/calculate_coverage.sh" 2>&1 | tee -a "$MASTER_LOG"
    echo "    Phase 4 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 4 (Coverage summary): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Phase 5: SNPGenie evolutionary selection (extra) ─────────
if [ "$RUN_PHASE_SNPGENIE" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 5: Running SNPGenie evolutionary selection analysis" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    bash "$SCRIPT_DIR/Scripts/run_snpgenie.sh" 2>&1 | tee -a "$MASTER_LOG"
    echo "    Phase 5 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 5 (SNPGenie): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Phase 6: Haplotype reconstruction (extra) ────────────────
if [ "$RUN_PHASE_HAPLOTYPE" = "true" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> PHASE 6: Reconstructing haplotypes (CliqueSNV & VILOCA)" | tee -a "$MASTER_LOG"
    echo "    $(date)" | tee -a "$MASTER_LOG"
    
    echo "  Running CliqueSNV..." | tee -a "$MASTER_LOG"
    bash "$SCRIPT_DIR/Scripts/run_cliquesnv.sh" 2>&1 | tee -a "$MASTER_LOG"
    
    echo "  Running VILOCA..." | tee -a "$MASTER_LOG"
    bash "$SCRIPT_DIR/Scripts/run_viloca.sh" 2>&1 | tee -a "$MASTER_LOG"
    
    echo "    Phase 6 complete: $(date)" | tee -a "$MASTER_LOG"
else
    echo ">>> Phase 6 (Haplotypes): SKIPPED" | tee -a "$MASTER_LOG"
fi

# ── Summary ────────────────────────────────────────────────
{
    echo ""
    echo "========================================================"
    echo "  PIPELINE RUN COMPLETE: $DATASET"
    echo "  Finished: $(date)"
    echo "========================================================"
} | tee -a "$MASTER_LOG"

# ── Sync results back to Windows workspace (lightweight files only) ────
source "$SCRIPT_DIR/config.sh"
if [ "$RESULTS_DIR" != "$VARIANT_ROOT/results/$DATASET" ]; then
    echo "" | tee -a "$MASTER_LOG"
    echo ">>> Generating run summary report..." | tee -a "$MASTER_LOG"
    (
        eval "$(conda shell.bash hook)"
        conda activate annotation-env
        python3 "$SCRIPT_DIR/Scripts/Helpers/generate_run_summary.py" \
            --results-dir "$RESULTS_DIR" \
            --manifest "$VARIANT_ROOT/config/samples_manifest.tsv" \
            --dataset "$DATASET"

        python3 "$SCRIPT_DIR/Scripts/Helpers/generate_variant_plots.py" \
            --results-dir "$RESULTS_DIR" \
            --manifest "$VARIANT_ROOT/config/samples_manifest.tsv" \
            --gff3 "${REFERENCE%.fasta}.gff3" \
            --dataset "$DATASET"

        python3 "$SCRIPT_DIR/Scripts/Helpers/generate_coverage_plots.py" \
            --results-dir "$RESULTS_DIR" \
            --manifest "$VARIANT_ROOT/config/samples_manifest.tsv" \
            --gff3 "${REFERENCE%.fasta}.gff3" \
            --dataset "$DATASET"

        python3 "$SCRIPT_DIR/Scripts/Helpers/generate_haplotype_plots.py" \
            --results-dir "$RESULTS_DIR" \
            --manifest "$VARIANT_ROOT/config/samples_manifest.tsv" \
            --gff3 "${REFERENCE%.fasta}.gff3" \
            --dataset "$DATASET"

        python3 "$SCRIPT_DIR/Scripts/export_excel_variants.py" \
            --results-dir "$RESULTS_DIR"
    ) 2>&1 | tee -a "$MASTER_LOG"

    echo "" | tee -a "$MASTER_LOG"
    echo ">>> Syncing lightweight results to Windows workspace..." | tee -a "$MASTER_LOG"
    echo "    Source: $RESULTS_DIR" | tee -a "$MASTER_LOG"
    echo "    Target: $VARIANT_ROOT/results/$DATASET" | tee -a "$MASTER_LOG"
    
    WINDOWS_RESULTS="$VARIANT_ROOT/results/$DATASET"
    mkdir -p "$WINDOWS_RESULTS/Annotated_variants" "$WINDOWS_RESULTS/Coverage" "$WINDOWS_RESULTS/SNPGenie" "$WINDOWS_RESULTS/LoFreq" "$WINDOWS_RESULTS/tables" "$WINDOWS_RESULTS/Plots"
    
    # Copy consolidated summaries
    cp -f "$RESULTS_DIR/consolidated_sample_summary.tsv" "$WINDOWS_RESULTS/" 2>/dev/null || true
    cp -f "$RESULTS_DIR/run_summary_report.md" "$WINDOWS_RESULTS/" 2>/dev/null || true

    # Copy tables and Plots
    cp -rf "$RESULTS_DIR/tables/"* "$WINDOWS_RESULTS/tables/" 2>/dev/null || true
    cp -rf "$RESULTS_DIR/Plots/"* "$WINDOWS_RESULTS/Plots/" 2>/dev/null || true

    # Copy Annotated variants, Coverage summaries, and SNPGenie selection analyses
    cp -rf "$RESULTS_DIR/Annotated_variants/"* "$WINDOWS_RESULTS/Annotated_variants/" 2>/dev/null || true
    cp -rf "$RESULTS_DIR/Coverage/"* "$WINDOWS_RESULTS/Coverage/" 2>/dev/null || true
    cp -rf "$RESULTS_DIR/SNPGenie/"* "$WINDOWS_RESULTS/SNPGenie/" 2>/dev/null || true
    
    # Copy LoFreq VCFs and QC stats (skipping the heavy BAM files)
    for sample_dir in "$RESULTS_DIR/LoFreq"/*/; do
        if [ -d "$sample_dir" ]; then
            sample=$(basename "$sample_dir")
            mkdir -p "$WINDOWS_RESULTS/LoFreq/$sample"
            cp -f "$sample_dir/variants.filtered.vcf.gz"* "$WINDOWS_RESULTS/LoFreq/$sample/" 2>/dev/null || true
            cp -f "$sample_dir/qc_stats.txt" "$WINDOWS_RESULTS/LoFreq/$sample/" 2>/dev/null || true
        fi
    done
    
    # Sync Haplotypes if executed
    if [ -d "$RESULTS_DIR/CliqueSNV" ]; then
        mkdir -p "$WINDOWS_RESULTS/CliqueSNV"
        cp -rf "$RESULTS_DIR/CliqueSNV/Analysis" "$WINDOWS_RESULTS/CliqueSNV/" 2>/dev/null || true
        for sample_dir in "$RESULTS_DIR/CliqueSNV"/*/; do
            if [ -d "$sample_dir" ] && [ "$(basename "$sample_dir")" != "logs" ] && [ "$(basename "$sample_dir")" != "Analysis" ]; then
                sample=$(basename "$sample_dir")
                mkdir -p "$WINDOWS_RESULTS/CliqueSNV/$sample"
                cp -rf "$sample_dir/tf_"* "$WINDOWS_RESULTS/CliqueSNV/$sample/" 2>/dev/null || true
            fi
        done
    fi
    
    if [ -d "$RESULTS_DIR/VILOCA" ]; then
        mkdir -p "$WINDOWS_RESULTS/VILOCA"
        for sample_dir in "$RESULTS_DIR/VILOCA"/*/; do
            if [ -d "$sample_dir" ] && [ "$(basename "$sample_dir")" != "logs" ]; then
                sample=$(basename "$sample_dir")
                mkdir -p "$WINDOWS_RESULTS/VILOCA/$sample"
                cp -f "$sample_dir/cooccurring_mutations.csv" "$WINDOWS_RESULTS/VILOCA/$sample/" 2>/dev/null || true
                cp -f "$sample_dir/coverage.txt" "$WINDOWS_RESULTS/VILOCA/$sample/" 2>/dev/null || true
            fi
        done
    fi
    
    # Copy master pipeline log
    cp -f "$MASTER_LOG" "$WINDOWS_RESULTS/"
    echo "    Sync complete! Final reports are now visible in Windows at:" | tee -a "$MASTER_LOG"
    echo "    $WINDOWS_RESULTS" | tee -a "$MASTER_LOG"
    echo "========================================================" | tee -a "$MASTER_LOG"
fi
