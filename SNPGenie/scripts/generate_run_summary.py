#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate coverage, variant calling, and selection statistics into a consolidated report."
    )
    parser.add_argument("--results-dir", required=True, help="Path to Windows or WSL results folder for the dataset")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--dataset", required=True, help="Active dataset name (e.g. mouse_veev)")
    return parser.parse_args()

def load_manifest(manifest_path: Path, dataset_filter: str) -> dict[str, dict]:
    """Loads metadata for samples matching the dataset."""
    samples = {}
    with manifest_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("dataset") == dataset_filter:
                bam_name = row.get("bam_name")
                if bam_name:
                    samples[bam_name] = {
                        "sample_id": row.get("sample_id"),
                        "dpi": row.get("dpi"),
                        "treatment": row.get("treatment"),
                        "study_code": row.get("study_code"),
                        "replicate_group": row.get("replicate_group"),
                    }
    return samples

def load_coverage(coverage_summary_path: Path) -> dict[str, dict]:
    """Loads coverage statistics per sample."""
    coverage = {}
    if coverage_summary_path.exists():
        with coverage_summary_path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample = row.get("sample")
                if sample:
                    coverage[sample] = {
                        "mean_depth": row.get("mean_depth", "0"),
                        "percent_above_1000x": row.get("percent_above_1000x", "0"),
                    }
    return coverage

def load_lofreq_qc(lofreq_dir: Path) -> dict[str, int]:
    """Loads called variant counts per sample from LoFreq qc_stats.txt."""
    counts = {}
    if lofreq_dir.exists():
        for sample_dir in lofreq_dir.glob("*/"):
            if not sample_dir.is_dir():
                continue
            qc_file = sample_dir / "qc_stats.txt"
            if qc_file.exists():
                sample = sample_dir.name
                counts[sample] = 0
                with qc_file.open("r", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("Variants (filtered):"):
                            try:
                                counts[sample] = int(line.split(":")[1].strip())
                            except ValueError:
                                pass
    return counts

def load_ivar_qc(ivar_dir: Path) -> dict[str, int]:
    """Loads called variant counts per sample from iVar qc_stats.txt."""
    counts = {}
    if ivar_dir.exists():
        for sample_dir in ivar_dir.glob("*/"):
            if not sample_dir.is_dir():
                continue
            qc_file = sample_dir / "qc_stats.txt"
            if qc_file.exists():
                sample = sample_dir.name
                counts[sample] = 0
                with qc_file.open("r", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("Variants:"):
                            try:
                                counts[sample] = int(line.split(":")[1].strip())
                            except ValueError:
                                pass
    return counts

def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    manifest_path = Path(args.manifest)
    dataset = args.dataset

    # Input paths
    coverage_summary_file = results_dir / "Coverage" / "coverage_summary.tsv"
    lofreq_dir = results_dir / "LoFreq"
    ivar_dir = results_dir / "Ivar"
    selection_key_file = results_dir / "SNPGenie" / "analysis" / "summary" / "selection_gene_key_table.tsv"
    methods_file = results_dir / "SNPGenie" / "analysis" / "summary" / "methods_key_parameters.tsv"

    # Outputs
    consolidated_tsv = results_dir / "consolidated_sample_summary.tsv"
    consolidated_md = results_dir / "run_summary_report.md"

    # Load data
    samples_meta = load_manifest(manifest_path, dataset)
    coverage_data = load_coverage(coverage_summary_file)
    lofreq_variants = load_lofreq_qc(lofreq_dir)
    ivar_variants = load_ivar_qc(ivar_dir)

    if not samples_meta:
        print(f"WARNING: No samples found in manifest for dataset '{dataset}'")
        return

    # Aggregate consolidated rows
    consolidated_rows = []
    for sample, meta in sorted(samples_meta.items()):
        cov = coverage_data.get(sample, {"mean_depth": "0", "percent_above_1000x": "0"})
        lf_vars = lofreq_variants.get(sample, 0)
        iv_vars = ivar_variants.get(sample, 0)
        
        # Check if sample was skipped (no virus)
        no_virus = "No"
        if lf_vars == 0 and float(cov["mean_depth"]) < 1:
            no_virus = "Yes"

        consolidated_rows.append({
            "sample_name": sample,
            "sample_id": meta["sample_id"],
            "dpi": meta["dpi"],
            "treatment": meta["treatment"],
            "replicate_group": meta["replicate_group"],
            "mean_coverage": cov["mean_depth"],
            "percent_above_1000x": cov["percent_above_1000x"],
            "lofreq_variants": str(lf_vars),
            "ivar_variants": str(iv_vars),
            "no_viral_reads": no_virus
        })

    # Write consolidated TSV
    fieldnames = [
        "sample_name", "sample_id", "dpi", "treatment", "replicate_group", 
        "mean_coverage", "percent_above_1000x", "lofreq_variants", "ivar_variants", "no_viral_reads"
    ]
    with consolidated_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(consolidated_rows)
    print(f"Wrote consolidated sample spreadsheet to: {consolidated_tsv}")

    # Generate Markdown Report
    with consolidated_md.open("w", encoding="utf-8") as f:
        f.write(f"# Variant Analysis Run Summary Report: `{dataset}`\n\n")
        
        f.write("## 1. Run Metadata\n")
        f.write(f"* **Dataset**: {dataset}\n")
        f.write(f"* **Samples manifest**: `{manifest_path.name}`\n")
        f.write(f"* **Total samples listed**: {len(samples_meta)}\n")
        f.write(f"* **Samples with viral reads**: {sum(1 for r in consolidated_rows if r['no_viral_reads'] == 'No')}\n")
        f.write(f"* **Samples skipped (Mock / no virus)**: {sum(1 for r in consolidated_rows if r['no_viral_reads'] == 'Yes')}\n\n")

        # Evolutionary selection findings
        f.write("## 2. Evolutionary Selection Findings (SNPGenie)\n")
        if selection_key_file.exists():
            f.write("The table below details evolutionary selection pressure ($\pi_N - \pi_S$) calculated for each of the 9 mature viral proteins:\n\n")
            with selection_key_file.open("r", encoding="utf-8") as sf:
                reader = csv.reader(sf, delimiter="\t")
                header = next(reader)
                
                # Format headers for markdown
                f.write("| " + " | ".join(header) + " |\n")
                f.write("| " + " | ".join(["---"] * len(header)) + " |\n")
                for row in reader:
                    f.write("| " + " | ".join(row) + " |\n")
            f.write("\n*Note: Limma FDR indicates Benjamini-Hochberg corrected p-values for differential selection across timepoints.*\n\n")
        else:
            f.write("*No selection analysis results found. Run the SNPGenie phase to generate selection findings.*\n\n")

        # Methods and key parameters
        f.write("## 3. Key Parameters Used\n")
        if methods_file.exists():
            f.write("| Parameter | Setting |\n| --- | --- |\n")
            with methods_file.open("r", encoding="utf-8") as mf:
                reader = csv.reader(mf, delimiter="\t")
                next(reader) # skip header
                for row in reader:
                    if len(row) >= 2:
                        f.write(f"| {row[0]} | {row[1]} |\n")
            f.write("\n")
        else:
            f.write("*Methods parameter file not found.*\n\n")

        # Sample summary table
        f.write("## 4. Consolidated Sample Performance\n")
        f.write("The table below consolidates sequencing coverage and variant calling performance for each sample:\n\n")
        f.write("| Sample | ID | DPI | Treatment | Mean Coverage | % covered $\ge 1000\\times$ | LoFreq Variants | iVar Variants | Skipped (No Virus) |\n")
        f.write("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for r in consolidated_rows:
            f.write(f"| {r['sample_name']} | {r['sample_id']} | {r['dpi']} | {r['treatment']} | {r['mean_coverage']}x | {r['percent_above_1000x']}% | {r['lofreq_variants']} | {r['ivar_variants']} | {r['no_viral_reads']} |\n")

    print(f"Wrote summary Markdown report to: {consolidated_md}")

if __name__ == "__main__":
    main()
