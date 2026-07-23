#!/usr/bin/env python3
import argparse
import csv
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract coverage per position across replicates and generate publication-ready average coverage plots per DPI."
    )
    parser.add_argument("--results-dir", required=True, help="Path to WSL or Windows results folder for the dataset")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--gff3", required=True, help="Path to the viral_only.gff3 reference file")
    parser.add_argument("--dataset", required=True, help="Active dataset name (e.g. mouse_veev)")
    return parser.parse_args()

def load_gene_coordinates(gff3_path: Path) -> dict[str, dict]:
    """Dynamically parses mature peptide boundaries from the active GFF3 reference."""
    genes = {}
    if not gff3_path.exists():
        print(f"WARNING: GFF3 reference not found at {gff3_path}. Using fallback coordinates.")
        return {}
    
    with gff3_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 9 and parts[2] == "CDS":
                attrs = parts[8].strip()
                name_match = re.search(r"Name=([^;]+)", attrs)
                if name_match:
                    name = name_match.group(1)
                    genes[name] = {
                        "start": int(parts[3]),
                        "end": int(parts[4])
                    }
    return genes

def load_manifest(manifest_path: Path, dataset_filter: str) -> dict[str, dict]:
    """Loads experimental metadata (DPI, replicate) for active dataset, dynamically calculating replicates."""
    samples = {}
    group_counts = {}
    with manifest_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("dataset") == dataset_filter:
                bam_name = row.get("bam_name")
                if bam_name:
                    rep_group = row.get("replicate_group", "default")
                    group_counts[rep_group] = group_counts.get(rep_group, 0) + 1
                    rep_num = group_counts[rep_group]
                    
                    samples[bam_name] = {
                        "dpi": int(row.get("dpi", 0)),
                        "replicate": rep_num,
                        "treatment": row.get("treatment", "unknown")
                    }
    return samples

def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    manifest_path = Path(args.manifest)
    gff3_path = Path(args.gff3)
    dataset = args.dataset
    
    # 1. Load metadata
    gene_positions = load_gene_coordinates(gff3_path)
    samples_meta = load_manifest(manifest_path, dataset)
    
    if not samples_meta:
        print(f"ERROR: No samples found in manifest for dataset '{dataset}'")
        return
        
    print(f"Loaded {len(gene_positions)} gene coordinates and {len(samples_meta)} samples.")
    
    # 2. Read per-position coverage txt files
    coverage_dir = results_dir / "Coverage"
    coverage_list = []
    
    for sample, meta in sorted(samples_meta.items()):
        cov_file = coverage_dir / f"{sample}_coverage.txt"
        if cov_file.exists():
            # Read samtools depth format: CHROM  POS  DEPTH
            try:
                df = pd.read_csv(cov_file, sep="\t", header=None, names=["CHROM", "Position", "Depth"])
                if not df.empty:
                    df["sample_name"] = sample
                    df["dpi"] = meta["dpi"]
                    df["replicate"] = meta["replicate"]
                    df["treatment"] = meta["treatment"]
                    coverage_list.append(df)
            except Exception as e:
                print(f"  ⚠ Error reading {cov_file.name}: {e}")
                
    if not coverage_list:
        print("No coverage files found. Skipping coverage plot.")
        return
        
    combined_coverage = pd.concat(coverage_list, ignore_index=True)
    
    # Only average replicates for Challenge/Infected samples (exclude Mock where depth is ~0)
    # E.g. treatment = "Challenge"
    challenge_coverage = combined_coverage[combined_coverage["treatment"] == "Challenge"]
    if challenge_coverage.empty:
        print("No Challenge/Infected sample coverage found. Falling back to all samples.")
        challenge_coverage = combined_coverage
        
    # Group by DPI and Position to calculate Mean and SEM across replicates
    print("Calculating mean coverage and SEM across replicates...")
    agg_coverage = (
        challenge_coverage.groupby(['dpi', 'Position'])['Depth']
        .agg(['mean', 'sem'])
        .reset_index()
    )
    
    # Create tables and plots directories
    tables_dir = results_dir / "tables"
    plots_dir = results_dir / "Plots"
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Export aggregated coverage table
    agg_coverage.to_csv(tables_dir / "genome_coverage_by_position.tsv", sep="\t", index=False)
    print(f"Wrote coverage stats to: {tables_dir / 'genome_coverage_by_position.tsv'}")
    
    # 3. Visualization
    print("Generating publication coverage plots...")
    
    # Create SVG subfolder
    svg_dir = plots_dir / "SVG"
    svg_dir.mkdir(parents=True, exist_ok=True)
    
    # Okabe-Ito color palette definitions
    OKABE_ITO = {
        "orange": "#E69F00",
        "skyblue": "#56B4E9",
        "green": "#009E73",
        "yellow": "#F0E442",
        "blue": "#0072B2",
        "vermilion": "#D55E00",
        "purple": "#CC79A7",
        "grey": "#999999"
    }
    
    # Setup styles
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 8
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.rcParams['savefig.dpi'] = 300
    
    def apply_plot_theme(ax):
        ax.grid(False)
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_color('black')
            ax.spines[spine].set_linewidth(0.8)
        ax.tick_params(axis='both', which='both', top=False, right=False, width=0.8, labelsize=8)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_facecolor('white')
    
    unique_dpis = sorted(list(agg_coverage["dpi"].unique()))
    genome_end = int(agg_coverage["Position"].max())
    
    def format_lab_title(ds_name: str, metric: str) -> str:
        ds = ds_name.lower()
        if "veev" in ds:
            study_str = "VEEV Studies 045 & 047"
        elif "eeev" in ds:
            study_str = "EEEV Studies 046 & 048"
        else:
            study_str = ds.upper()
        return f"{metric} — {study_str}"

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=300)
    
    # Map DPI values to Okabe-Ito colors (Orange, Green, Blue, Vermilion)
    dpi_colors = {
        1: OKABE_ITO["orange"],
        2: OKABE_ITO["green"],
        3: OKABE_ITO["blue"],
        4: OKABE_ITO["vermilion"]
    }
    
    max_depth_plotted = 0.0
    
    for idx, dpi_val in enumerate(unique_dpis):
        subset = agg_coverage[agg_coverage['dpi'] == dpi_val]
        if subset.empty:
            continue
            
        color = dpi_colors.get(dpi_val, list(OKABE_ITO.values())[idx % len(OKABE_ITO)])
        ax.plot(subset['Position'], subset['mean'].clip(lower=1.0), label=f'DPI {dpi_val}',
                color=color, lw=1.2, alpha=0.9, zorder=3)
        # Handle nan SEMs (when only 1 replicate exists)
        sem_vals = subset['sem'].fillna(0.0)
        ax.fill_between(subset['Position'],
                        (subset['mean'] - sem_vals).clip(lower=1.0),
                        (subset['mean'] + sem_vals).clip(lower=1.0),
                        color=color, alpha=0.18, lw=0, zorder=2)
                        
        max_depth_plotted = max(max_depth_plotted, (subset['mean'] + sem_vals).max())
        
    # Overlay gene boundaries in top header space above maximum depth
    if gene_positions:
        safe_max = max(1.0, max_depth_plotted)
        ax.axhline(safe_max, color='gray', linestyle='--', linewidth=0.8, zorder=2)
        stagger_map = {'E3': safe_max * 1.15, '6K': safe_max * 1.15, 'TF': safe_max * 1.15}
        
        for gene, info in gene_positions.items():
            ax.axvline(info['start'], color='gray', linewidth=0.5, linestyle=':', alpha=0.5, zorder=1)
            ax.axvline(info['end'], color='gray', linewidth=0.5, linestyle=':', alpha=0.5, zorder=1)
            
            mid = (info['start'] + info['end']) / 2
            y_pos = stagger_map.get(gene, safe_max * 1.05)
            
            ax.text(mid, y_pos, gene, ha='center', va='center', fontsize=7, fontweight='bold', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.12', facecolor='white', alpha=0.9, edgecolor='#94A3B8', linewidth=0.5))
                    
    ax.set_xlabel('Nucleotide Position', fontsize=8, fontweight='bold')
    ax.set_ylabel('Sequencing Depth', fontsize=8, fontweight='bold')
    ax.set_title(format_lab_title(dataset, "Viral Genome Coverage Depth"), fontsize=9, fontweight='bold')
    ax.set_xlim(0, genome_end + 100)
    ax.set_yscale('log')
    ax.set_ylim(0.8, max_depth_plotted * 1.25)
    
    import matplotlib.ticker as ticker
    # Setup custom log ticks for depth (e.g. 1, 10, 100, 1000, 10000, 100000)
    # We choose formatter that handles numbers nicely
    ax.get_yaxis().set_major_locator(ticker.LogLocator(base=10.0))
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y)}" if y >= 1 else f"{y}"))
    ax.get_yaxis().set_minor_formatter(ticker.NullFormatter())
    ax.legend(fontsize=8, loc='lower right', frameon=True, edgecolor='black', fancybox=False)
    
    apply_plot_theme(ax)
    
    plt.tight_layout()
    plt.savefig(svg_dir / "genome_coverage_distribution.svg", format='svg', bbox_inches='tight')
    plt.savefig(plots_dir / "genome_coverage_distribution.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Coverage plots successfully saved to: {plots_dir}")

if __name__ == "__main__":
    main()
