#!/usr/bin/env python3
"""
generate_variant_plots.py — Publication-grade iSNV density & frequency spectrum visualization.
Parses machine VCF outputs, calculates percentage metrics (0-100%), writes summary TSV,
and generates variant density and allele frequency spectrum figures.
"""

import sys
import os
import argparse
import gzip
import math
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

__version__ = "1.3.0"

def parse_vcf(vcf_path: Path) -> list[dict]:
    """Parse a single VCF file (plain text or gzip) and extract variant records."""
    records = []
    sample_name = vcf_path.name.replace(".csq.vcf", "").replace(".variants.filtered.vcf.gz", "").replace(".vcf.gz", "").replace(".vcf", "")
    
    open_fn = gzip.open if str(vcf_path).endswith('.gz') else open
    
    try:
        with open_fn(vcf_path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                info_str = parts[7]
                
                dp = 0
                af = 0.0
                for item in info_str.split(';'):
                    if item.startswith('DP=') and not item.startswith('DP4='):
                        try:
                            dp = int(item.split('=')[1])
                        except ValueError:
                            dp = 0
                    elif item.startswith('AF='):
                        try:
                            af = float(item.split('=')[1].split(',')[0])
                        except ValueError:
                            af = 0.0
                
                af_pct = af * 100.0
                if af_pct >= 1.0:
                    tier = ">=1%"
                elif af_pct >= 0.1:
                    tier = "0.1-1% (candidate)"
                else:
                    tier = "<0.1% (exploratory)"
                
                records.append({
                    'sample': sample_name,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'depth': dp,
                    'af_proportion': af,
                    'af_percent': af_pct,
                    'tier': tier
                })
    except Exception as e:
        print(f"Warning: Failed to parse {vcf_path}: {e}", file=sys.stderr)
        
    return records

def generate_plots(df: pd.DataFrame, outdir: Path, dataset: str):
    """Generate publication-grade PNG plots for variant density and frequency spectrum."""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    # 1. Variant Density Plot (Positional Distribution)
    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    
    if df.empty:
        ax.text(0.5, 0.5, "No variants called in dataset", horizontalalignment='center', verticalalignment='center', fontsize=14)
        ax.set_axis_off()
    else:
        palette = {">=1%": "#2b5c8f", "0.1-1% (candidate)": "#e06d53", "<0.1% (exploratory)": "#8c8c8c"}
        
        sns.scatterplot(
            data=df,
            x="pos",
            y="af_percent",
            hue="tier",
            style="sample",
            palette=palette,
            s=70,
            alpha=0.8,
            ax=ax
        )
        
        ax.set_title(f"Intra-host Variant Frequency Distribution ({dataset})", fontsize=14, fontweight="bold", pad=12)
        ax.set_xlabel("Genomic Position (bp)", fontsize=12)
        ax.set_ylabel("Allele Frequency (%)", fontsize=12)
        ax.set_ylim(-1, max(105, df['af_percent'].max() + 5))
        
        # Reference thresholds
        ax.axhline(1.0, color="#2b5c8f", linestyle="--", linewidth=1, alpha=0.7, label="1.0% Threshold")
        ax.axhline(0.1, color="#e06d53", linestyle=":", linewidth=1, alpha=0.7, label="0.1% Noise Floor")
        
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)
    
    plt.tight_layout()
    density_plot_path = outdir / f"{dataset}_variant_density.png"
    fig.savefig(density_plot_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Variant density plot written to {density_plot_path}")
    
    # 2. Allele Frequency Spectrum Histogram
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    if df.empty:
        ax.text(0.5, 0.5, "No variants called in dataset", horizontalalignment='center', verticalalignment='center', fontsize=14)
        ax.set_axis_off()
    else:
        sns.histplot(
            data=df,
            x="af_percent",
            bins=30,
            hue="tier",
            palette=palette,
            multiple="stack",
            edgecolor="black",
            linewidth=0.5,
            ax=ax
        )
        
        ax.set_title(f"Allele Frequency Spectrum ({dataset})", fontsize=14, fontweight="bold", pad=12)
        ax.set_xlabel("Allele Frequency (%)", fontsize=12)
        ax.set_ylabel("Variant Count", fontsize=12)
    
    plt.tight_layout()
    spectrum_plot_path = outdir / f"{dataset}_allele_frequency_spectrum.png"
    fig.savefig(spectrum_plot_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Allele frequency spectrum plot written to {spectrum_plot_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate publication-grade variant density and frequency plots.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--vcf-dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--outdir", required=True, help="Output directory for plots and summary TSV")
    parser.add_argument("--dataset", default="viral_analysis", help="Dataset name label")

    args = parser.parse_args()

    vcf_dir = Path(args.vcf_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    vcf_files = sorted(list(vcf_dir.glob("*.vcf")) + list(vcf_dir.glob("*.vcf.gz")))
    
    all_records = []
    for vcf_file in vcf_files:
        all_records.extend(parse_vcf(vcf_file))
        
    df = pd.DataFrame(all_records)
    
    # Save summary TSV with explicit % units
    summary_tsv_path = outdir / "variant_frequency_summary_pct.tsv"
    if not df.empty:
        df.to_csv(summary_tsv_path, sep="\t", index=False)
    else:
        with summary_tsv_path.open("w", encoding="utf-8") as f:
            f.write("sample\tchrom\tpos\tref\talt\tdepth\taf_proportion\taf_percent\ttier\n")
            
    print(f"Processed {len(df)} variants across {len(vcf_files)} VCFs. Output TSV: {summary_tsv_path}")

    # Generate plots
    generate_plots(df, outdir, args.dataset)

if __name__ == "__main__":
    main()
