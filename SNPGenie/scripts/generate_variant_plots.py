#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import entropy

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract variants from VCFs, calculate Shannon Entropy, and generate publication-ready variant distribution and diversity plots."
    )
    parser.add_argument("--results-dir", required=True, help="Path to WSL or Windows results folder for the dataset")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--gff3", required=True, help="Path to the viral_only.gff3 reference file")
    parser.add_argument("--dataset", required=True, help="Active dataset name (e.g. mouse_veev)")
    parser.add_argument("--min-coverage", type=int, default=1000, help="Coverage depth filtering threshold (default: 1000)")
    parser.add_argument("--min-frequency", type=float, default=0.01, help="Allele frequency filtering threshold (default: 0.01)")
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

def parse_bcsq(bcsq_field: str) -> dict:
    """Parses BCSQ consequence field: consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change"""
    default = {
        'consequence': None, 'gene': None, 'transcript': None,
        'biotype': None, 'strand': None, 'amino_acid_change': None, 'dna_change': None
    }
    if not bcsq_field or bcsq_field == "." or str(bcsq_field).startswith("@"):
        return default
        
    bcsq_str = bcsq_field
    if isinstance(bcsq_field, tuple):
        if len(bcsq_field) > 0:
            bcsq_str = bcsq_field[0]
        else:
            return default
            
    parts = str(bcsq_str).split('|')
    return {
        'consequence': parts[0] if len(parts) > 0 else None,
        'gene': parts[1] if len(parts) > 1 else None,
        'transcript': parts[2] if len(parts) > 2 else None,
        'biotype': parts[3] if len(parts) > 3 else None,
        'strand': parts[4] if len(parts) > 4 else None,
        'amino_acid_change': parts[5] if len(parts) > 5 else None,
        'dna_change': parts[6] if len(parts) > 6 else None
    }

def parse_vcf_file(vcf_path: Path) -> pd.DataFrame:
    """Parses a VCF file manually to extract variants and annotations."""
    variants = []
    if not vcf_path.exists():
        return pd.DataFrame()
        
    with vcf_path.open('r', encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 8:
                variant_data = {
                    'CHROM': parts[0],
                    'POS': int(parts[1]),
                    'REF': parts[3],
                    'ALT': parts[4],
                    'QUAL': parts[5] if parts[5] != '.' else None,
                    'FILTER': parts[6] if parts[6] != '.' else 'PASS',
                    'DP': None,
                    'AF': None,
                    'BCSQ_raw': None
                }
                
                # Parse INFO field
                info_parts = parts[7].split(';')
                for info in info_parts:
                    if '=' in info:
                        key, value = info.split('=', 1)
                        if key == 'DP':
                            try:
                                variant_data['DP'] = int(value)
                            except ValueError:
                                pass
                        elif key == 'AF':
                            try:
                                variant_data['AF'] = float(value)
                            except ValueError:
                                pass
                        elif key == 'BCSQ':
                            variant_data['BCSQ_raw'] = value
                
                # Parse FORMAT field for iVar alternate frequency fallback
                if len(parts) >= 10:
                    format_keys = parts[8].split(':')
                    sample_values = parts[9].split(':')
                    if 'ALT_FREQ' in format_keys:
                        idx = format_keys.index('ALT_FREQ')
                        if len(sample_values) > idx:
                            try:
                                variant_data['AF'] = float(sample_values[idx])
                            except ValueError:
                                pass
                    if 'DP' in format_keys and variant_data['DP'] is None:
                        idx = format_keys.index('DP')
                        if len(sample_values) > idx:
                            try:
                                variant_data['DP'] = int(sample_values[idx])
                            except ValueError:
                                pass

                # Parse BCSQ consequence
                bcsq_parsed = parse_bcsq(variant_data['BCSQ_raw'])
                variant_data.update(bcsq_parsed)
                variants.append(variant_data)
                
    return pd.DataFrame(variants) if variants else pd.DataFrame()

def af_vector_from_group(group: pd.DataFrame) -> np.ndarray:
    """Builds a normalized allele-frequency vector representing reference and alternative alleles."""
    afs = []
    for _, row in group.iterrows():
        af = row.get("AF")
        if pd.notna(af):
            afs.append(float(af))
            
    if not afs:
        return np.nan
        
    max_alt_total = sum(afs)
    p_ref = max(0.0, 1.0 - max_alt_total)
    freqs = np.array([p_ref] + afs)
    
    # Normalize vector
    s = freqs.sum()
    return freqs / s if s > 0 else np.nan

def entropy_from_group_AF(group: pd.DataFrame) -> float:
    """Calculates Shannon Entropy from an allele frequency vector."""
    freqs = af_vector_from_group(group)
    if isinstance(freqs, float) and np.isnan(freqs):
        return np.nan
    freqs = freqs[freqs > 0]
    return float(entropy(freqs, base=np.e))

def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    manifest_path = Path(args.manifest)
    gff3_path = Path(args.gff3)
    dataset = args.dataset
    
    # 1. Load coordinates and manifest
    gene_positions = load_gene_coordinates(gff3_path)
    samples_meta = load_manifest(manifest_path, dataset)
    
    if not samples_meta:
        print(f"ERROR: No samples found in manifest for dataset '{dataset}'")
        return
        
    print(f"Loaded {len(gene_positions)} gene coordinates and {len(samples_meta)} samples.")
    
    # 2. Extract variants from annotated directories
    annotated_dir = results_dir / "Annotated_variants"
    lofreq_dir = annotated_dir / "LoFreq"
    ivar_dir = annotated_dir / "Ivar"
    
    lofreq_list = []
    ivar_list = []
    
    for sample, meta in samples_meta.items():
        # LoFreq VCF
        lf_vcf = lofreq_dir / f"{sample}_filtered.vcf"
        if not lf_vcf.exists():
            lf_vcf = lofreq_dir / f"{sample}.vcf"
        if lf_vcf.exists():
            df = parse_vcf_file(lf_vcf)
            if not df.empty:
                df["sample_name"] = sample
                df["dpi"] = meta["dpi"]
                df["replicate"] = meta["replicate"]
                df["treatment"] = meta["treatment"]
                df["caller"] = "LoFreq"
                lofreq_list.append(df)
                
        # iVar VCF
        iv_vcf = ivar_dir / f"{sample}_filtered.vcf"
        if not iv_vcf.exists():
            iv_vcf = ivar_dir / f"{sample}.vcf"
        if iv_vcf.exists():
            df = parse_vcf_file(iv_vcf)
            if not df.empty:
                df["sample_name"] = sample
                df["dpi"] = meta["dpi"]
                df["replicate"] = meta["replicate"]
                df["treatment"] = meta["treatment"]
                df["caller"] = "iVar"
                ivar_list.append(df)
                
    # Combine caller variants
    lf_all = pd.concat(lofreq_list, ignore_index=True) if lofreq_list else pd.DataFrame()
    iv_all = pd.concat(ivar_list, ignore_index=True) if ivar_list else pd.DataFrame()
    
    # Create tables and plots directories
    tables_dir = results_dir / "tables"
    plots_dir = results_dir / "Plots"
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 3. Filter variants
    min_dp = args.min_coverage
    min_af = args.min_frequency
    
    print(f"Filtering variants: Coverage >= {min_dp}x, Frequency >= {min_af * 100}%")
    
    lf_filt = pd.DataFrame()
    iv_filt = pd.DataFrame()
    
    if not lf_all.empty:
        # Fill missing DP and AF defaults
        lf_all["DP"] = lf_all["DP"].fillna(0)
        lf_all["AF"] = lf_all["AF"].fillna(0.0)
        lf_filt = lf_all[(lf_all["DP"] >= min_dp) & (lf_all["AF"] >= min_af)]
        
    if not iv_all.empty:
        iv_all["DP"] = iv_all["DP"].fillna(0)
        iv_all["AF"] = iv_all["AF"].fillna(0.0)
        iv_filt = iv_all[(iv_all["DP"] >= min_dp) & (iv_all["AF"] >= min_af) & (iv_all["FILTER"] == "PASS")]

    combined_variants = pd.concat([lf_filt, iv_filt], ignore_index=True)
    if combined_variants.empty:
        print("No variants left after filtering. Skipping downstream plots.")
        return
        
    # Save combined variants table
    combined_variants.to_csv(tables_dir / "combined_filtered_variants.tsv", sep="\t", index=False)
    print(f"Wrote filtered variants to: {tables_dir / 'combined_filtered_variants.tsv'}")
    
    # 4. Generate Shannon Entropy per site
    print("Calculating Shannon Entropy (population diversity) per site...")
    # Use all detected variants down to 0.1% (0.001) frequency to capture true population diversity
    lf_entropy_source = lf_all[(lf_all["DP"] >= min_dp) & (lf_all["AF"] >= 0.001)] if not lf_all.empty else pd.DataFrame()
    iv_entropy_source = iv_all[(iv_all["DP"] >= min_dp) & (iv_all["AF"] >= 0.001)] if not iv_all.empty else pd.DataFrame()
    entropy_source = lf_entropy_source if not lf_entropy_source.empty else iv_entropy_source
    
    group_keys = ["sample_name", "dpi", "POS"]
    entropy_per_site = (
        entropy_source
        .groupby(group_keys, dropna=False)
        .apply(entropy_from_group_AF)
        .reset_index(name="Entropy")
    )
    
    entropy_per_site["Entropy"] = entropy_per_site["Entropy"].fillna(0.0)
    
    # Aggregate by DPI and Position
    agg_entropy = (
        entropy_per_site.groupby(['dpi', 'POS'])['Entropy']
        .agg(['mean', 'sem'])
        .reset_index()
    )
    
    # Export full genome position grid for DPIs present
    unique_dpis = sorted(list(entropy_per_site["dpi"].unique()))
    genome_end = int(max(v['end'] for v in gene_positions.values())) if gene_positions else int(entropy_source["POS"].max())
    
    position_grid = pd.MultiIndex.from_product(
        [unique_dpis, np.arange(1, genome_end + 1)],
        names=['DPI', 'Position']
    ).to_frame(index=False)
    
    entropy_genome = (
        position_grid
        .merge(agg_entropy.rename(columns={'dpi': 'DPI', 'POS': 'Position', 'mean': 'Shannon_Entropy_Mean', 'sem': 'Shannon_Entropy_SEM'}), 
               on=['DPI', 'Position'], how='left')
        .sort_values(['DPI', 'Position'])
    )
    
    entropy_genome['Shannon_Entropy_Mean'] = entropy_genome['Shannon_Entropy_Mean'].fillna(0.0)
    entropy_genome['Shannon_Entropy_SEM'] = entropy_genome['Shannon_Entropy_SEM'].fillna(0.0)
    entropy_genome.to_csv(tables_dir / "shannon_entropy_by_position.tsv", sep="\t", index=False)
    print(f"Wrote genome-wide entropy values to: {tables_dir / 'shannon_entropy_by_position.tsv'}")
    
    # 5. Visualizations
    print("Generating publication plots...")
    
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
    
    # Colors for mutation types
    color_map = {
        'missense': OKABE_ITO["vermilion"],   # Nonsynonymous
        'synonymous': OKABE_ITO["green"],    # Synonymous
        'other': OKABE_ITO["grey"]           # Non-coding or other
    }
    
    # Plot 1: Variant Distribution Scatter Plot
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=300)
    plot_data = lf_filt.copy() if not lf_filt.empty else iv_filt.copy()
    
    def get_mut_type(c):
        c = str(c).lower()
        if 'missense' in c:
            return 'missense'
        elif 'synonymous' in c:
            return 'synonymous'
        return 'other'
        
    plot_data['mutation_type'] = plot_data['consequence'].apply(get_mut_type)
    
    for m_type in ['synonymous', 'missense', 'other']:
        subset = plot_data[plot_data['mutation_type'] == m_type]
        if not subset.empty:
            label_name = 'Nonsynonymous' if m_type == 'missense' else ('Synonymous' if m_type == 'synonymous' else 'Other')
            ax.scatter(subset['POS'], subset['AF'] * 100,
                       c=color_map[m_type],
                       label=f"{label_name} ({len(subset)})",
                       s=12, alpha=0.85, edgecolors='black', linewidths=0.3)
                       
    # Overlay gene boundaries
    if gene_positions:
        for gene, info in gene_positions.items():
            ax.axvline(info['start'], color='gray', linewidth=0.8, alpha=0.4, zorder=1)
            ax.axvline(info['end'], color='gray', linewidth=0.8, alpha=0.4, zorder=1)
            
            mid = (info['start'] + info['end']) / 2
            ax.text(mid, 85, gene, ha='center', va='top', fontsize=8, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85, edgecolor='none'))
                    
    ax.set_xlabel('Nucleotide Position', fontsize=8, fontweight='bold')
    ax.set_ylabel('Variant Frequency (%)', fontsize=8, fontweight='bold')
    ax.set_title(f'Variant Frequency Distribution ({dataset})', fontsize=9, fontweight='bold')
    ax.set_xlim(0, genome_end + 100)
    ax.set_yscale('log')
    ax.set_ylim(0.8, 105)
    
    import matplotlib.ticker as ticker
    ax.set_yticks([1, 10, 100])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(ticker.NullFormatter())
    ax.legend(fontsize=8, frameon=True, fancybox=False, edgecolor='black', loc='lower left')
    
    apply_plot_theme(ax)
    
    plt.tight_layout()
    plt.savefig(svg_dir / "variant_frequency_distribution.svg", format='svg', bbox_inches='tight')
    plt.savefig(plots_dir / "variant_frequency_distribution.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Shannon Entropy Plot
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=300)
    
    # Map DPI values to Okabe-Ito colors (Orange, Green, Blue, Vermilion)
    dpi_colors = {
        1: OKABE_ITO["orange"],
        2: OKABE_ITO["green"],
        3: OKABE_ITO["blue"],
        4: OKABE_ITO["vermilion"]
    }
    
    for idx, dpi_val in enumerate(unique_dpis):
        subset = agg_entropy[agg_entropy['dpi'] == dpi_val]
        if subset.empty:
            continue
            
        color = dpi_colors.get(dpi_val, list(OKABE_ITO.values())[idx % len(OKABE_ITO)])
        ax.plot(subset['POS'], subset['mean'], label=f'DPI {dpi_val}',
                color=color, lw=1.2, alpha=0.9, zorder=3)
        ax.fill_between(subset['POS'],
                        subset['mean'] - subset['sem'],
                        subset['mean'] + subset['sem'],
                        color=color, alpha=0.18, lw=0, zorder=2)
                        
    # Overlay gene boundaries
    if gene_positions:
        for gene, info in gene_positions.items():
            ax.axvline(info['start'], color='gray', linewidth=0.5, alpha=0.3, zorder=1)
            mid = (info['start'] + info['end']) / 2
            ax.text(mid, 0.65, gene, ha='center', va='top', fontsize=8, fontweight='bold', zorder=4,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.8, edgecolor='none'))
                    
    ax.set_xlabel('Nucleotide Position', fontsize=8, fontweight='bold')
    ax.set_ylabel('Shannon Entropy (nats)', fontsize=8, fontweight='bold')
    ax.set_title(f'Viral Population Diversity (Shannon Entropy) ({dataset})', fontsize=9, fontweight='bold')
    ax.set_xlim(0, genome_end + 100)
    ax.set_ylim(0, 0.7)
    ax.legend(fontsize=8, loc='upper left')
    
    apply_plot_theme(ax)
    
    plt.tight_layout()
    plt.savefig(svg_dir / "shannon_entropy_distribution.svg", format='svg', bbox_inches='tight')
    plt.savefig(plots_dir / "shannon_entropy_distribution.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: 2-way Venn Diagram (if package is available and both callers have data)
    if not lf_filt.empty and not iv_filt.empty:
        try:
            from matplotlib_venn import venn2
            lf_muts = set(lf_filt['POS'].astype(str) + '_' + lf_filt['REF'] + '>' + lf_filt['ALT'])
            iv_muts = set(iv_filt['POS'].astype(str) + '_' + iv_filt['REF'] + '>' + iv_filt['ALT'])
            
            fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
            venn2([lf_muts, iv_muts], set_labels=('LoFreq', 'iVar'), set_colors=(OKABE_ITO["blue"], OKABE_ITO["orange"]))
            plt.title('Shared Mutations Between LoFreq and iVar', fontsize=12, fontweight='bold')
            plt.tight_layout()
            plt.savefig(plots_dir / "shared_mutations_venn.png", format='png', dpi=300, bbox_inches='tight')
            plt.close()
            print("Generated shared mutations Venn diagram.")
        except ImportError:
            print("matplotlib_venn not installed. Skipping Venn diagram plot.")
            
    print(f"Plots successfully saved to: {plots_dir}")

if __name__ == "__main__":
    main()
