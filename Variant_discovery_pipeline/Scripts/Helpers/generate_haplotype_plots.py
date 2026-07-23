#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import sem

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate publication-ready plots for CliqueSNV, VILOCA, and SNPGenie selection analyses."
    )
    parser.add_argument("--results-dir", required=True, help="Path to WSL or Windows results folder for the dataset")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--gff3", required=True, help="Path to the viral_only.gff3 reference file")
    parser.add_argument("--dataset", required=True, help="Active dataset name (e.g. mouse_veev)")
    parser.add_argument("--no-normalize", action="store_true", help="Do not normalize CliqueSNV haplotype frequencies to sum to 1.0")
    return parser.parse_args()

def load_gene_coordinates(gff3_path: Path) -> dict[str, dict]:
    """Parses CDS boundaries from the active GFF3 reference."""
    genes = {}
    if not gff3_path.exists():
        print(f"WARNING: GFF3 reference not found at {gff3_path}.")
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
    """Loads experimental metadata (DPI, treatment) for active dataset."""
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
                        "treatment": row.get("treatment", "unknown"),
                        "route": row.get("route", "unknown"),
                        "study_code": row.get("study_code", "")
                    }
    return samples

# =====================================================================
# matlotlib theme and colors
# =====================================================================
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

# =====================================================================
# 1. CliqueSNV stacked bar plot
# =====================================================================
def load_cliquesnv_data(cliquesnv_dir: Path, manifest_samples: dict) -> pd.DataFrame:
    rows = []
    for sample_id, meta in manifest_samples.items():
        tf_dir = cliquesnv_dir / sample_id / "tf_0p01"
        json_path = tf_dir / "primary_only.json"
        
        if not json_path.exists() and tf_dir.exists():
            jsons = list(tf_dir.glob("*.json"))
            if jsons:
                json_path = jsons[0]
                
        if json_path.exists():
            try:
                with open(json_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                haplotypes = data.get("haplotypes", [])
                for idx, hap in enumerate(haplotypes):
                    freq = float(hap.get("frequency", 0.0))
                    snps_field = hap.get("snps", "")
                    if not snps_field or snps_field == "[]" or snps_field == "-":
                        snp_positions = tuple()
                    else:
                        snp_positions = tuple(sorted([int(x) for x in re.findall(r'\d+', snps_field)]))
                    
                    rows.append({
                        "sample_id": sample_id,
                        "dpi": meta["dpi"],
                        "treatment": meta["treatment"],
                        "route": meta.get("route", "unknown"),
                        "study_code": meta.get("study_code", ""),
                        "frequency": freq,
                        "snp_positions": snp_positions
                    })
            except Exception as e:
                print(f"  Warning: failed to load {json_path}: {e}")
    return pd.DataFrame(rows)

def cluster_haplotypes(df: pd.DataFrame) -> dict[tuple[int, ...], str]:
    if df.empty:
        return {}
    
    unique_sigs = list(set(df["snp_positions"].unique()))
    sig_to_cluster = {}
    
    # Always keep Ref separate
    if tuple() in unique_sigs:
        sig_to_cluster[tuple()] = "Ref (Wild-type)"
        other_sigs = [s for s in unique_sigs if s != tuple()]
    else:
        other_sigs = list(unique_sigs)
        
    if not other_sigs:
        return sig_to_cluster
        
    if len(other_sigs) <= 5:
        for idx, sig in enumerate(other_sigs):
            sig_to_cluster[sig] = f"Haplotype_{idx+1}"
        return sig_to_cluster
        
    # Compute symmetric difference distance matrix
    n = len(other_sigs)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            s1 = set(other_sigs[i])
            s2 = set(other_sigs[j])
            dist = len(s1 ^ s2)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
            
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    
    try:
        condensed_dist = squareform(dist_matrix)
        Z = linkage(condensed_dist, method='complete')
        labels = fcluster(Z, t=5, criterion='maxclust')
        for idx, sig in enumerate(other_sigs):
            cluster_num = labels[idx]
            sig_to_cluster[sig] = f"Haplotype_Cluster_{cluster_num}"
    except Exception as e:
        print(f"  Warning: cliquesnv signature clustering failed, using fallback labeling: {e}")
        for idx, sig in enumerate(other_sigs):
            sig_to_cluster[sig] = f"Haplotype_{idx+1}"
            
    return sig_to_cluster

def draw_cliquesnv_plot(df: pd.DataFrame, plots_dir: Path, svg_dir: Path, dataset: str, no_normalize: bool = False):
    if df.empty:
        print("  Skipping CliqueSNV plot (no haplotype data loaded).")
        return
        
    # Remove DPI 1 samples as requested
    df = df[df["dpi"] > 1].copy()
    if df.empty:
        print("  Skipping CliqueSNV plot (no samples left after DPI > 1 filter).")
        return
        
    sig_to_cluster = cluster_haplotypes(df)
    df["Haplotype"] = [sig_to_cluster.get(x, "Unknown") for x in df["snp_positions"]]
    
    # Aggregate duplicate groups if any
    pivot_df = df.groupby(["sample_id", "dpi", "treatment", "route", "study_code", "Haplotype"])["frequency"].sum().unstack(fill_value=0.0)
    
    if not no_normalize:
        row_sums = pivot_df.sum(axis=1)
        pivot_df = pivot_df.div(row_sums.where(row_sums > 0, 1.0), axis=0)
        
    pivot_df = pivot_df.reset_index()

    def get_route_rank(r_str):
        r = str(r_str).lower()
        if "intra" in r or r == "in":
            return 0
        return 1

    pivot_df["route_rank"] = pivot_df["route"].apply(get_route_rank)
    pivot_df = pivot_df.sort_values(by=["route_rank", "dpi", "treatment", "sample_id"]).reset_index(drop=True)
    
    # Sort Haplotypes so Ref is first, then rest
    all_hap_cols = [c for c in pivot_df.columns if c not in ["sample_id", "dpi", "treatment", "route", "study_code", "route_rank"]]
    ref_col = [c for c in all_hap_cols if "ref" in c.lower()]
    other_cols = sorted([c for c in all_hap_cols if c not in ref_col])
    ordered_haps = ref_col + other_cols
    
    # Color scheme mapping
    palette = [OKABE_ITO["skyblue"], OKABE_ITO["orange"], OKABE_ITO["green"], 
               OKABE_ITO["blue"], OKABE_ITO["vermilion"], OKABE_ITO["purple"], 
               OKABE_ITO["yellow"], OKABE_ITO["grey"]]
    color_map = {hap: palette[idx % len(palette)] for idx, hap in enumerate(ordered_haps)}
    
    def format_lab_title(ds_name: str, metric: str) -> str:
        ds = ds_name.lower()
        if "veev" in ds:
            study_str = "VEEV Studies 045 & 047"
        elif "eeev" in ds:
            study_str = "EEEV Studies 046 & 048"
        else:
            study_str = ds.upper()
        return f"{metric} — {study_str}"

    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=300)
    x_pos = np.arange(len(pivot_df))
    bottom = np.zeros(len(pivot_df))
    
    for hap in ordered_haps:
        heights = pivot_df[hap].values
        ax.bar(x_pos, heights, bottom=bottom, label=hap, 
               color=color_map[hap], width=0.7, edgecolor='black', linewidth=0.3)
        bottom += heights
        
    if not no_normalize:
        ax.set_ylabel("Haplotype Frequency (reconstructed >= 1%)", fontsize=8, fontweight='bold')
    else:
        ax.set_ylabel("Haplotype Frequency", fontsize=8, fontweight='bold')
    ax.set_xlabel("Sample ID", fontsize=8, fontweight='bold')
    ax.set_title(format_lab_title(dataset, "CliqueSNV Haplotype Composition"), fontsize=9, fontweight='bold')
    ax.set_ylim(0, 1.35)
    
    # Clean X-axis tick labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels([s.replace("s", "") for s in pivot_df["sample_id"]], rotation=90, fontsize=6)
    
    # Draw top brackets for Route / Study & sub-brackets for DPI
    routes = pivot_df.groupby("route_rank", sort=False)
    for r_rank, r_group in routes:
        r_start = r_group.index.min()
        r_end = r_group.index.max()
        r_mid = (r_start + r_end) / 2
        
        raw_route = r_group["route"].iloc[0]
        study_id = r_group["study_code"].iloc[0]
        r_name_clean = "Intranasal" if r_rank == 0 else "Subcutaneous"
        route_label = f"{r_name_clean} (Study {study_id})" if study_id else r_name_clean
        
        # Top Bracket for Route
        y_r = 1.20
        ax.plot([r_start - 0.35, r_end + 0.35], [y_r, y_r], color='black', lw=1.0)
        ax.plot([r_start - 0.35, r_start - 0.35], [y_r - 0.02, y_r], color='black', lw=1.0)
        ax.plot([r_end + 0.35, r_end + 0.35], [y_r - 0.02, y_r], color='black', lw=1.0)
        ax.text(r_mid, y_r + 0.02, route_label, ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # Group by DPI within route
        dpis = r_group.groupby("dpi", sort=False)
        for d_val, d_group in dpis:
            d_start = d_group.index.min()
            d_end = d_group.index.max()
            d_mid = (d_start + d_end) / 2
            
            # Sub-bracket for DPI
            y_d = 1.07
            ax.plot([d_start - 0.35, d_end + 0.35], [y_d, y_d], color='gray', lw=0.8)
            ax.plot([d_start - 0.35, d_start - 0.35], [y_d - 0.015, y_d], color='gray', lw=0.8)
            ax.plot([d_end + 0.35, d_end + 0.35], [y_d - 0.015, y_d], color='gray', lw=0.8)
            ax.text(d_mid, y_d + 0.015, f"DPI {d_val}", ha='center', va='bottom', fontsize=7, fontweight='bold')
                
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, edgecolor='black', fancybox=False, fontsize=8)
    apply_plot_theme(ax)
    
    plt.tight_layout()
    plt.savefig(svg_dir / "cliquesnv_haplotype_composition.svg", format='svg', bbox_inches='tight')
    plt.savefig(plots_dir / "cliquesnv_haplotype_composition.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved CliqueSNV Stacked Bar Plot")

# =====================================================================
# 2. VILOCA Track + Linkage Arcs Plot
# =====================================================================
def parse_vcf_variants(vcf_path: Path) -> pd.DataFrame:
    rows = []
    if not vcf_path.exists():
        return pd.DataFrame(rows)
    with vcf_path.open("r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 8:
                pos = int(parts[1])
                info = parts[7]
                af_match = re.search(r"\b(?:AF|Freq)=([0-9.]+)", info, re.IGNORECASE)
                if af_match:
                    af = float(af_match.group(1))
                    rows.append({"POS": pos, "AF": af})
    return pd.DataFrame(rows)

def draw_viloca_arc(ax, x1, x2, height, alpha=0.5, color='gray'):
    t = np.linspace(0, np.pi, 50)
    x = (x1 + x2)/2 + (x2 - x1)/2 * np.cos(t)
    y = height * np.sin(t)
    ax.plot(x, y, color=color, alpha=alpha, lw=0.8)

def generate_viloca_plots(viloca_dir: Path, plots_dir: Path, svg_dir: Path, gene_coords: dict, dataset: str, manifest_samples: dict):
    # Find samples with cooccurring_mutations.csv files using manifest
    samples_found = []
    for sample_id in manifest_samples.keys():
        sample_dir = viloca_dir / sample_id
        vcf = sample_dir / "snv.vcf"
        cooc = sample_dir / "cooccurring_mutations.csv"
        if vcf.exists() and cooc.exists() and os.path.getsize(cooc) > 10:
            samples_found.append(sample_id)
                
    if not samples_found:
        print("  Skipping VILOCA plots (no completed samples with co-occurring mutation arcs found).")
        return
        
    # Select up to 3 representative samples
    samples_to_plot = sorted(samples_found)[:3]
    print(f"  Generating VILOCA plots for representative samples: {', '.join(samples_to_plot)}")
    
    # Identify gene boundaries color mapping
    gene_colors = [OKABE_ITO["orange"], OKABE_ITO["green"], OKABE_ITO["blue"], 
                   OKABE_ITO["vermilion"], OKABE_ITO["purple"], OKABE_ITO["yellow"]]
    gene_color_map = {gene: gene_colors[idx % len(gene_colors)] for idx, gene in enumerate(gene_coords.keys())}
    
    for sample in samples_to_plot:
        vcf_path = viloca_dir / sample / "snv.vcf"
        cooc_path = viloca_dir / sample / "cooccurring_mutations.csv"
        
        snvs = parse_vcf_variants(vcf_path)
        try:
            links = pd.read_csv(cooc_path)
        except Exception:
            continue
            
        if snvs.empty or links.empty:
            continue
            
        # Parse links position columns dynamically
        pos_cols = [c for c in links.columns if "pos" in c.lower() or "position" in c.lower()]
        freq_cols = [c for c in links.columns if "freq" in c.lower() or "prob" in c.lower()]
        
        if len(pos_cols) < 2 or not freq_cols:
            continue
            
        pos1_col, pos2_col = pos_cols[0], pos_cols[1]
        freq_col = freq_cols[0]
        
        # Filter links
        links[freq_col] = pd.to_numeric(links[freq_col], errors="coerce").fillna(0.0)
        filtered_links = links[(links[freq_col] >= 0.05) & (abs(links[pos1_col] - links[pos2_col]) < 2000)]
        
        # Build 2-panel figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 5), sharex=True, dpi=300)
        
        # Panel 1: Variant frequencies
        # Map SNV positions to genes for color coding
        def get_gene_color(pos):
            for gene, info in gene_coords.items():
                if info["start"] <= pos <= info["end"]:
                    return gene_color_map[gene]
            return OKABE_ITO["grey"]
            
        snvs["color"] = snvs["POS"].apply(get_gene_color)
        
        for color, group in snvs.groupby("color"):
            ax1.scatter(group["POS"], group["AF"] * 100, c=color, s=15, edgecolors='black', linewidths=0.3, alpha=0.85)
            
        ax1.set_ylabel("Minor Allele Freq (%)", fontsize=8, fontweight='bold')
        ax1.set_title(f"VILOCA Local Haplotype SNVs & Linkage: {sample} ({dataset})", fontsize=9, fontweight='bold')
        ax1.set_yscale('log')
        ax1.set_ylim(0.4, 105)
        apply_plot_theme(ax1)
        
        # Panel 2: Linkage arcs
        max_dist = 0
        for _, row in filtered_links.iterrows():
            x1, x2 = int(row[pos1_col]), int(row[pos2_col])
            w = float(row[freq_col])
            dist = abs(x2 - x1)
            max_dist = max(max_dist, dist)
            draw_viloca_arc(ax2, x1, x2, height=dist/2, alpha=min(1.0, w + 0.1), color=OKABE_ITO["vermilion"])
            
        ax2.set_ylabel("Linkage Distance (bp)", fontsize=8, fontweight='bold')
        ax2.set_xlabel("Genome Position (nt)", fontsize=8, fontweight='bold')
        ax2.set_ylim(0, max(500, max_dist / 2 * 1.1))
        apply_plot_theme(ax2)
        
        # Draw GFF3 gene annotations on the top plot
        if gene_coords:
            for gene, info in gene_coords.items():
                ax1.axvline(info['start'], color='gray', linewidth=0.5, alpha=0.3, zorder=1)
                ax1.axvline(info['end'], color='gray', linewidth=0.5, alpha=0.3, zorder=1)
                mid = (info['start'] + info['end']) / 2
                ax1.text(mid, 85, gene, ha='center', va='top', fontsize=6, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.8, edgecolor='none'))
                         
        plt.tight_layout()
        plt.savefig(svg_dir / f"viloca_linkage_{sample}.svg", format='svg', bbox_inches='tight')
        plt.savefig(plots_dir / f"viloca_linkage_{sample}.png", format='png', dpi=300, bbox_inches='tight')
        plt.close()
        
    print("  ✓ Saved VILOCA Genome Linkage Track Plots")

# =====================================================================
# 3. SNPGenie sliding window selection plot
# =====================================================================
def run_sliding_window(df: pd.DataFrame, window_size=40, step=1) -> pd.DataFrame:
    """Calculates sliding-window piN/piS in Python."""
    results = []
    
    # Sort by site (genomic coordinate), which is numeric
    df = df.sort_values("site")
    
    n_codons = len(df)
    if n_codons < window_size:
        return pd.DataFrame()
        
    for start_idx in range(0, n_codons - window_size + 1, step):
        window_df = df.iloc[start_idx : start_idx + window_size]
        
        mid_site = int(window_df["site"].mean())
        
        # Calculate raw piN and piS
        sum_n_diffs = window_df["N_diffs"].sum()
        sum_n_sites = window_df["N_sites"].sum()
        sum_s_diffs = window_df["S_diffs"].sum()
        sum_s_sites = window_df["S_sites"].sum()
        
        pi_n = sum_n_diffs / sum_n_sites if sum_n_sites > 0 else 0.0
        pi_s = sum_s_diffs / sum_s_sites if sum_s_sites > 0 else 0.0
        
        # Calculate piN - piS selection difference metric
        pi_diff = pi_n - pi_s
        
        results.append({
            "site": mid_site,
            "pi_n": pi_n,
            "pi_s": pi_s,
            "ratio": pi_diff
        })
        
    return pd.DataFrame(results)

def generate_snpgenie_plots(snpgenie_dir: Path, plots_dir: Path, svg_dir: Path, manifest_samples: dict, dataset: str):
    # Group samples by DPI and treatment for plotting
    samples_found = []
    for sample, meta in manifest_samples.items():
        codon_file = snpgenie_dir / "output" / "minfreq_0p01" / sample / "codon_results.txt"
        if codon_file.exists() and os.path.getsize(codon_file) > 100:
            samples_found.append((sample, meta["dpi"], meta["treatment"]))
            
    if not samples_found:
        print("  Skipping SNPGenie plots (no codon results found).")
        return
        
    print(f"  Calculating SNPGenie sliding windows for {len(samples_found)} samples...")
    
    # Process each sample
    all_window_data = []
    for sample, dpi, treatment in samples_found:
        codon_file = snpgenie_dir / "output" / "minfreq_0p01" / sample / "codon_results.txt"
        try:
            df = pd.read_csv(codon_file, sep="\t")
        except Exception:
            continue
            
        # Standard columns required
        required = ["codon", "site", "product", "N_diffs", "N_sites", "S_diffs", "S_sites"]
        if any(col not in df.columns for col in required):
            continue
            
        # Filter out overlap regions or double count products, separate by non-structural and structural
        df["region"] = df["product"].apply(lambda p: "Non-structural (nsP1-4)" if str(p).startswith("nsP") else "Structural")
        
        for region, group in df.groupby("region"):
            # run sliding windows (size=40 codons, step=2)
            win_df = run_sliding_window(group, window_size=40, step=2)
            if not win_df.empty:
                win_df["sample_id"] = sample
                win_df["dpi"] = dpi
                win_df["treatment"] = treatment
                win_df["region"] = region
                all_window_data.append(win_df)
                
    if not all_window_data:
        print("  Skipping SNPGenie plots (could not compile window data).")
        return
        
    master_win_df = pd.concat(all_window_data, ignore_index=True)
    
    # Aggregate by DPI, Region, and Position across replicates (mean +/- SEM)
    agg_df = (
        master_win_df.groupby(["region", "dpi", "site"])
        .agg(
            mean_ratio=("ratio", "mean"),
            mean_se=("ratio", "sem")
        )
        .reset_index()
    )
    
    agg_df["mean_se"] = agg_df["mean_se"].fillna(0.0)
    
    # Draw selection lines faceted by Region
    regions = [r for r in ["Non-structural (nsP1-4)", "Structural"] if r in agg_df["region"].unique()]
    if not regions:
        print("  Skipping SNPGenie plots (no matching regions found).")
        return
        
    unique_dpis = sorted(list(agg_df["dpi"].unique()))
    dpi_colors = {
        1: OKABE_ITO["orange"],
        2: OKABE_ITO["green"],
        3: OKABE_ITO["blue"],
        4: OKABE_ITO["vermilion"]
    }
    
    def format_lab_title(ds_name: str, metric: str) -> str:
        ds = ds_name.lower()
        if "veev" in ds:
            study_str = "VEEV Studies 045 & 047"
        elif "eeev" in ds:
            study_str = "EEEV Studies 046 & 048"
        else:
            study_str = ds.upper()
        return f"{metric} — {study_str}"

    n_regions = len(regions)
    fig, axes = plt.subplots(n_regions, 1, figsize=(7, 2.5 * n_regions), squeeze=False, sharex=False, dpi=300)
    
    for idx, region in enumerate(regions):
        ax = axes[idx, 0]
        region_subset = agg_df[agg_df["region"] == region]
        
        # Add neutrality line
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.8, label="Neutrality (piN = piS)")
        
        for d_idx, dpi_val in enumerate(unique_dpis):
            dpi_subset = region_subset[region_subset["dpi"] == dpi_val].sort_values("site")
            if dpi_subset.empty:
                continue
                
            color = dpi_colors.get(dpi_val, list(OKABE_ITO.values())[d_idx % len(OKABE_ITO)])
            
            # Smooth window lines using a rolling average to reduce noisy spikes
            dpi_subset["smoothed_ratio"] = dpi_subset["mean_ratio"].rolling(window=5, center=True, min_periods=1).mean()
            dpi_subset["smoothed_se"] = dpi_subset["mean_se"].rolling(window=5, center=True, min_periods=1).mean()
            
            ax.plot(dpi_subset["site"], dpi_subset["smoothed_ratio"], label=f"DPI {dpi_val}", 
                    color=color, lw=1.2, alpha=0.9)
            
            # Fill standard error ribbon
            ax.fill_between(
                dpi_subset["site"],
                (dpi_subset["smoothed_ratio"] - dpi_subset["smoothed_se"]),
                (dpi_subset["smoothed_ratio"] + dpi_subset["smoothed_se"]),
                color=color, alpha=0.15
            )
            
        ax.set_title(f"{region} ({format_lab_title(dataset, 'Selection Pressure')})", fontsize=8, fontweight="bold")
        ax.set_ylabel(r"Selection ($\pi_N - \pi_S$)", fontsize=8, fontweight="bold")
        ax.set_xlabel("Genome Position (nt)", fontsize=8, fontweight="bold")
        apply_plot_theme(ax)
        if idx == 0:
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=True, edgecolor='black', fancybox=False, fontsize=7)
            
            # Add explanatory box below legend explaining positive, neutral, and purifying selection
            explanation = (
                r"$\mathbf{Selection\ Metric\ (\pi_N - \pi_S):}$" "\n"
                r"• $\mathbf{> 0}$: Positive Selection" "\n"
                r"   (adaptive / diversifying)" "\n"
                r"• $\mathbf{= 0}$: Neutral Evolution" "\n"
                r"• $\mathbf{< 0}$: Purifying Selection" "\n"
                r"   (functional constraint)"
            )
            ax.text(
                1.02, 0.42, explanation, transform=ax.transAxes, fontsize=6.5,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='#F8F9FA', edgecolor='black', lw=0.6)
            )
            
    plt.tight_layout()
    plt.savefig(svg_dir / "snpgenie_selection_sliding_window.svg", format='svg', bbox_inches='tight')
    plt.savefig(plots_dir / "snpgenie_selection_sliding_window.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved SNPGenie Selection Pressure sliding-window plot")

# =====================================================================
# Main execution
# =====================================================================
def main():
    args = parse_args()
    results_dir = Path(args.results_dir)
    manifest_path = Path(args.manifest)
    gff3_path = Path(args.gff3)
    dataset = args.dataset
    
    plots_dir = results_dir / "Plots"
    svg_dir = plots_dir / "SVG"
    svg_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating advanced haplotype visualizations for {dataset}...")
    
    gene_coords = load_gene_coordinates(gff3_path)
    manifest_samples = load_manifest(manifest_path, dataset)
    
    # 1. CliqueSNV stacked bar plot
    cliquesnv_dir = results_dir / "CliqueSNV"
    if cliquesnv_dir.exists():
        cs_df = load_cliquesnv_data(cliquesnv_dir, manifest_samples)
        draw_cliquesnv_plot(cs_df, plots_dir, svg_dir, dataset, no_normalize=args.no_normalize)
        
    # 2. VILOCA Local tracks
    viloca_dir = results_dir / "VILOCA"
    if viloca_dir.exists():
        generate_viloca_plots(viloca_dir, plots_dir, svg_dir, gene_coords, dataset, manifest_samples)
        
    # 3. SNPGenie Selection
    snpgenie_dir = results_dir / "SNPGenie"
    if snpgenie_dir.exists():
        generate_snpgenie_plots(snpgenie_dir, plots_dir, svg_dir, manifest_samples, dataset)

if __name__ == "__main__":
    main()
