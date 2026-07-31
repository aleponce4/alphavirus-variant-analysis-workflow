#!/usr/bin/env python3
"""
plot_haplotypes.py
Generate quasispecies haplotype reconstruction plots:
1. {dataset}_haplotype_frequencies.png (stacked bar chart, % y-axis)
2. {dataset}_haplotype_network.png (Minimum Spanning Network with treatment pie nodes)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

__version__ = "1.0.0"

## Okabe-Ito colorblind-friendly palette
OKABE_ITO = [
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish Purple
    "#999999",  # Grey
    "#000000",  # Black
]

PLT_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.title_fontsize": 8,
    "figure.dpi": 300,
}
plt.rcParams.update(PLT_STYLE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate publication-ready haplotype reconstruction plots."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--frequency-csv",
        help="Path to haplotype_frequency_by_sample.csv table.",
    )
    parser.add_argument(
        "--sequences-fasta",
        help="Path to haplotype_sequences.fasta file.",
    )
    parser.add_argument(
        "--samplesheet",
        help="Path to samplesheet CSV (sample,fastq_1,fastq_2,treatment).",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for generated PNG figures.",
    )
    parser.add_argument(
        "--dataset",
        default="viral_analysis",
        help="Dataset name prefix for output figure filenames.",
    )
    return parser.parse_args()


def load_samplesheet_metadata(samplesheet_path: Path | None) -> Dict[str, dict]:
    mapping: Dict[str, dict] = {}
    if not samplesheet_path or not samplesheet_path.exists():
        return mapping
    with samplesheet_path.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample_id = (row.get("sample") or row.get("sample_id") or "").strip()
            treatment = (row.get("treatment") or "unknown").strip()
            dpi_str = (row.get("dpi") or "0").strip()
            
            # Extract numeric dpi value
            m = re.search(r'\d+', dpi_str)
            dpi_val = int(m.group()) if m else 0
            
            if sample_id:
                mapping[sample_id] = {"treatment": treatment, "dpi": dpi_val}
    return mapping


def read_fasta_sequences(fasta_path: Path | None) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    if not fasta_path or not fasta_path.exists():
        return sequences

    current_id = ""
    seq_parts: List[str] = []
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(seq_parts).upper()
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_id:
            sequences[current_id] = "".join(seq_parts).upper()

    return sequences


def hamming_distance(seq1: str, seq2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2) if c1 != "N" and c2 != "N")



def generate_haplotype_network(
    df: pd.DataFrame,
    sequences: Dict[str, str],
    out_path: Path,
    dataset: str,
) -> None:
    """Generate Minimum Spanning Network (MSN) with Okabe-Ito treatment pie chart nodes (5 inches wide)."""
    if not HAS_NETWORKX:
        sys.stdout.write("Warning: networkx library not installed. Skipping haplotype network figure.\n")
        return

    if df.empty or not sequences:
        sys.stdout.write("Notice: No frequency data or sequences available. Skipping network plot.\n")
        return

    # Filter to >=1% tier in at least one sample
    qualifying_haps_set = set(df[df["frequency"] >= 0.01]["haplotype_id"].unique())
    qualifying_haps = [h for h in sorted(qualifying_haps_set) if h in sequences]

    # Explicitly anchor and label the Reference node
    ref_hap = None
    for h in qualifying_haps:
        sub = df[df["haplotype_id"] == h]
        if not sub.empty and "mutations" in sub.columns:
            mut_str = sub["mutations"].values[0]
            if pd.isna(mut_str) or mut_str == "REF" or mut_str == "":
                ref_hap = h
                break

    if ref_hap:
        sequences["Reference"] = sequences.pop(ref_hap)
        df.loc[df["haplotype_id"] == ref_hap, "haplotype_id"] = "Reference"
        qualifying_haps = [h if h != ref_hap else "Reference" for h in qualifying_haps]
    else:
        # If the reference isn't in qualifying_haps but exists in df with 0 freq, we could inject it
        pass

    if len(qualifying_haps) < 2:
        sys.stdout.write(f"Notice: Fewer than 2 qualifying haplotypes (>=1% tier) found ({len(qualifying_haps)}). Skipping network plot.\n")
        return

    def hamming_distance(s1: str, s2: str) -> int:
        min_l = min(len(s1), len(s2))
        return sum(1 for a, b in zip(s1[:min_l], s2[:min_l]) if a != 'N' and b != 'N' and a != b)

    # Calculate pairwise distances
    n_nodes = len(qualifying_haps)
    dist_matrix: Dict[Tuple[str, str], int] = {}
    graph = nx.Graph()

    for h in qualifying_haps:
        graph.add_node(h)

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            h1 = qualifying_haps[i]
            h2 = qualifying_haps[j]
            dist = hamming_distance(sequences[h1], sequences[h2])
            dist_matrix[(h1, h2)] = dist
            graph.add_edge(h1, h2, weight=dist)

    # Compute Minimum Spanning Tree (Kruskal) with deterministic tie-breaking
    mst = nx.minimum_spanning_tree(graph, weight="weight")

    # Compute layout
    pos = nx.spring_layout(mst, seed=42)

    # Determine grouping column (DPI or Treatment)
    if "dpi" in df.columns and len(df["dpi"].unique()) > 1:
        group_col = "dpi"
        groups = sorted(df[group_col].unique())
        legend_title = "DPI"
        group_labels = [f"DPI {g}" for g in groups]
    else:
        group_col = "treatment"
        groups = sorted(df[group_col].unique())
        legend_title = "Group"
        group_labels = [str(g) for g in groups]

    treat_color_map = {g: OKABE_ITO[i % len(OKABE_ITO)] for i, g in enumerate(groups)}

    node_treat_comp: Dict[str, List[float]] = {}
    node_abundance: Dict[str, float] = {}

    for h in qualifying_haps:
        sub = df[df["haplotype_id"] == h]
        
        # Calculate TOTAL abundance across all samples
        total_freq = sub["frequency"].sum() if not sub.empty else 0.01
        node_abundance[h] = total_freq

        comp = []
        for g in groups:
            # Sum of frequency for this specific group across samples
            g_sub = sub[sub[group_col] == g]
            g_freq = g_sub["frequency"].sum() if not g_sub.empty else 0.0
            comp.append(g_freq)
            
        sum_comp = sum(comp)
        if sum_comp > 0:
            comp = [c / sum_comp for c in comp]
        else:
            comp = [1.0 / len(groups)] * len(groups)
        node_treat_comp[h] = comp

    # Create 5-inch wide figure
    fig, ax = plt.subplots(figsize=(5, 4.5))

    # Draw edges with distance labels
    nx.draw_networkx_edges(mst, pos, ax=ax, edge_color="#94a3b8", width=1.2)
    edge_labels = {(u, v): f"{d['weight']} SNV" for u, v, d in mst.edges(data=True)}
    nx.draw_networkx_edge_labels(mst, pos, edge_labels=edge_labels, ax=ax, font_size=7, font_color="#475569")

    # Draw pie chart nodes using inset_axes
    base_radius = 0.05
    for node, (x, y) in pos.items():
        abundance = node_abundance[node]
        radius = base_radius * math.sqrt(abundance / 0.5 + 0.2)
        comp = node_treat_comp[node]

        # Draw inset pie
        ax_pie = ax.inset_axes([x - radius, y - radius, 2 * radius, 2 * radius], transform=ax.transData)
        ax_pie.pie(comp, colors=[treat_color_map[g] for g in groups], startangle=90)
        ax_pie.set_aspect("equal")

        # Label node
        ax.text(x, y + radius + 0.025, node, fontsize=7, fontweight="bold", ha="center", va="bottom", bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="#cbd5e1", alpha=0.85))

    # Legend for groups
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=treat_color_map[g]) for g in groups]
    ax.legend(legend_handles, group_labels, title=legend_title, loc="upper left", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=7, title_fontsize=8)

    ax.set_title(f"{dataset} — Haplotype Minimum Spanning Network (>=1% Tier)", pad=10, fontweight="bold", fontsize=9)
    ax.axis("off")

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    sys.stdout.write(f"Wrote network plot: {out_path}\n")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dataset = args.dataset
    freq_plot_path = out_dir / f"{dataset}_haplotype_frequencies.png"
    network_plot_path = out_dir / f"{dataset}_haplotype_network.png"

    # Load data
    df = pd.DataFrame()
    if args.frequency_csv and Path(args.frequency_csv).exists():
        try:
            df = pd.read_csv(args.frequency_csv)
        except Exception as e:
            sys.stderr.write(f"Warning: Failed to read frequency CSV {args.frequency_csv}: {e}\n")

    sequences = {}
    if args.sequences_fasta and Path(args.sequences_fasta).exists():
        sequences = read_fasta_sequences(Path(args.sequences_fasta))

    metadata = load_samplesheet_metadata(Path(args.samplesheet) if args.samplesheet else None)
    if not df.empty and metadata and "treatment" in df.columns:
        df["treatment"] = df["sample_id"].map(lambda s: metadata.get(s, {}).get("treatment", df.loc[df["sample_id"] == s, "treatment"].values[0] if "treatment" in df.columns else "unknown"))
        df["dpi"] = df["sample_id"].map(lambda s: metadata.get(s, {}).get("dpi", 0))

    generate_haplotype_network(df, sequences, network_plot_path, dataset)

    sys.stdout.write("Finished haplotype plotting.\n")


if __name__ == "__main__":
    main()
