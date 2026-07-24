#!/usr/bin/env python3
"""
generate_viral_pct_vertical_plots.py

Generates simple publication-quality vertical barplots of % Viral Reads across datasets:
1. VEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B]
2. EEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B] (Excluding BDGR)
3. EEEV: BDGR 251 Antiviral Study (Untreated EEEV [Panel A] vs. BDGR 251 + EEEV [Panel B])

Features:
- Normal vertical bars (samples along X-axis, % Viral Reads on Y-axis).
- Standard Vermilion Red (#D55E00) viral palette.
- Percentage text callouts above viral bars.
- Clear DPI and Treatment group brackets along X-axis.
- 100% Black Arial typography (publication grade).
- Clean plot area with zero internal gridlines.
- Saves high-res PNG (300 DPI) and vector SVG formats.
"""

import os
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# Setup directories
BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent
RESULTS_DIR = BASE_DIR / "results"
MANIFEST_PATH = BASE_DIR / "config" / "samples_manifest.tsv"
LAB_PLOTS_DIR = RESULTS_DIR / "Lab_Meeting_Plots"
SVG_DIR = LAB_PLOTS_DIR / "SVG"

LAB_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
SVG_DIR.mkdir(parents=True, exist_ok=True)

# Matplotlib global style configuration
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9.0
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['savefig.dpi'] = 300

COLOR_VIRUS = '#D55E00'   # Standard Okabe-Ito Vermilion/Red
COLOR_MOCK_BG = '#F0F4F8'  # Soft Blue-Gray Shading for Mock
COLOR_ALT_BG = '#F8F9FA'   # Light Gray Shading for Alternating Groups

def parse_idxstats(idx_file_path, viral_chr):
    """Parses samtools idxstats MultiQC file to extract host and viral mapped reads per sample."""
    if not os.path.exists(idx_file_path):
        print(f"[WARNING] idxstats file missing: {idx_file_path}")
        return pd.DataFrame()
        
    df = pd.read_csv(idx_file_path, sep='\t')
    rows = []
    for _, r in df.iterrows():
        sample_raw = str(r['Sample'])
        bam_name = sample_raw.split('_')[0]
        
        host_mapped = 0
        viral_mapped = 0
        
        for col in df.columns:
            if col in ['Sample', '*']:
                continue
            val_str = str(r[col])
            mapped_count = ast.literal_eval(val_str)[0] if val_str.startswith('[') else int(r[col])
            
            if col == viral_chr:
                viral_mapped += mapped_count
            else:
                host_mapped += mapped_count
                
        total_mapped = host_mapped + viral_mapped
        viral_pct = (viral_mapped / total_mapped * 100) if total_mapped > 0 else 0.0
        
        rows.append({
            'bam_name': bam_name,
            'host_mapped': host_mapped,
            'viral_mapped': viral_mapped,
            'total_mapped': total_mapped,
            'viral_pct': viral_pct
        })
    return pd.DataFrame(rows)

def load_all_mapping_data():
    """Loads manifest and merges mapping stats for all datasets."""
    manifest = pd.read_csv(MANIFEST_PATH, sep='\t')
    
    veev_idx = BASE_DIR.parent / "nfcore_results" / "mouse_veev" / "multiqc" / "star_salmon" / "multiqc_report_data" / "multiqc_samtools_idxstats.txt"
    eeev_idx = BASE_DIR.parent / "nfcore_results" / "mouse_eeev" / "multiqc" / "star_salmon" / "multiqc_report_data" / "multiqc_samtools_idxstats.txt"
    
    df_veev = parse_idxstats(veev_idx, 'KP282671.1')
    df_eeev = parse_idxstats(eeev_idx, 'KP282670.1')
    
    mapping_all = pd.concat([df_veev, df_eeev], ignore_index=True)
    merged = pd.merge(manifest, mapping_all, on='bam_name', how='inner')
    return merged

def add_bottom_margin_brackets(ax, df_sorted):
    """Adds group background shading and bottom X-axis brackets for DPI and Treatment."""
    groups = []
    df_sorted['dpi_clean'] = "DPI " + df_sorted['dpi'].astype(str)
    df_sorted['group_key'] = df_sorted['treatment'].astype(str) + "___" + df_sorted['dpi_clean']
    
    current_key = None
    start_idx = 0
    
    for i, row in df_sorted.iterrows():
        k = row['group_key']
        if k != current_key:
            if current_key is not None:
                trt, dpi_str = current_key.split('___')
                groups.append({'treatment': trt, 'dpi': dpi_str, 'start_x': start_idx, 'end_x': i - 1})
            current_key = k
            start_idx = i
    if current_key is not None:
        trt, dpi_str = current_key.split('___')
        groups.append({'treatment': trt, 'dpi': dpi_str, 'start_x': start_idx, 'end_x': len(df_sorted) - 1})

    # 1. Background shading per DPI group
    for idx, g in enumerate(groups):
        start_x = g['start_x'] - 0.48
        end_x = g['end_x'] + 0.48
        bg_color = COLOR_MOCK_BG if 'Mock' in g['treatment'] else (COLOR_ALT_BG if idx % 2 == 0 else '#FFFFFF')
        ax.axvspan(start_x, end_x, color=bg_color, alpha=0.65, zorder=0)

    # 2. Add DPI Brackets below X-axis
    trans = ax.get_xaxis_transform()
    y_dpi_line = -0.17
    y_dpi_text = -0.20
    
    for g in groups:
        mid_x = (g['start_x'] + g['end_x']) / 2.0
        x_left = g['start_x'] - 0.35
        x_right = g['end_x'] + 0.35
        
        # Bracket line
        ax.plot([x_left, x_left, x_right, x_right],
                [y_dpi_line + 0.015, y_dpi_line, y_dpi_line, y_dpi_line + 0.015],
                color='black', lw=1.1, transform=trans, clip_on=False, zorder=10)
        
        # DPI Label
        ax.text(mid_x, y_dpi_text, g['dpi'], ha='center', va='top', fontsize=8.5, fontweight='bold', color='black', transform=trans, clip_on=False)

    # 3. Add Major Treatment Brackets further down
    treatment_spans = []
    cur_trt = None
    t_start = 0
    for i, g in enumerate(groups):
        if g['treatment'] != cur_trt:
            if cur_trt is not None:
                treatment_spans.append({'treatment': cur_trt, 'start_x': groups[t_start]['start_x'], 'end_x': groups[i-1]['end_x']})
            cur_trt = g['treatment']
            t_start = i
    if cur_trt is not None:
        treatment_spans.append({'treatment': cur_trt, 'start_x': groups[t_start]['start_x'], 'end_x': groups[-1]['end_x']})

    y_trt_line = -0.30
    y_trt_text = -0.34
    
    for t_span in treatment_spans:
        mid_x = (t_span['start_x'] + t_span['end_x']) / 2.0
        x_left = t_span['start_x'] - 0.40
        x_right = t_span['end_x'] + 0.40
        
        # Major Bracket Line
        ax.plot([x_left, x_left, x_right, x_right],
                [y_trt_line + 0.015, y_trt_line, y_trt_line, y_trt_line + 0.015],
                color='black', lw=1.4, transform=trans, clip_on=False, zorder=10)
        
        # Major Treatment Label
        ax.text(mid_x, y_trt_text, t_span['treatment'], ha='center', va='top', fontsize=9.5, fontweight='bold', color='black', transform=trans, clip_on=False)

def plot_viral_pct_vertical_dual(df_panel_a, df_panel_b, title_a, title_b, main_title, output_prefix, global_y_max=30.0):
    """Generates dual vertical barplots showing exclusively % Viral Reads."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 7.8), sharey=True)
    fig.suptitle(main_title, fontsize=15, fontweight='bold', y=0.98, color='black')
    
    panels = [(ax1, df_panel_a, title_a, 'A'), (ax2, df_panel_b, title_b, 'B')]
    
    for ax, df_sub, title_text, panel_letter in panels:
        if df_sub.empty:
            ax.text(0.5, 0.5, "No Data Available", ha='center', va='center', transform=ax.transAxes)
            continue
            
        # Sort samples logically from Left to Right
        df_sorted = df_sub.copy()
        df_sorted['treatment_rank'] = df_sorted['treatment'].apply(lambda x: 0 if 'Mock' in str(x) else (1 if 'alone' in str(x).lower() or 'challenge' in str(x).lower() or 'eeev' in str(x).lower() else 2))
        df_sorted = df_sorted.sort_values(by=['treatment_rank', 'dpi', 'sample_id'], ascending=[True, True, True]).reset_index(drop=True)
        
        x_positions = np.arange(len(df_sorted))
        viral_pcts = df_sorted['viral_pct'].values
        
        # Vertical Vermilion Bars
        bars = ax.bar(x_positions, viral_pcts, color=COLOR_VIRUS, edgecolor='black', linewidth=0.5, width=0.8, zorder=3)
        
        # Annotate percentages above vertical bars (only for detected viral reads > 0.01%)
        for i, v_pct in enumerate(viral_pcts):
            if v_pct >= 0.1:
                ax.text(x_positions[i], v_pct + (global_y_max * 0.02), f"{v_pct:.1f}%", ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='black')
            elif v_pct > 0.01:
                ax.text(x_positions[i], v_pct + (global_y_max * 0.02), f"{v_pct:.2f}%", ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='black')
                
        # Add Bottom Margin Brackets
        add_bottom_margin_brackets(ax, df_sorted)
        
        # Titles & Formatting
        ax.set_title(f"Panel {panel_letter}: {title_text}", fontsize=12, fontweight='bold', pad=14, color='black')
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"s{s}" for s in df_sorted['sample_id']], rotation=90, fontsize=8.5, fontweight='bold', color='black')
        ax.set_xlabel("Sample ID", fontsize=9.5, fontweight='bold', labelpad=2, color='black')
        
        if ax == ax1:
            ax.set_ylabel("Viral Reads (% of Mapped Output)", fontsize=10.5, fontweight='bold', labelpad=8, color='black')
            
        ax.set_ylim(0, global_y_max)
        ax.grid(False)
        
        # Style Spines in Pure Black
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_color('black')
            ax.spines[spine].set_linewidth(0.8)

    plt.tight_layout()
    fig.subplots_adjust(top=0.87, bottom=0.25, left=0.08, right=0.98, wspace=0.20)
    
    # Save PNG and SVG
    png_path = LAB_PLOTS_DIR / f"{output_prefix}.png"
    svg_path = SVG_DIR / f"{output_prefix}.svg"
    
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"  [OK] Saved simple vertical viral % figure: {png_path}")

def main():
    print("Loading alignment mapping stats...")
    df = load_all_mapping_data()
    
    # Global Max Viral Pct across dataset for shared Y scale
    max_viral_pct = df['viral_pct'].max()
    y_max = float(np.ceil(max_viral_pct / 5.0) * 5.0 + 3.0)  # ~30% max limit
    print(f"Dataset Max Viral %: {max_viral_pct:.2f}% -> Shared Y-axis limit: {y_max:.0f}%")
    
    # -------------------------------------------------------------
    # Figure 1: VEEV IN vs SC (% Viral Reads)
    # -------------------------------------------------------------
    print("\n1. Generating Vertical Figure 1: VEEV Intranasal vs. Subcutaneous (% Viral Reads)...")
    veev_in = df[(df['dataset'] == 'mouse_veev') & (df['route'] == 'intranasal')]
    veev_sc = df[(df['dataset'] == 'mouse_veev') & (df['route'] == 'SC')]
    
    plot_viral_pct_vertical_dual(
        df_panel_a=veev_in,
        df_panel_b=veev_sc,
        title_a="VEEV Intranasal (IN) — Study 045",
        title_b="VEEV Subcutaneous (SC) — Study 047",
        main_title="VEEV Viral Load Breakdown (% Viral Reads of Total Mapped)",
        output_prefix="veev_viral_pct_vertical_IN_vs_SC",
        global_y_max=y_max
    )
    
    # -------------------------------------------------------------
    # Figure 2: EEEV IN vs SC (% Viral Reads)
    # -------------------------------------------------------------
    print("\n2. Generating Vertical Figure 2: EEEV Intranasal vs. Subcutaneous (% Viral Reads)...")
    eeev_in = df[(df['dataset'] == 'mouse_eeev') & (df['route'] == 'intranasal')]
    eeev_sc = df[(df['dataset'] == 'mouse_eeev') & (df['route'] == 'SC') & (df['study_code'] == '048')]
    
    plot_viral_pct_vertical_dual(
        df_panel_a=eeev_in,
        df_panel_b=eeev_sc,
        title_a="EEEV Intranasal (IN) — Study 046",
        title_b="EEEV Subcutaneous (SC) — Study 048",
        main_title="EEEV Viral Load Breakdown (% Viral Reads of Total Mapped)",
        output_prefix="eeev_viral_pct_vertical_IN_vs_SC",
        global_y_max=y_max
    )
    
    # -------------------------------------------------------------
    # Figure 3: EEEV BDGR Antiviral Study (% Viral Reads)
    # -------------------------------------------------------------
    print("\n3. Generating Vertical Figure 3: EEEV BDGR 251 Antiviral Study (% Viral Reads)...")
    bdgr_all = df[(df['dataset'] == 'mouse_eeev') & (df['sample_id'].astype(str).str.len() >= 5)]
    
    bdgr_panel_a = bdgr_all[bdgr_all['treatment'].isin(['Mock', 'EEEV'])].copy()
    bdgr_panel_b = bdgr_all[bdgr_all['treatment'].isin(['Mock', 'BDGR 251 + EEEV'])].copy()
    
    plot_viral_pct_vertical_dual(
        df_panel_a=bdgr_panel_a,
        df_panel_b=bdgr_panel_b,
        title_a="EEEV Infection Alone (Untreated)",
        title_b="BDGR 251 Antiviral Treatment + EEEV",
        main_title="EEEV BDGR 251 Antiviral Study — Viral Load Breakdown (% Viral Reads)",
        output_prefix="eeev_bdgr_study_viral_pct_vertical",
        global_y_max=y_max
    )
    
    print("\n[OK] All simple vertical viral % figures generated successfully in:")
    print(f"     PNG: {LAB_PLOTS_DIR}")
    print(f"     SVG: {SVG_DIR}")

if __name__ == '__main__':
    main()
