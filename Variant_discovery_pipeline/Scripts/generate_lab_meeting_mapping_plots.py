#!/usr/bin/env python3
"""
generate_lab_meeting_mapping_plots.py

Generates publication-quality horizontal stacked barplots (samples as rows) for lab meeting presentation:
1. VEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B]
2. EEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B] (Excluding BDGR)
3. EEEV: BDGR 251 Antiviral Study (Untreated EEEV [Panel A] vs. BDGR 251 + EEEV [Panel B])

Features:
- Samples as rows (horizontal bars, ax.barh) with bars close together (height=0.92, minimal gap).
- Bold Arial typography for titles, section headers, and percentage labels.
- Standard host (Okabe Blue #0072B2) vs. viral (Vermilion #D55E00) palette.
- Completely clean plot area with NO gridlines inside.
- Clear Y-axis group shading and right-hand annotations for Mock, DPI 1-4, Challenge, and Treatment.
- Saves high-res PNG (300 DPI) and vector SVG formats.
"""

import os
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import FuncFormatter
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
plt.rcParams['font.size'] = 8.5
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['savefig.dpi'] = 300

# Color Palette matching standard host vs virus theme
COLOR_HOST = '#0072B2'    # Standard Okabe-Ito Navy/Blue
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
        host_pct = (host_mapped / total_mapped * 100) if total_mapped > 0 else 0.0
        viral_pct = (viral_mapped / total_mapped * 100) if total_mapped > 0 else 0.0
        
        rows.append({
            'bam_name': bam_name,
            'host_mapped': host_mapped,
            'viral_mapped': viral_mapped,
            'total_mapped': total_mapped,
            'host_pct': host_pct,
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

def add_y_group_annotations(ax, df_sorted, x_max):
    """Adds background shading bands and group labels along the Y-axis rows."""
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
                groups.append({'treatment': trt, 'dpi': dpi_str, 'start_y': start_idx, 'end_y': i - 1})
            current_key = k
            start_idx = i
    if current_key is not None:
        trt, dpi_str = current_key.split('___')
        groups.append({'treatment': trt, 'dpi': dpi_str, 'start_y': start_idx, 'end_y': len(df_sorted) - 1})

    # 1. Background shading per group along Y rows
    for idx, g in enumerate(groups):
        start_y = g['start_y'] - 0.48
        end_y = g['end_y'] + 0.48
        bg_color = COLOR_MOCK_BG if 'Mock' in g['treatment'] else (COLOR_ALT_BG if idx % 2 == 0 else '#FFFFFF')
        ax.axhspan(start_y, end_y, color=bg_color, alpha=0.65, zorder=0)

    # 2. Add DPI Brackets & Annotations on the right margin
    for g in groups:
        mid_y = (g['start_y'] + g['end_y']) / 2.0
        ax.text(x_max * 1.01, mid_y, f"{g['treatment']} — {g['dpi']}", ha='left', va='center', fontsize=7.5, fontweight='bold', color='#1B365D')

def plot_horizontal_dual_panel(df_panel_a, df_panel_b, title_a, title_b, main_title, output_prefix):
    """Generates dual horizontal stacked barplots (1 row, 2 columns with samples as rows)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9.5), sharex=True)
    fig.suptitle(main_title, fontsize=15, fontweight='bold', y=0.985, color='#1B365D')
    
    panels = [(ax1, df_panel_a, title_a, 'A'), (ax2, df_panel_b, title_b, 'B')]
    
    # Calculate global max for X-axis scale
    max_reads_a = (df_panel_a['host_mapped'] + df_panel_a['viral_mapped']).max() if not df_panel_a.empty else 1e7
    max_reads_b = (df_panel_b['host_mapped'] + df_panel_b['viral_mapped']).max() if not df_panel_b.empty else 1e7
    global_x_max = max(max_reads_a, max_reads_b) * 1.25
    
    for ax, df_sub, title_text, panel_letter in panels:
        if df_sub.empty:
            ax.text(0.5, 0.5, "No Data Available", ha='center', va='center', transform=ax.transAxes)
            continue
            
        # Sort samples logically from TOP to BOTTOM (invert for barh so top sample is first)
        df_sorted = df_sub.copy()
        df_sorted['treatment_rank'] = df_sorted['treatment'].apply(lambda x: 0 if 'Mock' in str(x) else (1 if 'alone' in str(x).lower() or 'challenge' in str(x).lower() or 'eeev' in str(x).lower() else 2))
        df_sorted = df_sorted.sort_values(by=['treatment_rank', 'dpi', 'sample_id'], ascending=[True, True, True]).reset_index(drop=True)
        # Reverse dataframe so top sample is at y=0 (top of plot)
        df_sorted = df_sorted.iloc[::-1].reset_index(drop=True)
        
        y_positions = np.arange(len(df_sorted))
        host_reads = df_sorted['host_mapped'].values
        viral_reads = df_sorted['viral_mapped'].values
        viral_pcts = df_sorted['viral_pct'].values
        host_pcts = df_sorted['host_pct'].values
        
        # Plot Horizontal Stacked Bars (height=0.92 for minimal gap / tight rows)
        bars_host = ax.barh(y_positions, host_reads, color=COLOR_HOST, edgecolor='none', height=0.92, label='Mouse (Host) Reads', zorder=3)
        bars_viral = ax.barh(y_positions, viral_reads, left=host_reads, color=COLOR_VIRUS, edgecolor='none', height=0.92, label='Viral Reads', zorder=3)
        
        # Percentage Annotations inside & right of horizontal bars
        for i in range(len(df_sorted)):
            h = host_reads[i]
            v = viral_reads[i]
            v_pct = viral_pcts[i]
            h_pct = host_pcts[i]
            
            # Host % inside blue bar
            if h_pct > 70:
                ax.text(h * 0.5, y_positions[i], f"{h_pct:.1f}%", ha='center', va='center', color='white', fontsize=7.2, fontweight='bold')
                
            # Viral % inside or right of red bar
            if v_pct >= 0.5:
                x_pos = h + (v / 2.0) if v > global_x_max * 0.08 else (h + v + global_x_max * 0.012)
                text_col = 'white' if v > global_x_max * 0.08 else COLOR_VIRUS
                ax.text(x_pos, y_positions[i], f"{v_pct:.1f}%", ha='center' if text_col == 'white' else 'left', va='center', color=text_col, fontsize=7.2, fontweight='bold')
            elif v_pct > 0.01:
                ax.text(h + v + global_x_max * 0.01, y_positions[i], f"{v_pct:.2f}%", ha='left', va='center', color=COLOR_VIRUS, fontsize=6.8, fontweight='bold')
                
        # Add Y-axis background group shading
        add_y_group_annotations(ax, df_sorted, global_x_max)
        
        # Subplot Titles & Formatting
        ax.set_title(f"Panel {panel_letter}: {title_text}", fontsize=11.5, fontweight='bold', pad=18, color='#1B365D')
        ax.set_yticks(y_positions)
        ax.set_yticklabels([f"s{s}" for s in df_sorted['sample_id']], fontsize=8, fontweight='bold')
        ax.set_ylabel("Sample ID", fontsize=9.5, fontweight='bold')
        ax.set_xlabel("Mapped Read Count (Millions)", fontsize=9.5, fontweight='bold', labelpad=8)
        
        # REMOVE ALL GRIDLINES INSIDE PLOT AREA
        ax.grid(False)
        ax.set_xlim(0, global_x_max * 1.22)
        
        # Format X Axis to Millions (e.g., 20M, 40M)
        formatter = FuncFormatter(lambda x, _: f'{x*1e-6:.0f}M')
        ax.xaxis.set_major_formatter(formatter)

        # Style Spines
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_color('black')
            ax.spines[spine].set_linewidth(0.8)
    
    # Legend placed top right outside plot area
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.98, 0.96), frameon=True, facecolor='#F8F9FA', edgecolor='#CFD8DC', fontsize=8.8)
    
    plt.tight_layout()
    fig.subplots_adjust(top=0.88, bottom=0.10, wspace=0.38)
    
    # Save PNG and SVG
    png_path = LAB_PLOTS_DIR / f"{output_prefix}.png"
    svg_path = SVG_DIR / f"{output_prefix}.svg"
    
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"  [OK] Saved horizontal figure: {png_path}")

def main():
    print("Loading alignment mapping stats...")
    df = load_all_mapping_data()
    print(f"Total dataset records: {len(df)}")
    
    # -------------------------------------------------------------
    # Figure 1: VEEV IN vs SC (Horizontal Stacked Bars)
    # -------------------------------------------------------------
    print("\n1. Generating Horizontal Figure 1: VEEV Intranasal vs. Subcutaneous...")
    veev_in = df[(df['dataset'] == 'mouse_veev') & (df['route'] == 'intranasal')]
    veev_sc = df[(df['dataset'] == 'mouse_veev') & (df['route'] == 'SC')]
    
    plot_horizontal_dual_panel(
        df_panel_a=veev_in,
        df_panel_b=veev_sc,
        title_a="VEEV Intranasal (IN) — Study 045",
        title_b="VEEV Subcutaneous (SC) — Study 047",
        main_title="VEEV Host (Mouse) vs. Viral Read Mapping Breakdown",
        output_prefix="veev_reads_mapping_IN_vs_SC"
    )
    
    # -------------------------------------------------------------
    # Figure 2: EEEV IN vs SC (Horizontal Stacked Bars, Non-BDGR)
    # -------------------------------------------------------------
    print("\n2. Generating Horizontal Figure 2: EEEV Intranasal vs. Subcutaneous (Non-BDGR)...")
    eeev_in = df[(df['dataset'] == 'mouse_eeev') & (df['route'] == 'intranasal')]
    eeev_sc = df[(df['dataset'] == 'mouse_eeev') & (df['route'] == 'SC') & (df['study_code'] == '048')]
    
    plot_horizontal_dual_panel(
        df_panel_a=eeev_in,
        df_panel_b=eeev_sc,
        title_a="EEEV Intranasal (IN) — Study 046",
        title_b="EEEV Subcutaneous (SC) — Study 048",
        main_title="EEEV Host (Mouse) vs. Viral Read Mapping Breakdown",
        output_prefix="eeev_reads_mapping_IN_vs_SC"
    )
    
    # -------------------------------------------------------------
    # Figure 3: EEEV BDGR Antiviral Study (Horizontal Stacked Bars)
    # -------------------------------------------------------------
    print("\n3. Generating Horizontal Figure 3: EEEV BDGR 251 Antiviral Study...")
    bdgr_all = df[(df['dataset'] == 'mouse_eeev') & (df['sample_id'].astype(str).str.len() >= 5)]
    
    bdgr_panel_a = bdgr_all[bdgr_all['treatment'].isin(['Mock', 'EEEV'])].copy()
    bdgr_panel_b = bdgr_all[bdgr_all['treatment'].isin(['Mock', 'BDGR 251 + EEEV'])].copy()
    
    plot_horizontal_dual_panel(
        df_panel_a=bdgr_panel_a,
        df_panel_b=bdgr_panel_b,
        title_a="EEEV Infection Alone (Untreated)",
        title_b="BDGR 251 Antiviral Treatment + EEEV",
        main_title="EEEV BDGR 251 Antiviral Study — Host vs. Viral Read Mapping",
        output_prefix="eeev_bdgr_study_reads_mapping"
    )
    
    print("\n[OK] All horizontal lab meeting mapping figures generated successfully in:")
    print(f"     PNG: {LAB_PLOTS_DIR}")
    print(f"     SVG: {SVG_DIR}")

if __name__ == '__main__':
    main()
