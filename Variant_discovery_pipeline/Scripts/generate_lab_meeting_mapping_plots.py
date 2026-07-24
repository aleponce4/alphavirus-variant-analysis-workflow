#!/usr/bin/env python3
"""
generate_lab_meeting_mapping_plots.py

Generates publication-quality horizontal stacked barplots (samples as rows) for lab meeting presentation:
1. VEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B]
2. EEEV: Intranasal (IN) [Panel A] vs. Subcutaneous (SC) [Panel B] (Excluding BDGR)
3. EEEV: BDGR 251 Antiviral Study (Untreated EEEV [Panel A] vs. BDGR 251 + EEEV [Panel B])

Features:
- Samples as rows (horizontal bars, ax.barh) with bars close together (height=0.92, minimal gap).
- Horizontal "Sample ID" header placed cleanly at the top of the sample column (above top sample tick).
- Vertically ALIGNED host percentage numbers inside blue bars (fixed left-aligned position x=12M).
- Shared Read Count X-axis scale (0 to 150 Million Reads) across ALL figures for direct comparability.
- Left-side group brackets (outside plot margin) for DPI and Treatment categories.
- 100% Black Arial typography for all titles, labels, and ticks (publication grade).
- Standard host (Okabe Blue #0072B2) vs. viral (Vermilion #D55E00) palette.
- Completely clean plot area with NO gridlines inside.
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
plt.rcParams['font.size'] = 9.0
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

def add_left_margin_brackets(ax, df_sorted):
    """Adds background row shading and clean left-side outside brackets for DPI and Treatment."""
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

    # 2. Add DPI Brackets on Left Margin (outside Y axis)
    trans = ax.get_yaxis_transform()
    
    x_dpi_line = -0.16
    x_dpi_text = -0.18
    
    for g in groups:
        mid_y = (g['start_y'] + g['end_y']) / 2.0
        y_top = g['start_y'] - 0.35
        y_bot = g['end_y'] + 0.35
        
        # Bracket line
        ax.plot([x_dpi_line + 0.02, x_dpi_line, x_dpi_line, x_dpi_line + 0.02],
                [y_top, y_top, y_bot, y_bot],
                color='black', lw=1.1, transform=trans, clip_on=False, zorder=10)
        
        # DPI Label
        ax.text(x_dpi_text, mid_y, g['dpi'], ha='right', va='center', fontsize=8.5, fontweight='bold', color='black', transform=trans, clip_on=False)

    # 3. Add Major Treatment Brackets further left
    treatment_spans = []
    cur_trt = None
    t_start = 0
    for i, g in enumerate(groups):
        if g['treatment'] != cur_trt:
            if cur_trt is not None:
                treatment_spans.append({'treatment': cur_trt, 'start_y': groups[t_start]['start_y'], 'end_y': groups[i-1]['end_y']})
            cur_trt = g['treatment']
            t_start = i
    if cur_trt is not None:
        treatment_spans.append({'treatment': cur_trt, 'start_y': groups[t_start]['start_y'], 'end_y': groups[-1]['end_y']})

    x_trt_line = -0.32
    x_trt_text = -0.35
    
    for t_span in treatment_spans:
        mid_y = (t_span['start_y'] + t_span['end_y']) / 2.0
        y_top = t_span['start_y'] - 0.40
        y_bot = t_span['end_y'] + 0.40
        
        # Major Bracket Line
        ax.plot([x_trt_line + 0.02, x_trt_line, x_trt_line, x_trt_line + 0.02],
                [y_top, y_top, y_bot, y_bot],
                color='black', lw=1.4, transform=trans, clip_on=False, zorder=10)
        
        # Major Treatment Label (Rotated vertical for compact fit)
        ax.text(x_trt_text, mid_y, t_span['treatment'], ha='center', va='center', rotation=90, fontsize=9.5, fontweight='bold', color='black', transform=trans, clip_on=False)

def plot_horizontal_dual_panel(df_panel_a, df_panel_b, title_a, title_b, main_title, output_prefix, global_x_max=150e6):
    """Generates dual horizontal stacked barplots with identical shared X-axis max (150M reads) across all figures."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17.5, 9.5), sharex=True)
    fig.suptitle(main_title, fontsize=15, fontweight='bold', y=0.988, color='black')
    
    panels = [(ax1, df_panel_a, title_a, 'A'), (ax2, df_panel_b, title_b, 'B')]
    
    for ax, df_sub, title_text, panel_letter in panels:
        if df_sub.empty:
            ax.text(0.5, 0.5, "No Data Available", ha='center', va='center', transform=ax.transAxes)
            continue
            
        # Sort samples logically from TOP to BOTTOM
        df_sorted = df_sub.copy()
        df_sorted['treatment_rank'] = df_sorted['treatment'].apply(lambda x: 0 if 'Mock' in str(x) else (1 if 'alone' in str(x).lower() or 'challenge' in str(x).lower() or 'eeev' in str(x).lower() else 2))
        df_sorted = df_sorted.sort_values(by=['treatment_rank', 'dpi', 'sample_id'], ascending=[True, True, True]).reset_index(drop=True)
        # Reverse dataframe so top sample is at y=0
        df_sorted = df_sorted.iloc[::-1].reset_index(drop=True)
        
        y_positions = np.arange(len(df_sorted))
        host_reads = df_sorted['host_mapped'].values
        viral_reads = df_sorted['viral_mapped'].values
        viral_pcts = df_sorted['viral_pct'].values
        host_pcts = df_sorted['host_pct'].values
        
        # Plot Horizontal Stacked Bars (height=0.92 for minimal gap)
        bars_host = ax.barh(y_positions, host_reads, color=COLOR_HOST, edgecolor='none', height=0.92, label='Mouse (Host) Reads', zorder=3)
        bars_viral = ax.barh(y_positions, viral_reads, left=host_reads, color=COLOR_VIRUS, edgecolor='none', height=0.92, label='Viral Reads', zorder=3)
        
        # Percentage Annotations inside & right of horizontal bars (PERFECTLY ALIGNED Host % at x = 12M)
        ALIGNED_HOST_X = 12e6  # Fixed X coordinate for 100% vertical alignment inside blue bars
        
        for i in range(len(df_sorted)):
            h = host_reads[i]
            v = viral_reads[i]
            v_pct = viral_pcts[i]
            h_pct = host_pcts[i]
            
            # Host % PERFECTLY VERTICALLY ALIGNED inside blue bar at x=12M
            if h_pct > 50:
                ax.text(ALIGNED_HOST_X, y_positions[i], f"{h_pct:.1f}%", ha='left', va='center', color='white', fontsize=7.8, fontweight='bold')
                
            # Viral % inside or right of red bar
            if v_pct >= 0.5:
                x_pos = h + (v / 2.0) if v > global_x_max * 0.08 else (h + v + global_x_max * 0.015)
                text_col = 'white' if v > global_x_max * 0.08 else 'black'
                ax.text(x_pos, y_positions[i], f"{v_pct:.1f}%", ha='center' if text_col == 'white' else 'left', va='center', color=text_col, fontsize=7.5, fontweight='bold')
            elif v_pct > 0.01:
                ax.text(h + v + global_x_max * 0.012, y_positions[i], f"{v_pct:.2f}%", ha='left', va='center', color='black', fontsize=7.0, fontweight='bold')
                
        # Add Left Margin Brackets (outside Y-axis)
        add_left_margin_brackets(ax, df_sorted)
        
        # Subplot Titles & Formatting
        ax.set_title(f"Panel {panel_letter}: {title_text}", fontsize=12, fontweight='bold', pad=24, color='black')
        ax.set_yticks(y_positions)
        ax.set_yticklabels([f"s{s}" for s in df_sorted['sample_id']], fontsize=8.5, fontweight='bold', color='black')
        
        # Move "Sample ID" label to TOP, HORIZONTAL, directly above the top sample ID tick (s101 / s301)
        ax.set_ylabel("")  # Clear standard vertical Y label
        trans = ax.get_yaxis_transform()
        top_y_pos = len(df_sorted) - 0.15
        ax.text(-0.02, top_y_pos, "Sample ID", ha='right', va='bottom', fontsize=9.0, fontweight='bold', color='black', transform=trans, clip_on=False)
        
        ax.set_xlabel("Mapped Read Count (Millions)", fontsize=10, fontweight='bold', labelpad=8, color='black')
        
        # REMOVE ALL GRIDLINES INSIDE PLOT AREA
        ax.grid(False)
        ax.set_xlim(0, global_x_max)
        
        # Format X Axis to Millions (0M, 25M, 50M... 150M)
        formatter = FuncFormatter(lambda x, _: f'{x*1e-6:.0f}M')
        ax.xaxis.set_major_formatter(formatter)

        # Style Spines in Pure Black
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_color('black')
            ax.spines[spine].set_linewidth(0.8)
    
    # Global Legend placed top right above plot area
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.985, 0.98), frameon=True, facecolor='#F8F9FA', edgecolor='#CFD8DC', fontsize=9.0)
    
    plt.tight_layout()
    fig.subplots_adjust(top=0.87, bottom=0.10, left=0.18, wspace=0.45)
    
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
    
    # Calculate universal max read count across ALL 130 samples to set identical X-scale
    max_total_reads = (df['host_mapped'] + df['viral_mapped']).max()
    universal_x_max = 150e6  # 150 Million Reads shared scale for 100% direct comparability across figures
    print(f"Dataset Max Reads: {max_total_reads*1e-6:.2f}M -> Shared X-axis limit: {universal_x_max*1e-6:.0f}M")
    
    # -------------------------------------------------------------
    # Figure 1: VEEV IN vs SC (Shared 150M Scale, Left Brackets)
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
        output_prefix="veev_reads_mapping_IN_vs_SC",
        global_x_max=universal_x_max
    )
    
    # -------------------------------------------------------------
    # Figure 2: EEEV IN vs SC (Shared 150M Scale, Left Brackets)
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
        output_prefix="eeev_reads_mapping_IN_vs_SC",
        global_x_max=universal_x_max
    )
    
    # -------------------------------------------------------------
    # Figure 3: EEEV BDGR Antiviral Study (Shared 150M Scale, Left Brackets)
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
        output_prefix="eeev_bdgr_study_reads_mapping",
        global_x_max=universal_x_max
    )
    
    print("\n[OK] All horizontal lab meeting mapping figures generated successfully in:")
    print(f"     PNG: {LAB_PLOTS_DIR}")
    print(f"     SVG: {SVG_DIR}")

if __name__ == '__main__':
    main()
