#!/usr/bin/env python3
"""
plot_variants.py — Wrapper script delegating to generate_variant_plots.py
Ensures backwards compatibility while generating publication-grade plots and summary TSVs.
"""

import sys
import argparse
from pathlib import Path
from generate_variant_plots import parse_vcf, generate_plots, __version__
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Generate sliding-window variant plots and frequency tables.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--vcf-dir", required=True, help="Input directory containing annotated VCFs")
    parser.add_argument("--output-dir", "--outdir", dest="outdir", required=True, help="Output directory for plots")
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
    
    summary_tsv_path = outdir / "variant_frequency_summary_pct.tsv"
    if not df.empty:
        df.to_csv(summary_tsv_path, sep="\t", index=False)
    else:
        with summary_tsv_path.open("w", encoding="utf-8") as f:
            f.write("sample\tchrom\tpos\tref\talt\tdepth\taf_proportion\taf_percent\ttier\n")
            
    print(f"Processed {len(df)} variants across {len(vcf_files)} VCFs. Output TSV: {summary_tsv_path}")

    generate_plots(df, outdir, args.dataset)

if __name__ == "__main__":
    main()
