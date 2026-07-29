#!/usr/bin/env python3
import sys
import os
import argparse
import glob
from pathlib import Path

__version__ = "1.0.0"

def parse_vcf_af_percentages(vcf_path: str) -> list[dict]:
    """Parse VCF and convert machine allele frequency (proportion 0-1) to human percentage (0-100%)."""
    records = []
    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            info = parts[7]
            af = 0.0
            for item in info.split(';'):
                if item.startswith('AF='):
                    try:
                        af = float(item.split('=')[1].split(',')[0])
                    except ValueError:
                        af = 0.0
            records.append({
                'chrom': parts[0],
                'pos': int(parts[1]),
                'ref': parts[3],
                'alt': parts[4],
                'af_proportion': af,
                'af_percent': af * 100.0
            })
    return records

def main():
    parser = argparse.ArgumentParser(description="Generate sliding-window variant plots.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--vcf-dir", required=True, help="Input directory containing annotated VCFs")
    parser.add_argument("--output-dir", required=True, help="Output directory for plots")

    args = parser.parse_args()

    vcf_dir = Path(args.vcf_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    vcfs = list(vcf_dir.glob("*.vcf")) + list(vcf_dir.glob("*.vcf.gz"))
    summary_path = out_dir / "variant_frequency_summary_pct.tsv"

    total_records = 0
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("sample\tchrom\tpos\tref\talt\taf_proportion\taf_percent\ttier\n")
        for vcf in vcfs:
            sample = vcf.name.split('.')[0]
            records = parse_vcf_af_percentages(str(vcf))
            total_records += len(records)
            for r in records:
                af_pct = r['af_percent']
                if af_pct >= 1.0:
                    tier = ">=1%"
                elif af_pct >= 0.1:
                    tier = "0.1-1% (candidate)"
                else:
                    tier = "<0.1% (exploratory)"
                f.write(f"{sample}\t{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['af_proportion']:.6f}\t{af_pct:.2f}%\t{tier}\n")

    print(f"Processed {total_records} variants across {len(vcfs)} VCFs. Output: {summary_path}")

if __name__ == "__main__":
    main()
