#!/usr/bin/env python3
import sys
import os
import argparse
import glob
from pathlib import Path

__version__ = "1.0.0"

def main():
    parser = argparse.ArgumentParser(description="Export variant annotations and summary data to Excel workbook.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input-dir", required=True, help="Input directory with TSV/VCF files")
    parser.add_argument("--output-xlsx", required=True, help="Output Excel workbook file")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_path = Path(args.output_xlsx)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    files = list(input_dir.glob("*.tsv")) + list(input_dir.glob("*.vcf"))
    for file_path in files:
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 8:
                    info = parts[7]
                    af = 0.0
                    for item in info.split(";"):
                        if item.startswith("AF="):
                            try:
                                af = float(item.split("=")[1].split(",")[0])
                            except ValueError:
                                af = 0.0
                    af_pct = af * 100.0
                    if af_pct >= 1.0:
                        tier = ">=1%"
                    elif af_pct >= 0.1:
                        tier = "0.1-1% (candidate)"
                    else:
                        tier = "<0.1% (exploratory)"

                    rows.append({
                        "file": file_path.name,
                        "chrom": parts[0],
                        "pos": parts[1],
                        "ref": parts[3],
                        "alt": parts[4],
                        "af_proportion": f"{af:.6f}",
                        "af_percent": f"{af_pct:.2f}%",
                        "frequency_tier": tier
                    })

    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        df.to_excel(output_path, index=False)
        print(f"Exported {len(rows)} variant records to Excel workbook: {output_path}")
    except Exception:
        # Fallback if openpyxl / pandas not present
        csv_path = output_path.with_suffix(".csv")
        import csv
        with csv_path.open("w", encoding="utf-8", newline="") as f:
            if rows:
                writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)
        print(f"Exported {len(rows)} variant records to CSV fallback: {csv_path}")

if __name__ == "__main__":
    main()
