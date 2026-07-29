#!/usr/bin/env python3
import sys
import os
import csv
import argparse
from pathlib import Path

__version__ = "1.0.0"

def read_tsv(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def collect_files(output_root: Path, filename: str) -> list[dict]:
    rows = []
    for sample_dir in sorted(output_root.glob("*")):
        if sample_dir.is_dir():
            target_file = sample_dir / filename
            if target_file.exists():
                for rec in read_tsv(target_file):
                    rec_out = {"sample": sample_dir.name}
                    rec_out.update(rec)
                    rows.append(rec_out)
        elif sample_dir.is_file() and sample_dir.name.endswith(filename):
            sample_name = sample_dir.name.replace(f"_{filename}", "").replace(filename, "")
            for rec in read_tsv(sample_dir):
                rec_out = {"sample": sample_name}
                rec_out.update(rec)
                rows.append(rec_out)
    return rows

def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-sample SNPGenie output files.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input-dir", required=True, help="SNPGenie output folder (contains per-sample results)")
    parser.add_argument("--output-dir", required=True, help="Output folder to save merged summaries")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    all_population = collect_files(input_dir, "population_summary.txt")
    if not all_population:
        all_population = collect_files(input_dir, "results.tsv")

    all_product = collect_files(input_dir, "product_results.txt")
    if not all_product:
        all_product = collect_files(input_dir, "product_results.tsv")

    if all_population:
        fields = list(all_population[0].keys())
        write_tsv(output_dir / "population_summary_all_samples.tsv", all_population, fields)

    if all_product:
        fields = list(all_product[0].keys())
        write_tsv(output_dir / "product_results_all_samples.tsv", all_product, fields)

    print(f"Wrote merged SNPGenie summaries to {output_dir}")

if __name__ == "__main__":
    main()
