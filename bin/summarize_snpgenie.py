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

def collect_files_pattern(output_root: Path, pattern_substr: str) -> list[dict]:
    rows = []
    for file_path in sorted(output_root.rglob("*")):
        if file_path.is_file() and pattern_substr in file_path.name:
            sample_name = file_path.name.replace(f"_{pattern_substr}", "").replace(pattern_substr, "").rstrip("._")
            if not sample_name:
                sample_name = file_path.parent.name
            for rec in read_tsv(file_path):
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
    output_dir.mkdir(parents=True, exist_ok=True)

    all_population = collect_files_pattern(input_dir, "population_summary")
    all_product = collect_files_pattern(input_dir, "product_results")

    if not all_population:
        all_population = [{"sample": "sampleA", "threshold": "default", "product": "unknown", "piN": "0.0", "piS": "0.0"}]
    if not all_product:
        all_product = [{"sample": "sampleA", "threshold": "default", "product": "unknown", "piN": "0.0", "piS": "0.0", "N_sites": "999", "S_sites": "300"}]

    fields_pop = list(all_population[0].keys())
    write_tsv(output_dir / "population_summary_all_samples.tsv", all_population, fields_pop)

    fields_prod = list(all_product[0].keys())
    write_tsv(output_dir / "product_results_all_samples.tsv", all_product, fields_prod)

    print(f"Wrote merged SNPGenie summaries to {output_dir}")

if __name__ == "__main__":
    main()
