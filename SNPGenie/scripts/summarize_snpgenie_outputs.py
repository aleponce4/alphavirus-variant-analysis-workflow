#!/usr/bin/env python3

import csv
from pathlib import Path


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def collect_population(output_root: Path, threshold: str) -> list[dict]:
    rows = []
    for sample_dir in sorted(output_root.glob("*")):
        if not sample_dir.is_dir():
            continue
        pop_file = sample_dir / "population_summary.txt"
        if not pop_file.exists():
            continue
        for rec in read_tsv(pop_file):
            rec_out = {"threshold": threshold, "sample": sample_dir.name}
            rec_out.update(rec)
            rows.append(rec_out)
    return rows


def collect_product(output_root: Path, threshold: str) -> list[dict]:
    rows = []
    for sample_dir in sorted(output_root.glob("*")):
        if not sample_dir.is_dir():
            continue
        prod_file = sample_dir / "product_results.txt"
        if not prod_file.exists():
            continue
        for rec in read_tsv(prod_file):
            rec_out = {"threshold": threshold, "sample": sample_dir.name}
            rec_out.update(rec)
            rows.append(rec_out)
    return rows


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output_dir = root / "output"
    summary_dir = root / "summary"

    thresholds = {
        "minfreq_0p001": output_dir / "minfreq_0p001",
        "minfreq_0p01": output_dir / "minfreq_0p01",
    }

    all_population = []
    all_product = []
    for label, path in thresholds.items():
        all_population.extend(collect_population(path, label))
        all_product.extend(collect_product(path, label))

    if all_population:
        fields = list(all_population[0].keys())
        write_tsv(summary_dir / "population_summary_all_samples.tsv", all_population, fields)

    if all_product:
        fields = list(all_product[0].keys())
        write_tsv(summary_dir / "product_results_all_samples.tsv", all_product, fields)

    print(f"Wrote summaries to {summary_dir}")


if __name__ == "__main__":
    main()
