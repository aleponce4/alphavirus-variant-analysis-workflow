#!/usr/bin/env python3

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


try:
    from scipy.stats import kruskal
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires scipy. Install with: pip install scipy"
    ) from exc


DPI_PATTERN = re.compile(r"INH_(\d+)_DPI_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute delta (piN-piS) and gene-wise Kruskal-Wallis with BH-FDR by threshold."
    )
    parser.add_argument(
        "--input",
        default="SNPGenie/summary/product_results_all_samples.tsv",
        help="Input product summary TSV",
    )
    parser.add_argument(
        "--outdir",
        default="SNPGenie/analysis/delta_selection",
        help="Output directory for delta analysis tables",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR significance cutoff (default: 0.05)",
    )
    parser.add_argument(
        "--weak-epsilon",
        type=float,
        default=1e-6,
        help="Absolute delta threshold for piN≈piS classification (default: 1e-6)",
    )
    parser.add_argument(
        "--write-detailed",
        action="store_true",
        help="Also write per-threshold Kruskal tables.",
    )
    return parser.parse_args()


def parse_dpi(sample: str) -> str | None:
    match = DPI_PATTERN.search(sample)
    if not match:
        return None
    dpi = match.group(1)
    return f"dpi{dpi}"


def safe_float(value: str) -> float:
    if value in {"", "*"}:
        return float("nan")
    return float(value)


def bh_fdr(pvals: list[float]) -> list[float]:
    n = len(pvals)
    ranked = sorted(enumerate(pvals), key=lambda x: x[1])
    qvals = [1.0] * n

    min_adj = 1.0
    for rank in range(n, 0, -1):
        idx, pval = ranked[rank - 1]
        adj = pval * n / rank
        if adj < min_adj:
            min_adj = adj
        qvals[idx] = min(min_adj, 1.0)
    return qvals


def median(values: list[float]) -> float:
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2 == 1:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


def interpret_delta(delta_value: float, epsilon: float) -> str:
    if math.isnan(delta_value):
        return "NA"
    if delta_value < -epsilon:
        return "purifying_selection_signal"
    if abs(delta_value) <= epsilon:
        return "weak_constraint_signal"
    return "possible_adaptive_pressure_signal"


def format_float(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.12g}"


def run_kruskal_three_groups(g1: list[float], g2: list[float], g3: list[float]) -> tuple[float, float]:
    merged = g1 + g2 + g3
    if len(merged) == 0:
        return float("nan"), float("nan")

    if len(set(merged)) == 1:
        return 0.0, 1.0

    h_stat, pval = kruskal(g1, g2, g3)
    return float(h_stat), float(pval)


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_sample_rows = []
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    with in_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            if "pilot" in sample.lower():
                continue

            dpi = parse_dpi(sample)
            if dpi not in {"dpi1", "dpi3", "dpi5"}:
                continue

            threshold = row["threshold"]
            product = row["product"]
            pi_n = safe_float(row["piN"])
            pi_s = safe_float(row["piS"])

            if math.isnan(pi_n) or math.isnan(pi_s):
                continue

            delta = pi_n - pi_s
            per_sample_rows.append(
                {
                    "threshold": threshold,
                    "sample": sample,
                    "dpi": dpi,
                    "product": product,
                    "piN": format_float(pi_n),
                    "piS": format_float(pi_s),
                    "delta_piN_minus_piS": format_float(delta),
                    "selection_signal": interpret_delta(delta, args.weak_epsilon),
                }
            )
            grouped[threshold][product][dpi].append(delta)

    per_sample_rows.sort(key=lambda r: (r["threshold"], r["product"], r["dpi"], r["sample"]))

    per_sample_path = out_dir / "delta_per_sample.tsv"
    with per_sample_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "threshold",
            "sample",
            "dpi",
            "product",
            "piN",
            "piS",
            "delta_piN_minus_piS",
            "selection_signal",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_sample_rows)

    summary_rows = []
    for threshold, product_map in grouped.items():
        interim = []
        for product, dpi_map in product_map.items():
            g1 = dpi_map.get("dpi1", [])
            g3 = dpi_map.get("dpi3", [])
            g5 = dpi_map.get("dpi5", [])

            all_groups = [g1, g3, g5]
            if any(len(g) == 0 for g in all_groups):
                h_stat = float("nan")
                pval = float("nan")
            else:
                h_stat, pval = run_kruskal_three_groups(g1, g3, g5)

            row = {
                "threshold": threshold,
                "product": product,
                "n_dpi1": str(len(g1)),
                "n_dpi3": str(len(g3)),
                "n_dpi5": str(len(g5)),
                "median_delta_dpi1": format_float(median(g1)),
                "median_delta_dpi3": format_float(median(g3)),
                "median_delta_dpi5": format_float(median(g5)),
                "signal_dpi1": interpret_delta(median(g1), args.weak_epsilon),
                "signal_dpi3": interpret_delta(median(g3), args.weak_epsilon),
                "signal_dpi5": interpret_delta(median(g5), args.weak_epsilon),
                "kruskal_H": format_float(h_stat),
                "kruskal_p": format_float(pval),
            }
            interim.append((row, pval))

        valid_pvals = [p for _, p in interim if not math.isnan(p)]
        qvals = bh_fdr(valid_pvals) if valid_pvals else []
        qidx = 0
        for row, pval in interim:
            if math.isnan(pval):
                qval = float("nan")
                sig = "NA"
            else:
                qval = qvals[qidx]
                qidx += 1
                sig = "TRUE" if qval <= args.alpha else "FALSE"

            row["bh_fdr_q"] = format_float(qval)
            row["significant_fdr"] = sig
            summary_rows.append(row)

    summary_rows.sort(key=lambda r: (r["threshold"], r["product"]))

    summary_path = out_dir / "delta_kruskal_by_gene.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "threshold",
            "product",
            "n_dpi1",
            "n_dpi3",
            "n_dpi5",
            "median_delta_dpi1",
            "median_delta_dpi3",
            "median_delta_dpi5",
            "signal_dpi1",
            "signal_dpi3",
            "signal_dpi5",
            "kruskal_H",
            "kruskal_p",
            "bh_fdr_q",
            "significant_fdr",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    if args.write_detailed:
        for threshold in sorted({row["threshold"] for row in summary_rows}):
            threshold_rows = [r for r in summary_rows if r["threshold"] == threshold]
            per_threshold = out_dir / f"delta_kruskal_by_gene_{threshold}.tsv"
            with per_threshold.open("w", encoding="utf-8", newline="") as handle:
                fieldnames = list(threshold_rows[0].keys()) if threshold_rows else []
                if fieldnames:
                    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(threshold_rows)

    print(f"Wrote per-sample delta table: {per_sample_path}")
    print(f"Wrote Kruskal/BH table: {summary_path}")


if __name__ == "__main__":
    main()
