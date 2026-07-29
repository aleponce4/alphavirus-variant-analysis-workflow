#!/usr/bin/env python3
import sys
import os
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

__version__ = "1.0.0"

try:
    from scipy.stats import kruskal
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute delta (piN-piS) and gene-wise Kruskal-Wallis with BH-FDR using manifest DPIs."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input", required=True, help="Input product summary TSV")
    parser.add_argument("--outdir", required=True, help="Output directory for delta analysis tables")
    parser.add_argument("--manifest", required=True, help="Path to samples manifest CSV or TSV")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR significance cutoff (default: 0.05)")
    parser.add_argument("--weak-epsilon", type=float, default=1e-6, help="Absolute delta threshold for piN≈piS classification (default: 1e-6)")
    parser.add_argument("--write-detailed", action="store_true", help="Also write per-threshold Kruskal tables.")
    return parser.parse_args()

def load_manifest(manifest_path: Path) -> dict[str, str]:
    mapping = {}
    with manifest_path.open("r", encoding="utf-8") as f:
        delim = "," if manifest_path.suffix.lower() == ".csv" else "\t"
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            sample_id = row.get("sample") or row.get("bam_name")
            dpi = row.get("dpi") or row.get("treatment") or row.get("condition") or "group"
            if sample_id:
                mapping[sample_id] = f"dpi{dpi}" if str(dpi).isdigit() else str(dpi)
    return mapping

def safe_float(value: str) -> float:
    if value in {"", "*", "NA"}:
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")

def bh_fdr(pvals: list[float]) -> list[float]:
    n = len(pvals)
    if n == 0:
        return []
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

def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)

    if not manifest_path.exists():
        print(f"Manifest not found at {manifest_path}", file=sys.stderr)
        sys.exit(1)

    manifest_map = load_manifest(manifest_path)
    all_dpis = sorted(list(set(manifest_map.values())))
    if not all_dpis:
        all_dpis = ["dpi0"]

    per_sample_rows = []
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    if in_path.exists():
        with in_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample = row.get("sample", "")
                if sample not in manifest_map:
                    continue

                dpi = manifest_map[sample]
                threshold = row.get("threshold", "default")
                product = row.get("product", "unknown")
                pi_n = safe_float(row.get("piN", "NA"))
                pi_s = safe_float(row.get("piS", "NA"))

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
        fieldnames = ["threshold", "sample", "dpi", "product", "piN", "piS", "delta_piN_minus_piS", "selection_signal"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_sample_rows)

    summary_rows = []
    for threshold, product_map in grouped.items():
        interim = []
        for product, dpi_map in product_map.items():
            dpi_groups = [dpi_map[d] for d in all_dpis if len(dpi_map[d]) > 0]

            if len(dpi_groups) < 2 or not HAS_SCIPY:
                h_stat = float("nan")
                pval = float("nan")
            else:
                try:
                    all_vals = [v for g in dpi_groups for v in g]
                    if len(set(all_vals)) <= 1:
                        h_stat, pval = 0.0, 1.0
                    else:
                        h_stat, pval = kruskal(*dpi_groups)
                except Exception:
                    h_stat, pval = float("nan"), float("nan")

            row = {
                "threshold": threshold,
                "product": product,
                "kruskal_H": format_float(h_stat),
                "kruskal_p": format_float(pval),
            }
            for d in all_dpis:
                vals = dpi_map.get(d, [])
                row[f"n_{d}"] = str(len(vals))
                row[f"median_delta_{d}"] = format_float(median(vals))
                row[f"signal_{d}"] = interpret_delta(median(vals), args.weak_epsilon)

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
        fieldnames = ["threshold", "product"]
        for d in all_dpis:
            fieldnames.append(f"n_{d}")
        for d in all_dpis:
            fieldnames.append(f"median_delta_{d}")
        for d in all_dpis:
            fieldnames.append(f"signal_{d}")
        fieldnames.extend(["kruskal_H", "kruskal_p", "bh_fdr_q", "significant_fdr"])

        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote delta per-sample table: {per_sample_path}")
    print(f"Wrote Kruskal-Wallis summary table: {summary_path}")

if __name__ == "__main__":
    main()
