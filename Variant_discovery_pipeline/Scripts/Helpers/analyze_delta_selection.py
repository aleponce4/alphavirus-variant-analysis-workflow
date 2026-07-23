#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

try:
    from scipy.stats import kruskal
except ImportError as exc:
    raise SystemExit("This script requires scipy. Install with: pip install scipy") from exc

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute delta (piN-piS) and gene-wise Kruskal-Wallis with BH-FDR using manifest DPIs."
    )
    parser.add_argument("--input", required=True, help="Input product summary TSV")
    parser.add_argument("--outdir", required=True, help="Output directory for delta analysis tables")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR significance cutoff (default: 0.05)")
    parser.add_argument("--weak-epsilon", type=float, default=1e-6, help="Absolute delta threshold for piN≈piS classification (default: 1e-6)")
    parser.add_argument("--write-detailed", action="store_true", help="Also write per-threshold Kruskal tables.")
    return parser.parse_args()

def load_manifest(manifest_path: Path) -> dict[str, str]:
    """Maps bam_name -> dpi (e.g. s101 -> dpi1)"""
    mapping = {}
    with manifest_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            bam_name = row.get("bam_name")
            dpi = row.get("dpi")
            if bam_name and dpi:
                mapping[bam_name] = f"dpi{dpi}"
    return mapping

def safe_float(value: str) -> float:
    if value in {"", "*", "NA"}:
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

def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)

    if not manifest_path.exists():
        raise SystemExit(f"Manifest not found at {manifest_path}")

    # Load manifest mapping: bam_name -> dpi
    manifest_map = load_manifest(manifest_path)
    all_dpis = sorted(list(set(manifest_map.values())))
    print(f"Loaded manifest. Unique DPIs found: {', '.join(all_dpis)}")

    per_sample_rows = []
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    with in_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            if sample not in manifest_map:
                continue

            dpi = manifest_map[sample]
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
        fieldnames = ["threshold", "sample", "dpi", "product", "piN", "piS", "delta_piN_minus_piS", "selection_signal"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_sample_rows)

    summary_rows = []
    for threshold, product_map in grouped.items():
        interim = []
        for product, dpi_map in product_map.items():
            dpi_groups = [dpi_map[d] for d in all_dpis if len(dpi_map[d]) > 0]
            
            # Kruskal requires at least 2 groups and at least 2 samples per group for some implementations,
            # but scipy kruskal requires at least 2 groups with at least one element.
            if len(dpi_groups) < 2:
                h_stat = float("nan")
                pval = float("nan")
            else:
                try:
                    # check if all values are identical across all groups
                    all_vals = []
                    for g in dpi_groups:
                        all_vals.extend(g)
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
            # Dynamically add n, median, and signal for all DPIs
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
        # Fieldnames: threshold, product, n_dpi1, n_dpi2..., median_delta_dpi1..., signal_dpi1..., kruskal_H, kruskal_p, bh_fdr_q, significant_fdr
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
