#!/usr/bin/env python3

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


DPI_PATTERN = re.compile(r"INH_(\d+)_DPI_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build compact SNPGenie selection tables and a limma-focused key gene report."
    )
    parser.add_argument(
        "--base",
        default="SNPGenie/analysis/delta_selection",
        help="Directory containing delta/limma outputs (default: SNPGenie/analysis/delta_selection)",
    )
    parser.add_argument(
        "--product-summary",
        default="SNPGenie/summary/product_results_all_samples.tsv",
        help="Path to merged product summary TSV.",
    )
    parser.add_argument(
        "--threshold",
        default="minfreq_0p01",
        help="Threshold label to use for key table (default: minfreq_0p01)",
    )
    parser.add_argument(
        "--report-dir",
        default=None,
        help="Directory to write compact/key output tables (default: same as --base).",
    )
    parser.add_argument(
        "--direction-epsilon",
        type=float,
        default=1e-6,
        help="Absolute delta cutoff for neutral direction classification.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_dpi(sample: str) -> str | None:
    match = DPI_PATTERN.search(sample)
    if not match:
        return None
    return f"dpi{match.group(1)}"


def safe_float(value: str) -> float:
    if value in {"", "*", "NA"}:
        return float("nan")
    return float(value)


def fmt(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.12g}"


def mean(values: list[float]) -> float:
    if not values:
        return float("nan")
    return sum(values) / len(values)


def classify_direction(delta_value: float, epsilon: float) -> str:
    if math.isnan(delta_value):
        return "NA"
    if delta_value > epsilon:
        return "positive (piN>piS)"
    if delta_value < -epsilon:
        return "negative (piN<piS)"
    return "neutral (~0)"


def main() -> None:
    args = parse_args()
    base = Path(args.base)
    report_dir = Path(args.report_dir) if args.report_dir else base
    kruskal_file = base / "delta_kruskal_by_gene.tsv"
    limma_overall_file = base / "limma_overall_by_gene_all_thresholds.tsv"
    limma_contrasts_file = base / "limma_contrasts_by_gene_all_thresholds.tsv"
    product_summary_file = Path(args.product_summary)

    kr_rows = read_tsv(kruskal_file)
    lim_rows = read_tsv(limma_overall_file)
    con_rows = read_tsv(limma_contrasts_file)
    prod_rows = read_tsv(product_summary_file)

    lim_index = {(r["threshold"], r["product"]): r for r in lim_rows}

    gene_summary = []
    for r in kr_rows:
        key = (r["threshold"], r["product"])
        lim = lim_index.get(key, {})
        gene_summary.append(
            {
                "threshold": r["threshold"],
                "product": r["product"],
                "median_delta_dpi1": r["median_delta_dpi1"],
                "median_delta_dpi3": r["median_delta_dpi3"],
                "median_delta_dpi5": r["median_delta_dpi5"],
                "signal_dpi1": r["signal_dpi1"],
                "signal_dpi3": r["signal_dpi3"],
                "signal_dpi5": r["signal_dpi5"],
                "kruskal_p": r["kruskal_p"],
                "kruskal_bh_fdr": r["bh_fdr_q"],
                "kruskal_significant_fdr": r["significant_fdr"],
                "limma_overall_p": lim.get("p_value", "NA"),
                "limma_overall_bh_fdr": lim.get("bh_fdr", "NA"),
            }
        )

    gene_summary.sort(key=lambda x: (x["threshold"], x["product"]))
    write_tsv(
        report_dir / "selection_gene_summary.tsv",
        gene_summary,
        [
            "threshold",
            "product",
            "median_delta_dpi1",
            "median_delta_dpi3",
            "median_delta_dpi5",
            "signal_dpi1",
            "signal_dpi3",
            "signal_dpi5",
            "kruskal_p",
            "kruskal_bh_fdr",
            "kruskal_significant_fdr",
            "limma_overall_p",
            "limma_overall_bh_fdr",
        ],
    )

    con_rows.sort(key=lambda x: (x["threshold"], x["contrast"], x["product"]))
    write_tsv(
        report_dir / "selection_gene_contrasts.tsv",
        con_rows,
        [
            "threshold",
            "contrast",
            "product",
            "logFC",
            "moderated_t",
            "p_value",
            "bh_fdr",
            "B_stat",
        ],
    )

    selected_threshold = args.threshold

    limma_overall_index = {
        (r["threshold"], r["product"]): r for r in lim_rows
    }
    limma_dpi5_vs_dpi1_index = {
        (r["threshold"], r["product"]): r
        for r in con_rows
        if r.get("contrast") == "dpi5_vs_dpi1"
    }

    pi_by_gene_dpi = defaultdict(lambda: defaultdict(lambda: {"piN": [], "piS": [], "delta": []}))
    gene_codon_length: dict[str, int] = {}
    samples_by_dpi = defaultdict(set)

    for row in prod_rows:
        if row.get("threshold") != selected_threshold:
            continue

        sample = row.get("sample", "")
        if "pilot" in sample.lower():
            continue

        dpi = parse_dpi(sample)
        if dpi not in {"dpi1", "dpi3", "dpi5"}:
            continue

        gene = row["product"]
        pi_n = safe_float(row["piN"])
        pi_s = safe_float(row["piS"])
        if math.isnan(pi_n) or math.isnan(pi_s):
            continue

        pi_by_gene_dpi[gene][dpi]["piN"].append(pi_n)
        pi_by_gene_dpi[gene][dpi]["piS"].append(pi_s)
        pi_by_gene_dpi[gene][dpi]["delta"].append(pi_n - pi_s)
        samples_by_dpi[dpi].add(sample)

        if gene not in gene_codon_length:
            n_sites = safe_float(row.get("N_sites", "NA"))
            s_sites = safe_float(row.get("S_sites", "NA"))
            if not math.isnan(n_sites) and not math.isnan(s_sites):
                gene_codon_length[gene] = int(round((n_sites + s_sites) / 3.0))

    key_rows = []
    for gene in sorted(pi_by_gene_dpi.keys()):
        g = pi_by_gene_dpi[gene]
        mean_piN_dpi1 = mean(g["dpi1"]["piN"])
        mean_piN_dpi3 = mean(g["dpi3"]["piN"])
        mean_piN_dpi5 = mean(g["dpi5"]["piN"])
        mean_piS_dpi1 = mean(g["dpi1"]["piS"])
        mean_piS_dpi3 = mean(g["dpi3"]["piS"])
        mean_piS_dpi5 = mean(g["dpi5"]["piS"])
        mean_delta_dpi1 = mean(g["dpi1"]["delta"])
        mean_delta_dpi3 = mean(g["dpi3"]["delta"])
        mean_delta_dpi5 = mean(g["dpi5"]["delta"])

        overall = limma_overall_index.get((selected_threshold, gene), {})
        dpi5_vs_dpi1 = limma_dpi5_vs_dpi1_index.get((selected_threshold, gene), {})

        key_rows.append(
            {
                "Gene": gene,
                "Length (codons)": str(gene_codon_length.get(gene, "NA")),
                "Mean piN (dpi1)": fmt(mean_piN_dpi1),
                "Mean piN (dpi3)": fmt(mean_piN_dpi3),
                "Mean piN (dpi5)": fmt(mean_piN_dpi5),
                "Mean piS (dpi1)": fmt(mean_piS_dpi1),
                "Mean piS (dpi3)": fmt(mean_piS_dpi3),
                "Mean piS (dpi5)": fmt(mean_piS_dpi5),
                "Mean Delta (dpi1)": fmt(mean_delta_dpi1),
                "Mean Delta (dpi3)": fmt(mean_delta_dpi3),
                "Mean Delta (dpi5)": fmt(mean_delta_dpi5),
                "Limma FDR (overall)": overall.get("bh_fdr", "NA"),
                "Limma FDR (dpi5 vs dpi1)": dpi5_vs_dpi1.get("bh_fdr", "NA"),
                "Direction at dpi5": classify_direction(mean_delta_dpi5, args.direction_epsilon),
            }
        )

    key_table_path = report_dir / "selection_gene_key_table.tsv"
    write_tsv(
        key_table_path,
        key_rows,
        [
            "Gene",
            "Length (codons)",
            "Mean piN (dpi1)",
            "Mean piN (dpi3)",
            "Mean piN (dpi5)",
            "Mean piS (dpi1)",
            "Mean piS (dpi3)",
            "Mean piS (dpi5)",
            "Mean Delta (dpi1)",
            "Mean Delta (dpi3)",
            "Mean Delta (dpi5)",
            "Limma FDR (overall)",
            "Limma FDR (dpi5 vs dpi1)",
            "Direction at dpi5",
        ],
    )

    methods_rows = [
        {"parameter": "analysis_threshold", "value": selected_threshold},
        {"parameter": "timepoints", "value": "dpi1,dpi3,dpi5"},
        {
            "parameter": "samples_included",
            "value": str(sum(len(v) for v in samples_by_dpi.values())),
        },
        {
            "parameter": "samples_per_timepoint",
            "value": ",".join(f"{k}:{len(samples_by_dpi.get(k, []))}" for k in ["dpi1", "dpi3", "dpi5"]),
        },
        {"parameter": "exclude_rule", "value": "sample name contains 'pilot'"},
        {"parameter": "metric_piN", "value": "SNPGenie product_results piN"},
        {"parameter": "metric_piS", "value": "SNPGenie product_results piS"},
        {"parameter": "delta_definition", "value": "Delta = piN - piS"},
        {
            "parameter": "gene_length_codons",
            "value": "round((N_sites + S_sites) / 3)",
        },
        {
            "parameter": "limma_overall",
            "value": "moderated F-test for dpi effect from limma_overall_by_gene_all_thresholds.tsv",
        },
        {
            "parameter": "limma_contrast",
            "value": "dpi5_vs_dpi1 moderated contrast from limma_contrasts_by_gene_all_thresholds.tsv",
        },
        {"parameter": "multiple_testing", "value": "Benjamini-Hochberg FDR (column bh_fdr)"},
        {
            "parameter": "direction_rule_dpi5",
            "value": f"positive if Mean Delta (dpi5) > {args.direction_epsilon}; negative if < -{args.direction_epsilon}; neutral otherwise",
        },
    ]
    methods_path = report_dir / "methods_key_parameters.tsv"
    write_tsv(methods_path, methods_rows, ["parameter", "value"])

    print("Wrote compact tables:")
    print(report_dir / "selection_gene_summary.tsv")
    print(report_dir / "selection_gene_contrasts.tsv")
    print("Wrote key limma table:")
    print(key_table_path)
    print("Wrote methods parameter report:")
    print(methods_path)


if __name__ == "__main__":
    main()
