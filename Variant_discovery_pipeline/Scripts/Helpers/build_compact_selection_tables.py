#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build compact SNPGenie selection tables and a limma-focused key gene report (dataset-aware)."
    )
    parser.add_argument("--base", required=True, help="Directory containing delta/limma outputs")
    parser.add_argument("--product-summary", required=True, help="Path to merged product summary TSV.")
    parser.add_argument("--manifest", required=True, help="Path to samples_manifest.tsv")
    parser.add_argument("--threshold", default="minfreq_0p01", help="Threshold label to use for key table (default: minfreq_0p01)")
    parser.add_argument("--report-dir", default=None, help="Directory to write compact/key output tables")
    parser.add_argument("--direction-epsilon", type=float, default=1e-6, help="Absolute delta cutoff for neutral direction classification.")
    return parser.parse_args()

def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def load_manifest(manifest_path: Path) -> dict[str, str]:
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
    manifest_path = Path(args.manifest)

    kruskal_file = base / "delta_kruskal_by_gene.tsv"
    limma_overall_file = base / "limma_overall_by_gene_all_thresholds.tsv"
    limma_contrasts_file = base / "limma_contrasts_by_gene_all_thresholds.tsv"
    product_summary_file = Path(args.product_summary)

    kr_rows = read_tsv(kruskal_file)
    lim_rows = read_tsv(limma_overall_file)
    con_rows = read_tsv(limma_contrasts_file)
    prod_rows = read_tsv(product_summary_file)

    manifest_map = load_manifest(manifest_path)
    all_dpis = sorted(list(set(manifest_map.values())))
    last_dpi = all_dpis[-1] if all_dpis else "dpi1"

    lim_index = {(r["threshold"], r["product"]): r for r in lim_rows}

    # 1. Selection Gene Summary (merges Kruskal and Limma overall)
    gene_summary = []
    for r in kr_rows:
        key = (r["threshold"], r["product"])
        lim = lim_index.get(key, {})
        row = {
            "threshold": r["threshold"],
            "product": r["product"],
            "kruskal_p": r["kruskal_p"],
            "kruskal_bh_fdr": r["bh_fdr_q"],
            "kruskal_significant_fdr": r["significant_fdr"],
            "limma_overall_p": lim.get("p_value", "NA"),
            "limma_overall_bh_fdr": lim.get("bh_fdr", "NA"),
        }
        for d in all_dpis:
            row[f"median_delta_{d}"] = r.get(f"median_delta_{d}", "NA")
            row[f"signal_{d}"] = r.get(f"signal_{d}", "NA")
        gene_summary.append(row)

    gene_summary.sort(key=lambda x: (x["threshold"], x["product"]))

    summary_fields = ["threshold", "product"]
    for d in all_dpis:
        summary_fields.append(f"median_delta_{d}")
        summary_fields.append(f"signal_{d}")
    summary_fields.extend(["kruskal_p", "kruskal_bh_fdr", "kruskal_significant_fdr", "limma_overall_p", "limma_overall_bh_fdr"])

    write_tsv(report_dir / "selection_gene_summary.tsv", gene_summary, summary_fields)

    # 2. Selection Gene Contrasts
    if con_rows:
        con_rows.sort(key=lambda x: (x["threshold"], x["contrast"], x["product"]))
        write_tsv(
            report_dir / "selection_gene_contrasts.tsv",
            con_rows,
            ["threshold", "contrast", "product", "logFC", "moderated_t", "p_value", "bh_fdr", "B_stat"],
        )
    else:
        # Write empty contrasts file if none exists (single DPI dataset)
        write_tsv(
            report_dir / "selection_gene_contrasts.tsv",
            [],
            ["threshold", "contrast", "product", "logFC", "moderated_t", "p_value", "bh_fdr", "B_stat"],
        )

    # 3. Selection Gene Key Table (specific threshold)
    selected_threshold = args.threshold
    limma_dpi_vs_first_index = {}
    
    # Try to find a contrast comparing the last DPI to the first DPI
    if len(all_dpis) >= 2:
        target_contrast = f"{last_dpi}_vs_{all_dpis[0]}"
        limma_dpi_vs_first_index = {
            (r["threshold"], r["product"]): r
            for r in con_rows
            if r.get("contrast") == target_contrast
        }
    
    pi_by_gene_dpi = defaultdict(lambda: defaultdict(lambda: {"piN": [], "piS": [], "delta": []}))
    gene_codon_length = {}
    samples_by_dpi = defaultdict(set)

    for row in prod_rows:
        if row.get("threshold") != selected_threshold:
            continue

        sample = row.get("sample", "")
        if sample not in manifest_map:
            continue

        dpi = manifest_map[sample]
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
        overall = lim_index.get((selected_threshold, gene), {})
        dpi_vs_first = limma_dpi_vs_first_index.get((selected_threshold, gene), {})

        row_key = {
            "Gene": gene,
            "Length (codons)": str(gene_codon_length.get(gene, "NA")),
            "Limma FDR (overall)": overall.get("bh_fdr", "NA"),
        }
        
        # Add dynamic means
        for d in all_dpis:
            row_key[f"Mean piN ({d})"] = fmt(mean(g[d]["piN"]))
            row_key[f"Mean piS ({d})"] = fmt(mean(g[d]["piS"]))
            row_key[f"Mean Delta ({d})"] = fmt(mean(g[d]["delta"]))

        # Contrast and Direction info
        if len(all_dpis) >= 2:
            target_contrast = f"{last_dpi}_vs_{all_dpis[0]}"
            row_key[f"Limma FDR ({target_contrast})"] = dpi_vs_first.get("bh_fdr", "NA")
        
        mean_delta_last = mean(g[last_dpi]["delta"])
        row_key[f"Direction at {last_dpi}"] = classify_direction(mean_delta_last, args.direction_epsilon)
        
        key_rows.append(row_key)

    key_fields = ["Gene", "Length (codons)"]
    for d in all_dpis:
        key_fields.append(f"Mean piN ({d})")
        key_fields.append(f"Mean piS ({d})")
        key_fields.append(f"Mean Delta ({d})")
    
    key_fields.append("Limma FDR (overall)")
    if len(all_dpis) >= 2:
        target_contrast = f"{last_dpi}_vs_{all_dpis[0]}"
        key_fields.append(f"Limma FDR ({target_contrast})")
    key_fields.append(f"Direction at {last_dpi}")

    key_table_path = report_dir / "selection_gene_key_table.tsv"
    write_tsv(key_table_path, key_rows, key_fields)

    # 4. Methods Parameter Report
    methods_rows = [
        {"parameter": "analysis_threshold", "value": selected_threshold},
        {"parameter": "timepoints", "value": ",".join(all_dpis)},
        {
            "parameter": "samples_included",
            "value": str(sum(len(v) for v in samples_by_dpi.values())),
        },
        {
            "parameter": "samples_per_timepoint",
            "value": ",".join(f"{k}:{len(samples_by_dpi.get(k, []))}" for k in all_dpis),
        },
        {"parameter": "exclude_rule", "value": "sample name not in manifest or no variants"},
        {"parameter": "metric_piN", "value": "SNPGenie product_results piN"},
        {"parameter": "metric_piS", "value": "SNPGenie product_results piS"},
        {"parameter": "delta_definition", "value": "Delta = piN - piS"},
        {"parameter": "gene_length_codons", "value": "round((N_sites + S_sites) / 3)"},
        {"parameter": "multiple_testing", "value": "Benjamini-Hochberg FDR (column bh_fdr)"},
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
