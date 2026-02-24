#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors


T_CRIT_975 = {
    2: 12.706,
    3: 4.303,
    4: 3.182,
    5: 2.776,
    6: 2.571,
    7: 2.447,
    8: 2.365,
    9: 2.306,
    10: 2.262,
    11: 2.228,
    12: 2.201,
    13: 2.179,
    14: 2.160,
    15: 2.145,
    16: 2.131,
    17: 2.120,
    18: 2.110,
    19: 2.101,
    20: 2.093,
    30: 2.045,
}

DESIRED_GENE_ORDER = [
    "nsp1",
    "nsp2",
    "nsp3",
    "nsp4",
    "capsid",
    "e3",
    "e2",
    "6k",
    "e1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create compact publication-ready delta bar plots (threshold x dpi)."
    )
    parser.add_argument(
        "--delta-file",
        default="SNPGenie/analysis/delta_selection/detailed/delta_per_sample.tsv",
        help="Per-sample delta table from analyze_delta_selection.py",
    )
    parser.add_argument(
        "--summary-file",
        default="SNPGenie/analysis/summary/selection_gene_summary.tsv",
        help="Summary table containing limma overall BH-FDR per gene/threshold",
    )
    parser.add_argument(
        "--contrasts-file",
        default="SNPGenie/analysis/summary/selection_gene_contrasts.tsv",
        help="Contrast table containing day-specific limma BH-FDR",
    )
    parser.add_argument(
        "--outdir",
        default="SNPGenie/analysis/summary/figures",
        help="Output directory for PNG/PDF plots",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Top FDR threshold for significance markers (default: 0.05)",
    )
    parser.add_argument(
        "--color-cap",
        type=float,
        default=0.1,
        help="Upper FDR bound for color scale; larger values are shown in gray",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for PNG output",
    )
    parser.add_argument(
        "--yscale",
        choices=["linear", "symlog"],
        default="symlog",
        help="Y-axis scale type (default: symlog)",
    )
    parser.add_argument(
        "--linthresh",
        type=float,
        default=2e-5,
        help="Linear threshold around 0 for symlog scale",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: str) -> float:
    if value in {"", "NA", "*"}:
        return float("nan")
    return float(value)


def mean(values: list[float]) -> float:
    if not values:
        return float("nan")
    return sum(values) / len(values)


def sample_std(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    mu = mean(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (n - 1))


def t_critical_975(n: int) -> float:
    if n <= 1:
        return float("nan")
    if n in T_CRIT_975:
        return T_CRIT_975[n]
    if n < 30:
        nearest = max(k for k in T_CRIT_975 if k < n)
        return T_CRIT_975[nearest]
    return 1.96


def ci95(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    sd = sample_std(values)
    return t_critical_975(n) * sd / math.sqrt(n)


def sem(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    sd = sample_std(values)
    return sd / math.sqrt(n)


def threshold_label(threshold: str) -> str:
    if threshold == "minfreq_0p01":
        return "1%"
    if threshold == "minfreq_0p001":
        return "0.1%"
    return threshold


def ordered_genes_for_plot(available_genes: list[str]) -> list[str]:
    by_lower = {gene.lower(): gene for gene in available_genes}
    ordered = [by_lower[key] for key in DESIRED_GENE_ORDER if key in by_lower]
    remaining = [gene for gene in available_genes if gene not in ordered]
    return ordered + remaining


def fdr_to_stars(fdr: float, alpha: float) -> str:
    if math.isnan(fdr):
        return ""
    if fdr <= 0.001:
        return "***"
    if fdr <= 0.01:
        return "**"
    if fdr <= alpha:
        return "*"
    return ""


def make_plot(
    threshold: str,
    gene_order: list[str],
    dpis: list[str],
    means_by_dpi: dict[str, dict[str, float]],
    ci_by_gene_by_dpi: dict[str, dict[str, float]],
    fdr_by_gene_by_dpi: dict[str, dict[str, float]],
    alpha: float,
    color_cap: float,
    yscale: str,
    linthresh: float,
    outdir: Path,
    dpi: int,
) -> None:
    gene_order = ordered_genes_for_plot(gene_order)
    genes = [
        gene
        for gene in gene_order
        if any(gene in means_by_dpi.get(dpi_group, {}) for dpi_group in dpis)
    ]
    x = list(range(len(genes)))

    cmap = matplotlib.colormaps["viridis"]
    norm = colors.Normalize(vmin=0.0, vmax=color_cap)
    nonsig_color = "#cfcfcf"

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7.2, 4.2),
        sharex=True,
        constrained_layout=True,
    )

    max_abs = 1e-7
    for dpi_group in dpis:
        y_vals = [means_by_dpi.get(dpi_group, {}).get(gene, 0.0) for gene in genes]
        y_errs = [
            0.0
            if math.isnan(ci_by_gene_by_dpi.get(dpi_group, {}).get(gene, float("nan")))
            else ci_by_gene_by_dpi[dpi_group][gene]
            for gene in genes
        ]
        for value, err in zip(y_vals, y_errs):
            max_abs = max(max_abs, abs(value) + err)

    pad = max_abs * 0.15

    for axis_index, dpi_group in enumerate(dpis):
        ax = axes[axis_index]
        y = [means_by_dpi.get(dpi_group, {}).get(gene, 0.0) for gene in genes]
        yerr = [
            0.0
            if math.isnan(ci_by_gene_by_dpi.get(dpi_group, {}).get(gene, float("nan")))
            else ci_by_gene_by_dpi[dpi_group][gene]
            for gene in genes
        ]
        panel_colors = []
        for gene in genes:
            fdr = fdr_by_gene_by_dpi.get(dpi_group, {}).get(gene, float("nan"))
            if math.isnan(fdr) or fdr > color_cap:
                panel_colors.append(nonsig_color)
            else:
                panel_colors.append(cmap(norm(fdr)))

        ax.bar(
            x,
            y,
            yerr=yerr,
            color=panel_colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=2,
        )

        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_ylabel("Mean Δ")
        ax.set_title(f"Threshold {threshold_label(threshold)} | {dpi_group}")
        ax.set_ylim(-max_abs - pad, max_abs + pad)
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=linthresh)

        for idx, gene in enumerate(genes):
            fdr = fdr_by_gene_by_dpi.get(dpi_group, {}).get(gene, float("nan"))
            stars = fdr_to_stars(fdr, alpha)
            if stars:
                offset = max(max_abs * 0.03, linthresh * 0.3)
                ypos = y[idx] + yerr[idx] + offset if y[idx] >= 0 else y[idx] - yerr[idx] - offset
                va = "bottom" if y[idx] >= 0 else "top"
                ax.text(idx, ypos, stars, ha="center", va=va, fontsize=11, fontweight="bold")

        if axis_index == 0:
            ax.tick_params(axis="x", labelbottom=False)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label("FDR")

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(genes, rotation=35, ha="right")
    axes[-1].set_xlabel("Gene")

    stem = outdir / f"delta_bar_stacked_{threshold}"
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_all_dpi_plot(
    threshold: str,
    gene_order: list[str],
    pooled_values_by_gene: dict[str, list[float]],
    overall_fdr_by_gene: dict[str, float],
    filename_tag: str,
    alpha: float,
    color_cap: float,
    yscale: str,
    linthresh: float,
    outdir: Path,
    dpi: int,
) -> None:
    gene_order = ordered_genes_for_plot(gene_order)
    genes = [gene for gene in gene_order if gene in pooled_values_by_gene and pooled_values_by_gene[gene]]
    if not genes:
        return

    cmap = matplotlib.colormaps["viridis"]
    norm = colors.Normalize(vmin=0.0, vmax=color_cap)
    nonsig_color = "#cfcfcf"

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(7.2, 2.8), constrained_layout=True)
    x = list(range(len(genes)))
    y = [mean(pooled_values_by_gene[gene]) for gene in genes]
    yerr = [sem(pooled_values_by_gene[gene]) for gene in genes]
    max_abs = max([abs(v) + (0.0 if math.isnan(e) else e) for v, e in zip(y, yerr)] + [1e-7])

    colors_for_bars = []
    for gene in genes:
        fdr = overall_fdr_by_gene.get(gene, float("nan"))
        if math.isnan(fdr) or fdr > color_cap:
            colors_for_bars.append(nonsig_color)
        else:
            colors_for_bars.append(cmap(norm(fdr)))

    yerr_clean = [0.0 if math.isnan(e) else e for e in yerr]
    ax.bar(
        x,
        y,
        yerr=yerr_clean,
        color=colors_for_bars,
        edgecolor="black",
        linewidth=0.5,
        capsize=2,
    )

    for idx, gene in enumerate(genes):
        fdr = overall_fdr_by_gene.get(gene, float("nan"))
        stars = fdr_to_stars(fdr, alpha)
        if stars:
            offset = max(max_abs * 0.03, linthresh * 0.3)
            ypos = y[idx] + yerr_clean[idx] + offset if y[idx] >= 0 else y[idx] - yerr_clean[idx] - offset
            va = "bottom" if y[idx] >= 0 else "top"
            ax.text(x[idx], ypos, stars, ha="center", va=va, fontsize=10, fontweight="bold")

    pad = max_abs * 0.15
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylim(-max_abs - pad, max_abs + pad)
    if yscale == "symlog":
        ax.set_yscale("symlog", linthresh=linthresh)
    ax.set_title(f"Threshold {threshold_label(threshold)} | all dpi pooled ({filename_tag})")
    ax.set_ylabel("Mean Δ")
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=35, ha="right")
    ax.set_xlabel("Gene")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("FDR")

    stem = outdir / f"delta_bar_all_dpi_{filename_tag}_{threshold}"
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    delta_rows = read_tsv(Path(args.delta_file))
    summary_rows = read_tsv(Path(args.summary_file))
    contrast_rows = read_tsv(Path(args.contrasts_file))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    thresholds = ["minfreq_0p001", "minfreq_0p01"]
    stacked_dpis = ["dpi3", "dpi5"]
    all_dpis = ["dpi1", "dpi3", "dpi5"]

    limma_fdr_by_threshold_gene: dict[tuple[str, str], float] = {}
    kruskal_fdr_by_threshold_gene: dict[tuple[str, str], float] = {}
    gene_order_by_threshold: dict[str, list[str]] = defaultdict(list)
    for row in summary_rows:
        threshold = row["threshold"]
        gene = row["product"]
        limma_fdr = safe_float(row.get("limma_overall_bh_fdr", "NA"))
        kruskal_fdr = safe_float(row.get("kruskal_bh_fdr", "NA"))
        limma_fdr_by_threshold_gene[(threshold, gene)] = limma_fdr
        kruskal_fdr_by_threshold_gene[(threshold, gene)] = kruskal_fdr
        if gene not in gene_order_by_threshold[threshold]:
            gene_order_by_threshold[threshold].append(gene)

    fdr_by_threshold_dpi_gene: dict[tuple[str, str, str], float] = {}
    contrast_to_dpi = {
        "dpi3_vs_dpi1": "dpi3",
        "dpi5_vs_dpi1": "dpi5",
    }
    for row in contrast_rows:
        threshold = row.get("threshold", "")
        contrast = row.get("contrast", "")
        gene = row.get("product", "")
        dpi_group = contrast_to_dpi.get(contrast)
        if not dpi_group:
            continue
        fdr_by_threshold_dpi_gene[(threshold, dpi_group, gene)] = safe_float(row.get("bh_fdr", "NA"))

    replicate_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in delta_rows:
        threshold = row["threshold"]
        dpi_group = row["dpi"]
        gene = row["product"]
        sample = row["sample"]
        if "pilot" in sample.lower():
            continue
        delta_value = safe_float(row["delta_piN_minus_piS"])
        if math.isnan(delta_value):
            continue
        replicate_values[(threshold, dpi_group, gene)].append(delta_value)

    for threshold in thresholds:
        means_by_dpi: dict[str, dict[str, float]] = {}
        ci_by_gene_by_dpi: dict[str, dict[str, float]] = {}
        fdr_by_gene_by_dpi: dict[str, dict[str, float]] = {}
        pooled_values_by_gene: dict[str, list[float]] = {}
        limma_overall_fdr_by_gene: dict[str, float] = {}
        kruskal_overall_fdr_by_gene: dict[str, float] = {}

        genes = gene_order_by_threshold.get(threshold, [])
        for dpi_group in all_dpis:
            means_by_dpi[dpi_group] = {}
            ci_by_gene_by_dpi[dpi_group] = {}
            fdr_by_gene_by_dpi[dpi_group] = {}
            for gene in genes:
                vals = replicate_values.get((threshold, dpi_group, gene), [])
                if not vals:
                    continue
                means_by_dpi[dpi_group][gene] = mean(vals)
                ci_by_gene_by_dpi[dpi_group][gene] = sem(vals)
                if dpi_group in {"dpi3", "dpi5"}:
                    fdr_by_gene_by_dpi[dpi_group][gene] = fdr_by_threshold_dpi_gene.get(
                        (threshold, dpi_group, gene),
                        float("nan"),
                    )
                else:
                    fdr_by_gene_by_dpi[dpi_group][gene] = float("nan")

        for gene in genes:
            pooled_vals = []
            for dpi_group in all_dpis:
                pooled_vals.extend(replicate_values.get((threshold, dpi_group, gene), []))
            if pooled_vals:
                pooled_values_by_gene[gene] = pooled_vals
                limma_overall_fdr_by_gene[gene] = limma_fdr_by_threshold_gene.get((threshold, gene), float("nan"))
                kruskal_overall_fdr_by_gene[gene] = kruskal_fdr_by_threshold_gene.get((threshold, gene), float("nan"))

        has_any = any(means_by_dpi.get(dpi_group) for dpi_group in stacked_dpis)
        if has_any:
            make_plot(
                threshold=threshold,
                gene_order=genes,
                dpis=stacked_dpis,
                means_by_dpi=means_by_dpi,
                ci_by_gene_by_dpi=ci_by_gene_by_dpi,
                fdr_by_gene_by_dpi=fdr_by_gene_by_dpi,
                alpha=args.alpha,
                color_cap=args.color_cap,
                yscale=args.yscale,
                linthresh=args.linthresh,
                outdir=outdir,
                dpi=args.dpi,
            )

        has_any_all = any(means_by_dpi.get(dpi_group) for dpi_group in all_dpis)
        if has_any_all:
            make_all_dpi_plot(
                threshold=threshold,
                gene_order=genes,
                pooled_values_by_gene=pooled_values_by_gene,
                overall_fdr_by_gene=limma_overall_fdr_by_gene,
                filename_tag="limma",
                alpha=args.alpha,
                color_cap=args.color_cap,
                yscale=args.yscale,
                linthresh=args.linthresh,
                outdir=outdir,
                dpi=args.dpi,
            )
            make_all_dpi_plot(
                threshold=threshold,
                gene_order=genes,
                pooled_values_by_gene=pooled_values_by_gene,
                overall_fdr_by_gene=kruskal_overall_fdr_by_gene,
                filename_tag="KW",
                alpha=args.alpha,
                color_cap=args.color_cap,
                yscale=args.yscale,
                linthresh=args.linthresh,
                outdir=outdir,
                dpi=args.dpi,
            )

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
