#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a brief VILOCA report with one plot and one table."
    )
    parser.add_argument(
        "--haplotype-summary",
        default="Haplotype/viloca/Analysis/haplotype_level_summary.csv",
        help="Path to haplotype_level_summary.csv",
    )
    parser.add_argument(
        "--long-file",
        default="Haplotype/viloca/Analysis/filtered_linked_mutations_long.csv",
        help="Path to filtered_linked_mutations_long.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/viloca/Analysis",
        help="Output directory",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Top recurrent linked mutations to keep in output table",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def ensure_file(path: Path, label: str) -> None:
    if not path.exists() or not path.is_file():
        fail(f"{label} not found: {path}")


def require_columns(df: pd.DataFrame, required: list[str], label: str) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        fail(f"{label} missing required columns: {', '.join(missing)}")


def summarize_by_dpi(df: pd.DataFrame, value_column: str) -> pd.DataFrame:
    summary = (
        df.groupby("dpi", dropna=False)[value_column]
        .agg(["mean", "min", "max"])
        .reset_index()
    )
    summary["dpi"] = pd.to_numeric(summary["dpi"], errors="coerce")
    summary = summary.dropna(subset=["dpi"]).sort_values("dpi")
    return summary


def draw_panel(
    axis: plt.Axes,
    panel_df: pd.DataFrame,
    title: str,
    y_label: str,
) -> None:
    x = panel_df["dpi"].to_numpy()
    y = panel_df["mean"].to_numpy()
    lower = (panel_df["mean"] - panel_df["min"]).to_numpy()
    upper = (panel_df["max"] - panel_df["mean"]).to_numpy()

    axis.errorbar(
        x,
        y,
        yerr=[lower, upper],
        fmt="o",
        capsize=5,
        linewidth=1.5,
        markersize=6,
        color="#1f77b4",
    )
    axis.set_title(title)
    axis.set_xlabel("DPI")
    axis.set_ylabel(y_label)
    axis.set_xticks(x)
    axis.set_ylim(bottom=0)
    axis.grid(axis="y", alpha=0.3)


def make_summary_figure(haplotype_df: pd.DataFrame, out_path: Path) -> None:
    panel_a = summarize_by_dpi(haplotype_df, "linked_mutation_rows")
    panel_b = summarize_by_dpi(haplotype_df, "unique_mutations")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    draw_panel(
        axes[0],
        panel_a,
        "Panel A: Linked mutation rows by DPI",
        "linked_mutation_rows (mean with min-max)",
    )
    draw_panel(
        axes[1],
        panel_b,
        "Panel B: Unique mutations by DPI",
        "unique_mutations (mean with min-max)",
    )
    fig.suptitle("VILOCA linked mutation summary", fontsize=12)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_recurrent_table(long_df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    grouped = (
        long_df.groupby("mutation_key", dropna=False)
        .agg(
            samples_with_key=("sample_id", "nunique"),
            dpi3_samples=("sample_id", lambda s: s[long_df.loc[s.index, "dpi"] == 3].nunique()),
            dpi5_samples=("sample_id", lambda s: s[long_df.loc[s.index, "dpi"] == 5].nunique()),
            row_matches_in_long=("sample_id", "size"),
        )
        .reset_index()
    )

    if grouped.empty:
        return pd.DataFrame(
            columns=[
                "rank",
                "mutation_key",
                "samples_with_key",
                "prevalence_n_over_8",
                "dpi3_samples",
                "dpi5_samples",
                "row_matches_in_long",
            ]
        )

    grouped["prevalence_n_over_8"] = grouped["samples_with_key"].astype(int).astype(str) + "/8"

    grouped = grouped.sort_values(
        by=["samples_with_key", "row_matches_in_long", "mutation_key"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    grouped.insert(0, "rank", grouped.index + 1)
    grouped = grouped.head(top_n).copy()

    return grouped[
        [
            "rank",
            "mutation_key",
            "samples_with_key",
            "prevalence_n_over_8",
            "dpi3_samples",
            "dpi5_samples",
            "row_matches_in_long",
        ]
    ]


def main() -> None:
    args = parse_args()

    haplotype_summary_path = Path(args.haplotype_summary)
    long_file_path = Path(args.long_file)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ensure_file(haplotype_summary_path, "--haplotype-summary")
    ensure_file(long_file_path, "--long-file")

    haplotype_df = pd.read_csv(haplotype_summary_path)
    long_df = pd.read_csv(long_file_path)

    require_columns(
        haplotype_df,
        ["dpi", "linked_mutation_rows", "unique_mutations"],
        "Haplotype summary",
    )
    require_columns(long_df, ["mutation_key", "sample_id", "dpi"], "Long linked-mutation file")

    figure_path = out_dir / "dpi_linked_mutation_summary.png"
    table_path = out_dir / "top_recurrent_linked_mutations.csv"

    make_summary_figure(haplotype_df, figure_path)
    table_df = make_recurrent_table(long_df, args.top_n)
    table_df.to_csv(table_path, index=False)

    print("Brief report complete")
    print(f"Figure: {figure_path}")
    print(f"Table: {table_path}")
    print(f"Rows written to table: {len(table_df)}")


if __name__ == "__main__":
    main()
