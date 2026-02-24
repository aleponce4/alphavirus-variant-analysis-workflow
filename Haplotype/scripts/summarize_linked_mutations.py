#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "haplotype_id",
    "chrom",
    "start",
    "position",
    "ref",
    "var",
    "reads",
    "support",
    "coverage",
]

COOCCURRENCE_CANDIDATES = [
    "cooccurring_mutations.csv",
    "work/snv/cooccurring_mutations.csv",
    "snv/cooccurring_mutations.csv",
]


@dataclass(frozen=True)
class SampleInfo:
    sample_id: str
    dpi: int
    replicate_id: str
    sample_dir: Path


@dataclass
class RunStats:
    samples_found: int = 0
    samples_used: int = 0
    samples_skipped_no_file: int = 0
    samples_skipped_parse_error: int = 0
    rows_loaded: int = 0
    rows_after_cleaning: int = 0
    rows_after_thresholds: int = 0
    rows_after_haplotype_size_filter: int = 0


@dataclass(frozen=True)
class FilterThresholds:
    min_reads: float
    min_support: float
    min_coverage: float
    min_reads_ratio: float
    min_mutations_per_haplotype: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize linked mutations from VILOCA cooccurrence files."
    )
    parser.add_argument(
        "--viloca-dir",
        default="Haplotype/viloca",
        help="Directory containing VILOCA sample subdirectories.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/viloca/Analysis",
        help="Output directory for summary CSVs.",
    )
    parser.add_argument(
        "--sample-pattern",
        default=r"^INH_(?P<dpi>\d+)_DPI_(?P<replicate>R\d+)_.*$",
        help=(
            "Regex used to identify sample directories. Must include named groups "
            "'dpi' and 'replicate'."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        nargs="+",
        default=[3, 5],
        help="DPI values to include (default: 3 5).",
    )

    parser.add_argument("--min-reads", type=float, default=10.0)
    parser.add_argument("--min-support", type=float, default=0.80)
    parser.add_argument("--min-coverage", type=float, default=1000.0)
    parser.add_argument("--min-reads-ratio", type=float, default=0.01)
    parser.add_argument("--min-mutations-per-haplotype", type=int, default=2)

    return parser.parse_args()


def discover_samples(
    viloca_dir: Path,
    sample_regex: re.Pattern[str],
    allowed_dpi: set[int],
) -> list[SampleInfo]:
    samples: list[SampleInfo] = []
    for child in sorted(viloca_dir.iterdir()):
        if not child.is_dir():
            continue
        match = sample_regex.match(child.name)
        if not match:
            continue

        dpi = int(match.group("dpi"))
        if dpi not in allowed_dpi:
            continue

        replicate_id = match.group("replicate")
        samples.append(
            SampleInfo(
                sample_id=child.name,
                dpi=dpi,
                replicate_id=replicate_id,
                sample_dir=child,
            )
        )
    return samples


def pick_cooccurrence_file(sample_dir: Path) -> Path | None:
    for relative_path in COOCCURRENCE_CANDIDATES:
        candidate = sample_dir / relative_path
        if candidate.exists() and candidate.is_file():
            return candidate
    return None


def read_cooccurrence_csv(csv_path: Path) -> pd.DataFrame:
    dataframe = pd.read_csv(csv_path)

    unnamed_cols = [column for column in dataframe.columns if str(column).startswith("Unnamed")]
    if unnamed_cols:
        dataframe = dataframe.drop(columns=unnamed_cols)

    missing = [column for column in REQUIRED_COLUMNS if column not in dataframe.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    dataframe = dataframe.loc[:, REQUIRED_COLUMNS].copy()

    for column in ["position", "start", "reads", "support", "coverage"]:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")

    dataframe["position"] = to_integer_when_possible(dataframe["position"])
    dataframe["start"] = to_integer_when_possible(dataframe["start"])
    dataframe["reads"] = dataframe["reads"].astype("Float64")
    dataframe["support"] = dataframe["support"].astype("Float64")
    dataframe["coverage"] = dataframe["coverage"].astype("Float64")

    dataframe["haplotype_id"] = dataframe["haplotype_id"].astype("string").str.strip()
    dataframe["chrom"] = dataframe["chrom"].astype("string").str.strip()
    dataframe["ref"] = dataframe["ref"].astype("string").str.strip().str.upper()
    dataframe["var"] = dataframe["var"].astype("string").str.strip().str.upper()

    return dataframe


def to_integer_when_possible(series: pd.Series) -> pd.Series:
    rounded = series.round(0)
    integral_mask = series.notna() & (series == rounded)
    integer_series = rounded.where(integral_mask)
    return integer_series.astype("Int64")


def parse_window_bounds(haplotype_id: str) -> tuple[int | None, int | None]:
    match = re.search(r"-(\d+)-(\d+)$", haplotype_id)
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def add_metadata(df: pd.DataFrame, sample: SampleInfo, source_file: Path) -> pd.DataFrame:
    windows = df["haplotype_id"].apply(parse_window_bounds)
    df["sample_id"] = sample.sample_id
    df["dpi"] = sample.dpi
    df["replicate_id"] = sample.replicate_id
    df["source_file"] = str(source_file)
    df["window_start"] = windows.apply(lambda pair: pair[0]).astype("Int64")
    df["window_end"] = windows.apply(lambda pair: pair[1]).astype("Int64")

    df["haplotype_window_key"] = (
        df["sample_id"].astype("string") + "|" + df["haplotype_id"].astype("string")
    )
    df["mutation_key"] = (
        df["chrom"].astype("string")
        + ":"
        + df["position"].astype("string")
        + ":"
        + df["ref"].astype("string")
        + ">"
        + df["var"].astype("string")
    )
    df["reads_over_coverage"] = (df["reads"] / df["coverage"]).astype("Float64")
    return df


def remove_reference_or_non_mutation_rows(df: pd.DataFrame) -> pd.DataFrame:
    non_reference = ~df["haplotype_id"].str.contains("reference", case=False, na=False)
    non_blank_pos = df["position"].notna()
    non_blank_ref = df["ref"].notna() & (df["ref"] != "")
    non_blank_var = df["var"].notna() & (df["var"] != "")
    keep = non_reference & non_blank_pos & non_blank_ref & non_blank_var
    return df.loc[keep].copy()


def apply_filters(df: pd.DataFrame, thresholds: FilterThresholds) -> pd.DataFrame:
    filtered = df.loc[
        (df["reads"] >= thresholds.min_reads)
        & (df["support"] >= thresholds.min_support)
        & (df["coverage"] >= thresholds.min_coverage)
        & (df["reads_over_coverage"] >= thresholds.min_reads_ratio)
    ].copy()

    haplotype_counts = filtered.groupby("haplotype_window_key").size()
    keep_haplotypes = haplotype_counts[
        haplotype_counts >= thresholds.min_mutations_per_haplotype
    ].index

    return filtered.loc[filtered["haplotype_window_key"].isin(keep_haplotypes)].copy()


def make_collapsed_output(long_df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        long_df.sort_values(["sample_id", "haplotype_window_key", "position", "mutation_key"])
        .groupby(
            [
                "sample_id",
                "dpi",
                "replicate_id",
                "haplotype_id",
                "window_start",
                "window_end",
                "haplotype_window_key",
            ],
            dropna=False,
        )
        .agg(
            mutation_count=("mutation_key", "size"),
            mutation_keys=("mutation_key", lambda values: ";".join(values)),
            mean_reads=("reads", "mean"),
            mean_support=("support", "mean"),
            mean_coverage=("coverage", "mean"),
            mean_reads_over_coverage=("reads_over_coverage", "mean"),
        )
        .reset_index()
    )

    grouped["mean_reads"] = grouped["mean_reads"].round(4)
    grouped["mean_support"] = grouped["mean_support"].round(6)
    grouped["mean_coverage"] = grouped["mean_coverage"].round(4)
    grouped["mean_reads_over_coverage"] = grouped["mean_reads_over_coverage"].round(6)
    return grouped


def make_haplotype_summary(long_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        long_df.groupby(["sample_id", "dpi", "replicate_id"], dropna=False)
        .agg(
            haplotypes_passing=("haplotype_window_key", "nunique"),
            linked_mutation_rows=("mutation_key", "size"),
            unique_mutations=("mutation_key", "nunique"),
            median_support=("support", "median"),
            median_reads=("reads", "median"),
            median_coverage=("coverage", "median"),
        )
        .reset_index()
        .sort_values(["dpi", "sample_id"])
    )

    summary["median_support"] = summary["median_support"].round(6)
    summary["median_reads"] = summary["median_reads"].round(4)
    summary["median_coverage"] = summary["median_coverage"].round(4)
    return summary


def main() -> None:
    args = parse_args()
    viloca_dir = Path(args.viloca_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_regex = re.compile(args.sample_pattern)
    thresholds = FilterThresholds(
        min_reads=args.min_reads,
        min_support=args.min_support,
        min_coverage=args.min_coverage,
        min_reads_ratio=args.min_reads_ratio,
        min_mutations_per_haplotype=args.min_mutations_per_haplotype,
    )

    stats = RunStats()
    samples = discover_samples(viloca_dir, sample_regex, set(args.dpi))
    stats.samples_found = len(samples)

    all_frames: list[pd.DataFrame] = []

    for sample in samples:
        cooccurrence_file = pick_cooccurrence_file(sample.sample_dir)
        if cooccurrence_file is None:
            stats.samples_skipped_no_file += 1
            continue

        try:
            frame = read_cooccurrence_csv(cooccurrence_file)
        except Exception:
            stats.samples_skipped_parse_error += 1
            continue

        stats.samples_used += 1
        stats.rows_loaded += len(frame)

        frame = add_metadata(frame, sample, cooccurrence_file)
        frame = remove_reference_or_non_mutation_rows(frame)
        stats.rows_after_cleaning += len(frame)
        all_frames.append(frame)

    long_output_path = out_dir / "filtered_linked_mutations_long.csv"
    collapsed_output_path = out_dir / "filtered_linked_mutations_collapsed.csv"
    haplotype_summary_path = out_dir / "haplotype_level_summary.csv"

    if not all_frames:
        empty_long = pd.DataFrame(
            columns=REQUIRED_COLUMNS
            + [
                "sample_id",
                "dpi",
                "replicate_id",
                "source_file",
                "window_start",
                "window_end",
                "haplotype_window_key",
                "mutation_key",
                "reads_over_coverage",
            ]
        )
        empty_long.to_csv(long_output_path, index=False)
        pd.DataFrame().to_csv(collapsed_output_path, index=False)
        pd.DataFrame().to_csv(haplotype_summary_path, index=False)

        print(f"Samples found: {stats.samples_found}")
        print("Samples used: 0")
        print(
            "Samples skipped: "
            f"no_cooccurrence_file={stats.samples_skipped_no_file}, "
            f"parse_error={stats.samples_skipped_parse_error}"
        )
        print("Rows: loaded=0, after_cleaning=0, after_thresholds=0, after_haplotype_filter=0")
        print(f"Wrote: {long_output_path}")
        print(f"Wrote: {collapsed_output_path}")
        print(f"Wrote: {haplotype_summary_path}")
        return

    long_df = pd.concat(all_frames, ignore_index=True)

    threshold_df = long_df.loc[
        (long_df["reads"] >= thresholds.min_reads)
        & (long_df["support"] >= thresholds.min_support)
        & (long_df["coverage"] >= thresholds.min_coverage)
        & (long_df["reads_over_coverage"] >= thresholds.min_reads_ratio)
    ].copy()
    stats.rows_after_thresholds = len(threshold_df)

    filtered_df = apply_filters(long_df, thresholds)
    stats.rows_after_haplotype_size_filter = len(filtered_df)

    filtered_df = filtered_df.sort_values(
        ["sample_id", "window_start", "window_end", "haplotype_id", "position"]
    )
    filtered_df.to_csv(long_output_path, index=False)

    collapsed_df = make_collapsed_output(filtered_df)
    collapsed_df.to_csv(collapsed_output_path, index=False)

    haplotype_summary_df = make_haplotype_summary(filtered_df)
    haplotype_summary_df.to_csv(haplotype_summary_path, index=False)

    print(f"Samples found: {stats.samples_found}")
    print(f"Samples used: {stats.samples_used}")
    print(
        "Samples skipped: "
        f"no_cooccurrence_file={stats.samples_skipped_no_file}, "
        f"parse_error={stats.samples_skipped_parse_error}"
    )
    print(
        "Rows: "
        f"loaded={stats.rows_loaded}, "
        f"after_cleaning={stats.rows_after_cleaning}, "
        f"after_thresholds={stats.rows_after_thresholds}, "
        f"after_haplotype_filter={stats.rows_after_haplotype_size_filter}"
    )
    print(f"Wrote: {long_output_path}")
    print(f"Wrote: {collapsed_output_path}")
    print(f"Wrote: {haplotype_summary_path}")


if __name__ == "__main__":
    main()
