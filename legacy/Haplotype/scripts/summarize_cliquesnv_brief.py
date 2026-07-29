#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


SAMPLE_PATTERN = re.compile(r"^INH_(?P<dpi>\d+)_DPI_(?P<replicate>R\d+)_.*$")
TF_DIR_PATTERN = re.compile(r"^tf_(?P<tf>\d+p\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a brief per-sample summary table from CliqueSNV outputs."
    )
    parser.add_argument(
        "--cliquesnv-dir",
        default="Haplotype/cliquesnv",
        help="Directory containing CliqueSNV sample subdirectories.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/cliquesnv/Analysis",
        help="Directory where summary TSV is written.",
    )
    return parser.parse_args()


def parse_tf(tf_dir_name: str) -> float | None:
    match = TF_DIR_PATTERN.match(tf_dir_name)
    if not match:
        return None
    value = match.group("tf").replace("p", ".")
    try:
        return float(value)
    except ValueError:
        return None


def parse_sample(sample_id: str) -> tuple[int | None, str | None]:
    match = SAMPLE_PATTERN.match(sample_id)
    if not match:
        return None, None
    return int(match.group("dpi")), match.group("replicate")


def pick_single_file(folder: Path, suffix: str) -> Path | None:
    candidates = sorted(folder.glob(f"*{suffix}"))
    if not candidates:
        return None
    return candidates[-1]


def count_fasta_records(fasta_path: Path | None) -> int | None:
    if fasta_path is None or not fasta_path.exists():
        return None
    count = 0
    with fasta_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def summarize_json(json_path: Path | None) -> dict[str, object]:
    result: dict[str, object] = {
        "json_error": None,
        "found_haplotypes_json": None,
        "top_haplotype_frequency": None,
        "frequency_sum": None,
        "haplotype_seq_len_min_bp": None,
        "haplotype_seq_len_max_bp": None,
        "covered_start_min_1based": None,
        "covered_end_max_1based": None,
        "covered_span_min_bp": None,
        "covered_span_max_bp": None,
        "snp_pos_min": None,
        "snp_pos_max": None,
        "snp_pos_span_bp": None,
    }

    if json_path is None or not json_path.exists():
        result["json_error"] = "missing_json"
        return result

    try:
        payload = json.loads(json_path.read_text(encoding="utf-8", errors="replace"))
    except Exception as exc:
        result["json_error"] = f"json_parse_error:{exc}"
        return result

    error_field = payload.get("error", "none")
    result["json_error"] = error_field

    found_haplotypes = payload.get("foundHaplotypes")
    if isinstance(found_haplotypes, int):
        result["found_haplotypes_json"] = found_haplotypes

    haplotypes = payload.get("haplotypes")
    frequencies: list[float] = []
    seq_lengths: list[int] = []
    covered_starts: list[int] = []
    covered_ends: list[int] = []
    covered_spans: list[int] = []
    all_snp_positions: list[int] = []

    def parse_snp_positions(snps_value: object) -> list[int]:
        if not isinstance(snps_value, str):
            return []
        match = re.search(r"\[(.*)\]", snps_value)
        if not match:
            return []
        raw = match.group(1).strip()
        if raw == "":
            return []
        values: list[int] = []
        for token in raw.split(","):
            token = token.strip()
            if token == "":
                continue
            try:
                values.append(int(token))
            except ValueError:
                continue
        return values

    def non_n_span(sequence: str) -> tuple[int, int] | None:
        first = None
        last = None
        for idx, base in enumerate(sequence, start=1):
            if base.upper() != "N":
                if first is None:
                    first = idx
                last = idx
        if first is None or last is None:
            return None
        return first, last

    if isinstance(haplotypes, list):
        for haplotype in haplotypes:
            if not isinstance(haplotype, dict):
                continue
            frequency = haplotype.get("frequency")
            if isinstance(frequency, (int, float)):
                frequencies.append(float(frequency))

            sequence = haplotype.get("haplotype")
            if isinstance(sequence, str) and sequence != "":
                seq_lengths.append(len(sequence))
                span = non_n_span(sequence)
                if span is not None:
                    start_1b, end_1b = span
                    covered_starts.append(start_1b)
                    covered_ends.append(end_1b)
                    covered_spans.append(end_1b - start_1b + 1)

            all_snp_positions.extend(parse_snp_positions(haplotype.get("snps")))

    if frequencies:
        result["top_haplotype_frequency"] = max(frequencies)
        result["frequency_sum"] = sum(frequencies)

    if seq_lengths:
        result["haplotype_seq_len_min_bp"] = min(seq_lengths)
        result["haplotype_seq_len_max_bp"] = max(seq_lengths)

    if covered_starts and covered_ends and covered_spans:
        result["covered_start_min_1based"] = min(covered_starts)
        result["covered_end_max_1based"] = max(covered_ends)
        result["covered_span_min_bp"] = min(covered_spans)
        result["covered_span_max_bp"] = max(covered_spans)

    if all_snp_positions:
        snp_min = min(all_snp_positions)
        snp_max = max(all_snp_positions)
        result["snp_pos_min"] = snp_min
        result["snp_pos_max"] = snp_max
        result["snp_pos_span_bp"] = snp_max - snp_min + 1

    return result


def status_from_error(json_error: object) -> str:
    if json_error is None:
        return "unknown"
    if isinstance(json_error, str) and json_error in {"none", "null", ""}:
        return "ok"
    if isinstance(json_error, str) and json_error.startswith("missing_json"):
        return "missing"
    if isinstance(json_error, str) and json_error.startswith("json_parse_error"):
        return "parse_error"
    return "tool_error"


def main() -> None:
    args = parse_args()
    cliquesnv_dir = Path(args.cliquesnv_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    if not cliquesnv_dir.exists() or not cliquesnv_dir.is_dir():
        output_path = out_dir / "cliquesnv_brief_per_sample.tsv"
        pd.DataFrame().to_csv(output_path, sep="\t", index=False)
        print(f"Wrote empty summary (missing directory): {output_path}")
        return

    for sample_dir in sorted(cliquesnv_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        if sample_dir.name in {"logs", "Analysis"}:
            continue

        dpi, replicate_id = parse_sample(sample_dir.name)

        for tf_dir in sorted(sample_dir.iterdir()):
            if not tf_dir.is_dir():
                continue

            tf = parse_tf(tf_dir.name)
            if tf is None:
                continue

            json_file = pick_single_file(tf_dir, ".json")
            fasta_file = pick_single_file(tf_dir, ".fasta")

            json_summary = summarize_json(json_file)
            found_haplotypes_fasta = count_fasta_records(fasta_file)
            found_haplotypes_json = json_summary.get("found_haplotypes_json")

            found_haplotypes = found_haplotypes_json
            if found_haplotypes is None:
                found_haplotypes = found_haplotypes_fasta

            top_frequency = json_summary.get("top_haplotype_frequency")
            frequency_sum = json_summary.get("frequency_sum")

            if isinstance(top_frequency, float):
                top_frequency = round(top_frequency, 8)
            if isinstance(frequency_sum, float):
                frequency_sum = round(frequency_sum, 8)

            json_error = json_summary.get("json_error")

            rows.append(
                {
                    "sample_id": sample_dir.name,
                    "dpi": dpi,
                    "replicate_id": replicate_id,
                    "threshold": tf,
                    "threshold_label": tf_dir.name,
                    "status": status_from_error(json_error),
                    "json_error": json_error,
                    "found_haplotypes": found_haplotypes,
                    "top_haplotype_frequency": top_frequency,
                    "frequency_sum": frequency_sum,
                    "haplotype_seq_len_min_bp": json_summary.get("haplotype_seq_len_min_bp"),
                    "haplotype_seq_len_max_bp": json_summary.get("haplotype_seq_len_max_bp"),
                    "covered_start_min_1based": json_summary.get("covered_start_min_1based"),
                    "covered_end_max_1based": json_summary.get("covered_end_max_1based"),
                    "covered_span_min_bp": json_summary.get("covered_span_min_bp"),
                    "covered_span_max_bp": json_summary.get("covered_span_max_bp"),
                    "snp_pos_min": json_summary.get("snp_pos_min"),
                    "snp_pos_max": json_summary.get("snp_pos_max"),
                    "snp_pos_span_bp": json_summary.get("snp_pos_span_bp"),
                    "json_file": str(json_file) if json_file else None,
                    "fasta_file": str(fasta_file) if fasta_file else None,
                    "run_dir": str(tf_dir),
                }
            )

    output_path = out_dir / "cliquesnv_brief_per_sample.tsv"

    if not rows:
        pd.DataFrame(
            columns=[
                "sample_id",
                "dpi",
                "replicate_id",
                "threshold",
                "threshold_label",
                "status",
                "json_error",
                "found_haplotypes",
                "top_haplotype_frequency",
                "frequency_sum",
                "haplotype_seq_len_min_bp",
                "haplotype_seq_len_max_bp",
                "covered_start_min_1based",
                "covered_end_max_1based",
                "covered_span_min_bp",
                "covered_span_max_bp",
                "snp_pos_min",
                "snp_pos_max",
                "snp_pos_span_bp",
                "json_file",
                "fasta_file",
                "run_dir",
            ]
        ).to_csv(output_path, sep="\t", index=False)
        print(f"Wrote empty summary: {output_path}")
        return

    dataframe = pd.DataFrame(rows).sort_values(
        by=["threshold", "dpi", "sample_id", "replicate_id"],
        na_position="last",
    )
    dataframe.to_csv(output_path, sep="\t", index=False)

    print(f"Rows written: {len(dataframe)}")
    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
