#!/usr/bin/env python3
"""
summarize_linked_mutations.py
Parse VILOCA cooccurring_mutations.csv files, filter linked SNV pairs by depth and support,
and output cleaned long and recurrent linked mutation tables.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

__version__ = "1.0.0"

REQUIRED_COLUMNS = ["haplotype_id", "chrom", "position", "ref", "var", "reads", "support", "coverage"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize linked mutations from VILOCA cooccurrence outputs."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--viloca-input",
        nargs="*",
        help="Path(s) to VILOCA cooccurring_mutations.csv files or sample directories.",
    )
    parser.add_argument(
        "--samplesheet",
        help="Path to samplesheet CSV (sample,fastq_1,fastq_2,treatment) for metadata.",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for summary CSV files.",
    )
    parser.add_argument(
        "--min-reads",
        type=float,
        default=10.0,
        help="Minimum supporting read count for linked mutations (default: 10.0).",
    )
    parser.add_argument(
        "--min-support",
        type=float,
        default=0.80,
        help="Minimum posterior support for VILOCA haplotype window (default: 0.80).",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=100.0,
        help="Minimum local coverage depth in window (default: 100.0).",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="Minimum number of samples for recurrent linkage (default: 2).",
    )
    return parser.parse_args()


def load_samplesheet_metadata(samplesheet_path: Path | None) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not samplesheet_path or not samplesheet_path.exists():
        return mapping
    with samplesheet_path.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample_id = (row.get("sample") or row.get("sample_id") or "").strip()
            treatment = (row.get("treatment") or "unknown").strip()
            if sample_id:
                mapping[sample_id] = treatment
    return mapping


def find_cooccurrence_files(inputs: List[str] | None) -> List[Path]:
    files: List[Path] = []
    if not inputs:
        return files
    for item in inputs:
        p = Path(item)
        if p.is_file() and (p.name.endswith(".csv") or "cooccurring_mutations" in p.name):
            files.append(p)
        elif p.is_dir():
            files.extend(sorted(p.glob("**/*cooccurring_mutations*.csv")))
            files.extend(sorted(p.glob("**/*.csv")))
    return sorted(list(set(files)))


def parse_cooccurrence_csv(fpath: Path, sample_id: str, treatment: str, thresholds: Tuple[float, float, float]) -> List[dict]:
    min_reads, min_support, min_cov = thresholds
    rows: List[dict] = []

    if not fpath.exists() or fpath.stat().st_size == 0:
        return rows

    with fpath.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return rows

        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            # Skip non-matching CSVs silently
            return rows

        for row in reader:
            hap_id = (row.get("haplotype_id") or "").strip()
            if "reference" in hap_id.lower():
                continue

            try:
                chrom = (row.get("chrom") or "").strip()
                pos = int(float(row.get("position") or 0))
                ref = (row.get("ref") or "").strip().upper()
                var = (row.get("var") or "").strip().upper()
                reads = float(row.get("reads") or 0.0)
                support = float(row.get("support") or 0.0)
                coverage = float(row.get("coverage") or 0.0)
            except (ValueError, TypeError):
                continue

            if not ref or not var or ref == var:
                continue

            if reads < min_reads or support < min_support or coverage < min_cov:
                continue

            reads_over_cov = round(reads / coverage, 6) if coverage > 0 else 0.0

            rows.append({
                "sample_id": sample_id,
                "treatment": treatment,
                "haplotype_id": hap_id,
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "var": var,
                "reads": round(reads, 4),
                "support": round(support, 6),
                "coverage": round(coverage, 4),
                "reads_over_coverage": reads_over_cov,
                "mutation_key": f"{pos}{ref}>{var}",
            })

    return rows


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    long_out_path = out_dir / "linked_mutations_long.csv"
    recurrent_out_path = out_dir / "linked_mutations_recurrent.csv"

    metadata = load_samplesheet_metadata(Path(args.samplesheet) if args.samplesheet else None)
    files = find_cooccurrence_files(args.viloca_input)

    thresholds = (args.min_reads, args.min_support, args.min_coverage)

    long_rows: List[dict] = []
    for fpath in files:
        # Determine sample ID from filename if possible
        sample_id = fpath.stem.replace("_cooccurring_mutations", "").replace(".cooccurring_mutations", "")
        if sample_id in metadata:
            treatment = metadata[sample_id]
        else:
            # Fallback: inspect parent dir
            sample_id = fpath.parent.name if fpath.parent.name != "." else sample_id
            treatment = metadata.get(sample_id, "unknown")

        if treatment.strip().lower() == "mock":
            continue

        parsed = parse_cooccurrence_csv(fpath, sample_id, treatment, thresholds)
        long_rows.extend(parsed)

    # Empty handling
    if not long_rows:
        sys.stdout.write("Notice: No linked mutations passed thresholds. Writing empty output CSVs.\n")
        with long_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "treatment", "haplotype_id", "chrom", "position", "ref", "var", "reads", "support", "coverage", "reads_over_coverage"])
        with recurrent_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["mutation_pair", "n_samples", "sample_ids", "treatments", "mean_support", "mean_reads"])
        sys.exit(0)

    # Write linked_mutations_long.csv
    fieldnames_long = ["sample_id", "treatment", "haplotype_id", "chrom", "position", "ref", "var", "reads", "support", "coverage", "reads_over_coverage"]
    with long_out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_long, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(long_rows)

    # Identify linked mutation pairs per haplotype window
    # Group rows by (sample_id, haplotype_id)
    hap_groups: Dict[Tuple[str, str], List[dict]] = {}
    for r in long_rows:
        hap_groups.setdefault((r["sample_id"], r["haplotype_id"]), []).append(r)

    # Extract co-occurring pairs per window
    pair_samples: Dict[str, Dict[str, dict]] = {}  # pair_key -> sample_id -> {treatment, support_list, reads_list}

    for (sample_id, hap_id), muts in hap_groups.items():
        if len(muts) < 2:
            continue
        # Sort muts by position
        sorted_muts = sorted(muts, key=lambda x: x["position"])
        for i in range(len(sorted_muts)):
            for j in range(i + 1, len(sorted_muts)):
                m1 = sorted_muts[i]
                m2 = sorted_muts[j]
                pair_key = f"{m1['mutation_key']} + {m2['mutation_key']}"
                treatment = m1["treatment"]

                if pair_key not in pair_samples:
                    pair_samples[pair_key] = {}
                if sample_id not in pair_samples[pair_key]:
                    pair_samples[pair_key][sample_id] = {
                        "treatment": treatment,
                        "supports": [],
                        "reads": [],
                    }
                pair_samples[pair_key][sample_id]["supports"].append((m1["support"] + m2["support"]) / 2.0)
                pair_samples[pair_key][sample_id]["reads"].append((m1["reads"] + m2["reads"]) / 2.0)

    # Write linked_mutations_recurrent.csv
    recurrent_rows = []
    for pair_key, sample_map in pair_samples.items():
        n_samples = len(sample_map)
        if n_samples < args.min_samples:
            continue

        sample_ids = sorted(list(sample_map.keys()))
        treatments = sorted(list(set(info["treatment"] for info in sample_map.values())))

        all_supports = [s for info in sample_map.values() for s in info["supports"]]
        all_reads = [r for info in sample_map.values() for r in info["reads"]]

        mean_support = round(sum(all_supports) / max(1, len(all_supports)), 6)
        mean_reads = round(sum(all_reads) / max(1, len(all_reads)), 4)

        recurrent_rows.append({
            "mutation_pair": pair_key,
            "n_samples": n_samples,
            "sample_ids": ";".join(sample_ids),
            "treatments": ";".join(treatments),
            "mean_support": mean_support,
            "mean_reads": mean_reads,
        })

    recurrent_rows.sort(key=lambda x: (-x["n_samples"], -x["mean_support"], x["mutation_pair"]))

    fieldnames_rec = ["mutation_pair", "n_samples", "sample_ids", "treatments", "mean_support", "mean_reads"]
    with recurrent_out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_rec)
        writer.writeheader()
        writer.writerows(recurrent_rows)

    sys.stdout.write(f"Done. Processed {len(long_rows)} linked mutation records into {len(recurrent_rows)} recurrent pairs.\n")


if __name__ == "__main__":
    main()
