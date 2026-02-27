#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

VALID_BASES = {"A", "C", "G", "T"}


@dataclass
class SampleSummary:
    sample_id: str
    total_records: int = 0
    pass_records: int = 0
    af_records: int = 0
    snp_records: int = 0
    applied_snps: int = 0
    skipped_non_pass: int = 0
    skipped_low_af: int = 0
    skipped_non_snp: int = 0
    already_alt: int = 0
    skipped_ref_mismatch: int = 0
    skipped_out_of_range: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build PopART inputs from consensus FASTAs by applying LoFreq variants "
            "with AF > 1%% (PASS-only, SNP-only)."
        )
    )
    parser.add_argument(
        "--consensus-dir",
        default="Haplotype/consensus",
        help="Directory containing *.consensus.fasta files.",
    )
    parser.add_argument(
        "--vcf-dir",
        default="SNPGenie/input/vcf_by_sample",
        help="Directory containing per-sample LoFreq VCF files.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/consensus/PopART_lofreq_gt1pct",
        help="Output directory for PopART LoFreq files.",
    )
    parser.add_argument(
        "--min-af",
        type=float,
        default=0.01,
        help="Strictly greater-than AF threshold (default: 0.01).",
    )
    parser.add_argument(
        "--sample-regex",
        default=r"^INH_\d+_DPI_.*$",
        help="Regex to select sample IDs.",
    )
    parser.add_argument(
        "--exclude-dpi-days",
        default="1",
        help="Comma-separated DPI day numbers to exclude (default: 1).",
    )
    return parser.parse_args()


def parse_fasta_single(path: Path) -> tuple[str, str]:
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"Multiple FASTA records found in {path}")
                header = line[1:].strip()
                continue
            chunks.append(line)

    if header is None:
        raise ValueError(f"No FASTA header found in {path}")

    seq = "".join(chunks).upper().replace("U", "T")
    seq = "".join(base if base in VALID_BASES else "N" for base in seq)
    if not seq:
        raise ValueError(f"No sequence found in {path}")
    return header, seq


def format_fasta(records: list[tuple[str, str]], width: int = 80) -> str:
    blocks: list[str] = []
    for seq_id, sequence in records:
        lines = [f">{seq_id}"]
        for idx in range(0, len(sequence), width):
            lines.append(sequence[idx : idx + width])
        blocks.append("\n".join(lines))
    return "\n".join(blocks) + "\n"


def sanitize_trait_value(value: str) -> str:
    token = str(value).strip()
    token = token.replace(",", "_").replace("\t", "_").replace(" ", "_")
    return token if token else "NA"


def write_nexus(
    records: list[tuple[str, str]],
    path: Path,
    trait_labels: list[str] | None = None,
    trait_rows: dict[str, list[str]] | None = None,
) -> None:
    if not records:
        raise ValueError("Cannot write NEXUS: no sequences")

    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError("Cannot write NEXUS: sequence lengths are not aligned")

    nchar = lengths.pop()
    ntax = len(records)

    lines = [
        "#NEXUS",
        "",
        "Begin TAXA;",
        f"  Dimensions ntax={ntax};",
        "  TaxLabels",
    ]
    for seq_id, _ in records:
        lines.append(f"    {seq_id}")
    lines.extend([
        "  ;",
        "End;",
        "",
        "Begin DATA;",
        f"  Dimensions ntax={ntax} nchar={nchar};",
        "  Format datatype=dna missing=? gap=-;",
        "  Matrix",
    ])
    for seq_id, sequence in records:
        lines.append(f"  {seq_id} {sequence}")

    lines.extend(["  ;", "End;", ""])

    if trait_labels and trait_rows:
        lines.extend(
            [
                "Begin TRAITS;",
                f"  Dimensions NTRAITS={len(trait_labels)};",
                "  Format labels=yes missing=? separator=Comma;",
                f"  TraitLabels {' '.join(trait_labels)};",
                "  Matrix",
            ]
        )
        for seq_id, _ in records:
            values = trait_rows.get(seq_id, ["NA"] * len(trait_labels))
            sanitized = [sanitize_trait_value(val) for val in values]
            lines.append(f"  {seq_id} {','.join(sanitized)}")
        lines.extend(["  ;", "End;", ""])

    path.write_text("\n".join(lines), encoding="utf-8")


def parse_info_map(info_field: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for token in info_field.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            parsed[key] = value
    return parsed


def find_vcf(vcf_dir: Path, sample_id: str) -> Path | None:
    plain = vcf_dir / f"{sample_id}.vcf"
    gzipped = vcf_dir / f"{sample_id}.vcf.gz"
    if plain.exists():
        return plain
    if gzipped.exists():
        return gzipped
    return None


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def apply_lofreq_snps(
    sequence: str,
    vcf_path: Path,
    min_af: float,
    summary: SampleSummary,
) -> str:
    edited = list(sequence)
    with open_text(vcf_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            summary.total_records += 1
            pos_text = fields[1]
            ref = fields[3].upper()
            alt = fields[4].upper()
            filt = fields[6]
            info = fields[7]

            if filt != "PASS":
                summary.skipped_non_pass += 1
                continue
            summary.pass_records += 1

            info_map = parse_info_map(info)
            try:
                af = float(info_map.get("AF", "nan"))
            except ValueError:
                af = float("nan")

            if not (af > min_af):
                summary.skipped_low_af += 1
                continue
            summary.af_records += 1

            if "," in alt or len(ref) != 1 or len(alt) != 1:
                summary.skipped_non_snp += 1
                continue
            if ref not in VALID_BASES or alt not in VALID_BASES:
                summary.skipped_non_snp += 1
                continue
            summary.snp_records += 1

            try:
                pos = int(pos_text)
            except ValueError:
                summary.skipped_non_snp += 1
                continue

            idx = pos - 1
            if idx < 0 or idx >= len(edited):
                summary.skipped_out_of_range += 1
                continue

            if edited[idx] == alt:
                summary.already_alt += 1
                continue

            if edited[idx] != ref:
                summary.skipped_ref_mismatch += 1
                continue

            edited[idx] = alt
            summary.applied_snps += 1

    return "".join(edited)


def infer_sample_id(consensus_path: Path) -> str:
    name = consensus_path.name
    if name.endswith(".consensus.fasta"):
        return name[: -len(".consensus.fasta")]
    return consensus_path.stem


def infer_dpi_label(sample_id: str) -> str:
    match = re.match(r"^INH_(\d+)_DPI_", sample_id)
    if match:
        return f"DPI{match.group(1)}"
    return "NA"


def infer_dpi_day(sample_id: str) -> str:
    match = re.match(r"^INH_(\d+)_DPI_", sample_id)
    if match:
        return match.group(1)
    return "NA"


def infer_dpi_one_hot(sample_id: str) -> list[str]:
    day = infer_dpi_day(sample_id)
    return [
        "1" if day == "3" else "0",
        "1" if day == "5" else "0",
    ]


def parse_excluded_days(raw_days: str) -> set[str]:
    return {token.strip() for token in (raw_days or "").split(",") if token.strip()}


def main() -> None:
    args = parse_args()
    consensus_dir = Path(args.consensus_dir)
    vcf_dir = Path(args.vcf_dir)
    out_dir = Path(args.out_dir)
    sample_pattern = re.compile(args.sample_regex)
    excluded_days = parse_excluded_days(args.exclude_dpi_days)

    fasta_dir = out_dir / "fasta"
    nexus_dir = out_dir / "nexus"
    meta_dir = out_dir / "metadata"
    logs_dir = out_dir / "logs"
    for directory in (fasta_dir, nexus_dir, meta_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    consensus_paths = sorted(consensus_dir.glob("*.consensus.fasta"))
    if not consensus_paths:
        raise ValueError(f"No consensus FASTAs found in {consensus_dir}")

    records: list[tuple[str, str]] = []
    summaries: list[SampleSummary] = []
    missing_vcf: list[str] = []

    for consensus_path in consensus_paths:
        sample_id = infer_sample_id(consensus_path)
        if not sample_pattern.match(sample_id):
            continue
        if infer_dpi_day(sample_id) in excluded_days:
            continue

        _, sequence = parse_fasta_single(consensus_path)
        vcf_path = find_vcf(vcf_dir, sample_id)
        if vcf_path is None:
            missing_vcf.append(sample_id)
            continue

        summary = SampleSummary(sample_id=sample_id)
        edited = apply_lofreq_snps(sequence, vcf_path, args.min_af, summary)
        summaries.append(summary)
        records.append((sample_id, edited))

    if not records:
        raise ValueError("No PopART sequences generated. Check filters and inputs.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError("Generated sequences are not all same length")

    fasta_path = fasta_dir / "popart_lofreq_gt1pct.fasta"
    fasta_path.write_text(format_fasta(records), encoding="utf-8")

    nexus_path = nexus_dir / "popart_lofreq_gt1pct.nex"
    trait_labels = ["dpi3", "dpi5"]
    trait_rows = {
        seq_id: infer_dpi_one_hot(seq_id)
        for seq_id, _ in records
    }
    write_nexus(records, nexus_path, trait_labels=trait_labels, trait_rows=trait_rows)

    metadata_path = meta_dir / "sequence_metadata.csv"
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence_id", "sample_id", "dpi", "dpi_day"])
        for seq_id, _ in records:
            writer.writerow([seq_id, seq_id, infer_dpi_label(seq_id), infer_dpi_day(seq_id)])

    traits_path = meta_dir / "sequence_traits.csv"
    with traits_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence_id", "dpi3", "dpi5"])
        for seq_id, _ in records:
            dpi3, dpi5 = infer_dpi_one_hot(seq_id)
            writer.writerow([seq_id, dpi3, dpi5])

    summary_path = logs_dir / "run_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "total_records",
                "pass_records",
                "af_gt_threshold_records",
                "snp_records",
                "applied_snps",
                "already_alt",
                "skipped_non_pass",
                "skipped_low_af",
                "skipped_non_snp",
                "skipped_ref_mismatch",
                "skipped_out_of_range",
            ]
        )
        for row in summaries:
            writer.writerow(
                [
                    row.sample_id,
                    row.total_records,
                    row.pass_records,
                    row.af_records,
                    row.snp_records,
                    row.applied_snps,
                    row.already_alt,
                    row.skipped_non_pass,
                    row.skipped_low_af,
                    row.skipped_non_snp,
                    row.skipped_ref_mismatch,
                    row.skipped_out_of_range,
                ]
            )

    missing_path = logs_dir / "missing_vcf_samples.txt"
    missing_path.write_text("\n".join(missing_vcf) + ("\n" if missing_vcf else ""), encoding="utf-8")

    run_info = logs_dir / "run_info.txt"
    run_info.write_text(
        "\n".join(
            [
                f"timestamp={datetime.now().isoformat(timespec='seconds')}",
                f"consensus_dir={consensus_dir}",
                f"vcf_dir={vcf_dir}",
                f"out_dir={out_dir}",
                f"min_af_strict_gt={args.min_af}",
                f"exclude_dpi_days={','.join(sorted(excluded_days)) if excluded_days else 'none'}",
                f"sequence_count={len(records)}",
                f"missing_vcf_count={len(missing_vcf)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote FASTA: {fasta_path}")
    print(f"Wrote NEXUS: {nexus_path}")
    print(f"Wrote metadata: {metadata_path}")
    print(f"Wrote logs: {logs_dir}")


if __name__ == "__main__":
    main()
