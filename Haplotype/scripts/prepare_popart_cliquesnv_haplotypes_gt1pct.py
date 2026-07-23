#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

SNP_TOKEN_PATTERN = re.compile(r"^(?P<ref>[ACGT])(?P<pos>\d+)(?P<alt>[ACGT])$")
VALID_BASES = {"A", "C", "G", "T"}


@dataclass
class HaplotypeSummary:
    sample_id: str
    haplotype_id: str
    frequency_percent: float
    mutation_tokens_total: int
    mutation_tokens_applied: int
    skipped_non_snp: int
    skipped_ref_mismatch: int
    skipped_out_of_range: int


@dataclass
class HaplotypeRow:
    sample_id: str
    haplotype_id: str
    frequency_percent: float
    edited_sequence: str
    mutation_tokens_total: int
    mutation_tokens_applied: int
    skipped_non_snp: int
    skipped_ref_mismatch: int
    skipped_out_of_range: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build PopART inputs from CliqueSNV haplotypes by applying SNP mutation tokens "
            "to per-sample consensus sequences."
        )
    )
    parser.add_argument(
        "--haplotype-csv",
        default="Haplotype/cliquesnv/Analysis/cliquesnv_haplotype_frequency_by_sample.csv",
        help="CSV table containing per-sample haplotype frequencies and mutation lists.",
    )
    parser.add_argument(
        "--consensus-dir",
        default="Haplotype/consensus",
        help="Directory containing *.consensus.fasta files.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/consensus/PopART_cliquesnv_haplotypes_gt1pct",
        help="Output directory for PopART CliqueSNV haplotype files.",
    )
    parser.add_argument(
        "--min-frequency-percent",
        type=float,
        default=1.0,
        help="Strictly greater-than frequency cutoff in percent (default: 1.0).",
    )
    parser.add_argument(
        "--sample-regex",
        default=r"^INH_\d+_DPI_.*$",
        help="Regex to select sample IDs.",
    )
    parser.add_argument(
        "--taxon-id-mode",
        choices=["haplotype", "composite"],
        default="haplotype",
        help=(
            "How to name haplotype taxa: 'haplotype' uses HAP IDs; "
            "'composite' uses HAP + sample label."
        ),
    )
    parser.add_argument(
        "--weight-scale",
        type=int,
        default=0,
        help=(
            "If >0, replicate each haplotype taxon by round(freq_pct/100 * scale) "
            "to encode abundance for PopART node size. Default 0 disables replication."
        ),
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


def load_consensus_by_sample(consensus_dir: Path, sample_regex: str) -> dict[str, str]:
    sample_pattern = re.compile(sample_regex)
    mapping: dict[str, str] = {}
    for path in sorted(consensus_dir.glob("*.consensus.fasta")):
        sample_id = path.name[: -len(".consensus.fasta")]
        if not sample_pattern.match(sample_id):
            continue
        _, sequence = parse_fasta_single(path)
        mapping[sample_id] = sequence
    return mapping


def parse_mutation_tokens(raw: str) -> list[str]:
    text = (raw or "").strip()
    if not text or text == "REF":
        return []
    return [token.strip() for token in text.split(";") if token.strip()]


def normalize_haplotype_id(raw_haplotype: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_]+", "", raw_haplotype)
    token = token.replace("_", "")
    return token or "HAP"


def short_sample_label(sample_id: str) -> str:
    match = re.match(r"^INH_(\d+)_DPI_(R\d+)_.*$", sample_id)
    if match:
        dpi = match.group(1)
        replicate = match.group(2)
        return f"DPI{dpi}_{replicate}"
    return sample_id


def short_haplotype_label(haplotype_id: str) -> str:
    match = re.search(r"(\d+)", haplotype_id)
    if match:
        return f"H{int(match.group(1))}"
    return normalize_haplotype_id(haplotype_id)


def format_sequence_id(sample_id: str, haplotype_id: str, taxon_id_mode: str) -> str:
    hap_token = normalize_haplotype_id(haplotype_id)
    if taxon_id_mode == "composite":
        return f"{hap_token}__{short_sample_label(sample_id)}"
    return hap_token


def uniquify_taxon_id(base_id: str, used_ids: dict[str, int]) -> str:
    count = used_ids.get(base_id, 0) + 1
    used_ids[base_id] = count
    if count == 1:
        return base_id
    return f"{base_id}__dup{count - 1}"


def compute_weight_copies(frequency_percent: float, weight_scale: int) -> int:
    if weight_scale <= 0:
        return 1
    return max(1, int(round((frequency_percent / 100.0) * weight_scale)))


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


def apply_mutation_tokens(sequence: str, tokens: list[str]) -> tuple[str, int, int, int, int]:
    edited = list(sequence)
    applied = 0
    skipped_non_snp = 0
    skipped_ref_mismatch = 0
    skipped_out_of_range = 0

    for token in tokens:
        match = SNP_TOKEN_PATTERN.match(token)
        if not match:
            skipped_non_snp += 1
            continue

        ref = match.group("ref")
        alt = match.group("alt")
        pos = int(match.group("pos"))
        idx = pos - 1

        if idx < 0 or idx >= len(edited):
            skipped_out_of_range += 1
            continue
        if edited[idx] != ref:
            skipped_ref_mismatch += 1
            continue

        edited[idx] = alt
        applied += 1

    return "".join(edited), applied, skipped_non_snp, skipped_ref_mismatch, skipped_out_of_range


def main() -> None:
    args = parse_args()
    haplotype_csv = Path(args.haplotype_csv)
    consensus_dir = Path(args.consensus_dir)
    out_dir = Path(args.out_dir)

    fasta_dir = out_dir / "fasta"
    nexus_dir = out_dir / "nexus"
    meta_dir = out_dir / "metadata"
    logs_dir = out_dir / "logs"
    for directory in (fasta_dir, nexus_dir, meta_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    consensus_by_sample = load_consensus_by_sample(consensus_dir, args.sample_regex)
    if not consensus_by_sample:
        raise ValueError(f"No consensus FASTAs found in {consensus_dir}")

    excluded_days = parse_excluded_days(args.exclude_dpi_days)

    records: list[tuple[str, str]] = []
    haplotype_rows: list[HaplotypeRow] = []
    summaries: list[HaplotypeSummary] = []
    record_metadata: dict[str, tuple[str, str, str, float, int, int]] = {}
    used_taxon_ids: dict[str, int] = {}
    skipped_missing_consensus: list[str] = []
    skipped_below_threshold = 0

    with haplotype_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "haplotype_id", "frequency_percent", "mutations"}
        missing_cols = required - set(reader.fieldnames or [])
        if missing_cols:
            missing_str = ", ".join(sorted(missing_cols))
            raise ValueError(f"Missing required columns in {haplotype_csv}: {missing_str}")

        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            haplotype_id = (row.get("haplotype_id") or "").strip()
            frequency_text = (row.get("frequency_percent") or "").strip()
            mutations_text = (row.get("mutations") or "").strip()

            if not sample_id or not haplotype_id or not frequency_text:
                continue
            if infer_dpi_day(sample_id) in excluded_days:
                continue

            try:
                frequency_percent = float(frequency_text)
            except ValueError:
                continue

            if not (frequency_percent > args.min_frequency_percent):
                skipped_below_threshold += 1
                continue

            consensus = consensus_by_sample.get(sample_id)
            if consensus is None:
                skipped_missing_consensus.append(f"{sample_id}\t{haplotype_id}")
                continue

            tokens = parse_mutation_tokens(mutations_text)
            edited, applied, non_snp, ref_mismatch, out_of_range = apply_mutation_tokens(consensus, tokens)

            haplotype_rows.append(
                HaplotypeRow(
                    sample_id=sample_id,
                    haplotype_id=haplotype_id,
                    frequency_percent=frequency_percent,
                    edited_sequence=edited,
                    mutation_tokens_total=len(tokens),
                    mutation_tokens_applied=applied,
                    skipped_non_snp=non_snp,
                    skipped_ref_mismatch=ref_mismatch,
                    skipped_out_of_range=out_of_range,
                )
            )

    canonical_sequence_by_haplotype: dict[str, str] = {}
    canonical_rank_by_haplotype: dict[str, tuple[float, str]] = {}
    for row in haplotype_rows:
        rank = (row.frequency_percent, row.sample_id)
        previous = canonical_rank_by_haplotype.get(row.haplotype_id)
        if previous is None or rank > previous:
            canonical_rank_by_haplotype[row.haplotype_id] = rank
            canonical_sequence_by_haplotype[row.haplotype_id] = row.edited_sequence

    for row in haplotype_rows:
        canonical_sequence = canonical_sequence_by_haplotype[row.haplotype_id]
        base_id = format_sequence_id(row.sample_id, row.haplotype_id, args.taxon_id_mode)
        seq_id_base = uniquify_taxon_id(base_id, used_taxon_ids)
        copies = compute_weight_copies(row.frequency_percent, args.weight_scale)

        for copy_idx in range(1, copies + 1):
            seq_id = seq_id_base if copies == 1 else f"{seq_id_base}__r{copy_idx}"
            records.append((seq_id, canonical_sequence))
            record_metadata[seq_id] = (
                row.sample_id,
                infer_dpi_label(row.sample_id),
                row.haplotype_id,
                row.frequency_percent,
                copies,
                copy_idx,
            )

        summaries.append(
            HaplotypeSummary(
                sample_id=row.sample_id,
                haplotype_id=row.haplotype_id,
                frequency_percent=row.frequency_percent,
                mutation_tokens_total=row.mutation_tokens_total,
                mutation_tokens_applied=row.mutation_tokens_applied,
                skipped_non_snp=row.skipped_non_snp,
                skipped_ref_mismatch=row.skipped_ref_mismatch,
                skipped_out_of_range=row.skipped_out_of_range,
            )
        )

    if not records:
        raise ValueError("No haplotype sequences generated. Check filters and inputs.")

    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError("Generated haplotype sequences are not all same length")

    fasta_path = fasta_dir / "popart_cliquesnv_haplotypes_gt1pct.fasta"
    fasta_path.write_text(format_fasta(records), encoding="utf-8")

    nexus_path = nexus_dir / "popart_cliquesnv_haplotypes_gt1pct.nex"
    trait_labels = ["dpi3", "dpi5"]
    trait_rows = {
        seq_id: infer_dpi_one_hot(meta[0])
        for seq_id, meta in record_metadata.items()
    }
    write_nexus(records, nexus_path, trait_labels=trait_labels, trait_rows=trait_rows)

    metadata_path = meta_dir / "sequence_metadata.csv"
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sequence_id",
                "sample_id",
                "dpi",
                "haplotype_id",
                "frequency_percent",
                "weight_copies",
                "copy_index",
            ]
        )
        for seq_id, _ in records:
            sample_id, dpi, haplotype_id, freq_pct, copies, copy_idx = record_metadata[seq_id]
            writer.writerow(
                [
                    seq_id,
                    sample_id,
                    dpi,
                    haplotype_id,
                    f"{freq_pct:.3f}",
                    copies,
                    copy_idx,
                ]
            )

    traits_path = meta_dir / "sequence_traits.csv"
    with traits_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence_id", "dpi3", "dpi5"])
        for seq_id, _ in records:
            sample_id = record_metadata[seq_id][0]
            dpi3, dpi5 = infer_dpi_one_hot(sample_id)
            writer.writerow([seq_id, dpi3, dpi5])

    summary_path = logs_dir / "run_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "haplotype_id",
                "frequency_percent",
                "mutation_tokens_total",
                "mutation_tokens_applied",
                "skipped_non_snp",
                "skipped_ref_mismatch",
                "skipped_out_of_range",
            ]
        )
        for row in summaries:
            writer.writerow(
                [
                    row.sample_id,
                    row.haplotype_id,
                    f"{row.frequency_percent:.3f}",
                    row.mutation_tokens_total,
                    row.mutation_tokens_applied,
                    row.skipped_non_snp,
                    row.skipped_ref_mismatch,
                    row.skipped_out_of_range,
                ]
            )

    skipped_path = logs_dir / "skipped_missing_consensus.tsv"
    skipped_path.write_text(
        "sample_id\thaplotype_id\n"
        + ("\n".join(skipped_missing_consensus) + "\n" if skipped_missing_consensus else ""),
        encoding="utf-8",
    )

    run_info = logs_dir / "run_info.txt"
    run_info.write_text(
        "\n".join(
            [
                f"timestamp={datetime.now().isoformat(timespec='seconds')}",
                f"haplotype_csv={haplotype_csv}",
                f"consensus_dir={consensus_dir}",
                f"out_dir={out_dir}",
                f"min_frequency_percent_strict_gt={args.min_frequency_percent}",
                f"taxon_id_mode={args.taxon_id_mode}",
                f"weight_scale={args.weight_scale}",
                "sequence_merge_mode=canonical_by_haplotype_max_frequency",
                f"exclude_dpi_days={','.join(sorted(excluded_days)) if excluded_days else 'none'}",
                f"sequence_count={len(records)}",
                f"skipped_below_threshold={skipped_below_threshold}",
                f"skipped_missing_consensus={len(skipped_missing_consensus)}",
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
