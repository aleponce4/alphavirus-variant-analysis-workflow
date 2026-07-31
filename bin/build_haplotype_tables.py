#!/usr/bin/env python3
"""
build_haplotype_tables.py
Parse CliqueSNV haplotype outputs, assign shared global HAP_n identifiers across samples,
compute per-sample quasispecies metrics (diversity, sharing), and output proportion-based CSV tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set

__version__ = "1.0.0"

# Header frequency extraction regexes
FREQ_PATTERNS = [
    re.compile(r"(?:freq|frequency)[=_:\s]*([0-9]*\.[0-9]+(?:[eE][-+]?[0-9]+)?)", re.IGNORECASE),
    re.compile(r"_([0-9]+\.[0-9]+)$"),
    re.compile(r"_([0-9]+\.[0-9]+)_"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build CliqueSNV haplotype tables with shared IDs and sample frequencies."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--cliquesnv-input",
        nargs="*",
        help="Path(s) to CliqueSNV FASTA files or directories containing them.",
    )
    parser.add_argument(
        "--reference-fasta",
        required=True,
        help="Path to viral reference FASTA file.",
    )
    parser.add_argument(
        "--samplesheet",
        help="Path to samplesheet CSV (sample,fastq_1,fastq_2,treatment) for metadata.",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for generated CSV and FASTA files.",
    )
    parser.add_argument(
        "--tier-report",
        type=float,
        default=0.01,
        help="Frequency threshold for report tier (default: 0.01).",
    )
    parser.add_argument(
        "--tier-exploratory",
        type=float,
        default=0.001,
        help="Frequency threshold for candidate tier floor (default: 0.001).",
    )
    return parser.parse_args()


def load_reference_sequence(fasta_path: Path) -> Tuple[str, str]:
    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA missing: {fasta_path}")
    header = ""
    seq_parts: List[str] = []
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].split()[0]
            else:
                seq_parts.append(line)
    full_seq = "".join(seq_parts).upper()
    if not full_seq:
        raise ValueError(f"Reference FASTA {fasta_path} contains no sequence data.")
    return header, full_seq


def load_samplesheet_metadata(samplesheet_path: Path | None) -> Dict[str, Dict[str, str]]:
    mapping: Dict[str, Dict[str, str]] = {}
    if not samplesheet_path or not samplesheet_path.exists():
        return mapping
    with samplesheet_path.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample_id = (row.get("sample") or row.get("sample_id") or "").strip()
            treatment = (row.get("treatment") or "unknown").strip()
            dpi = (row.get("dpi") or row.get("day") or row.get("dpi_day") or "").strip()
            if sample_id:
                mapping[sample_id] = {"treatment": treatment, "dpi": dpi}
    return mapping


def parse_fasta_entries(fasta_path: Path) -> List[Tuple[str, float, str]]:
    """Return list of (header_id, frequency, sequence) tuples."""
    entries: List[Tuple[str, str]] = []
    current_header = ""
    current_seq: List[str] = []

    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    entries.append((current_header, "".join(current_seq).upper()))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            entries.append((current_header, "".join(current_seq).upper()))

    parsed: List[Tuple[str, float, str]] = []
    for header, seq in entries:
        freq = None
        for pat in FREQ_PATTERNS:
            m = pat.search(header)
            if m:
                try:
                    freq = float(m.group(1))
                    break
                except ValueError:
                    pass
        if freq is None:
            freq = 1.0 / max(1, len(entries))
        parsed.append((header, freq, seq))

    # Normalize frequencies per file if total is far from 1.0
    total_freq = sum(f for _, f, _ in parsed)
    if total_freq > 0 and abs(total_freq - 1.0) > 0.05:
        parsed = [(h, f / total_freq, s) for h, f, s in parsed]

    return parsed


def compute_mutation_signature(ref_seq: str, hap_seq: str) -> Tuple[str, List[str]]:
    n = min(len(ref_seq), len(hap_seq))
    tokens: List[str] = []
    for idx in range(n):
        r = ref_seq[idx]
        h = hap_seq[idx]
        if h == "N" or h == r:
            continue
        pos_1b = idx + 1
        if h == "-":
            tokens.append(f"{r}{pos_1b}del")
        else:
            tokens.append(f"{r}{pos_1b}{h}")
    sig = ";".join(tokens) if tokens else "REF"
    return sig, tokens


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    freq_out_path = out_dir / "haplotype_frequency_by_sample.csv"
    summary_out_path = out_dir / "haplotype_summary.csv"
    fasta_out_path = out_dir / "haplotype_sequences.fasta"

    # Load reference & metadata
    try:
        _, ref_seq = load_reference_sequence(Path(args.reference_fasta))
    except Exception as e:
        sys.stderr.write(f"ERROR reading reference FASTA: {e}\n")
        sys.exit(1)

    metadata = load_samplesheet_metadata(Path(args.samplesheet) if args.samplesheet else None)

    # Collect FASTA input files
    fasta_files: List[Path] = []
    if args.cliquesnv_input:
        for item in args.cliquesnv_input:
            p = Path(item)
            if p.is_file():
                fasta_files.append(p)
            elif p.is_dir():
                fasta_files.extend(sorted(p.glob("**/*.fasta")))

    # Graceful degradation if no files or empty inputs
    if not fasta_files:
        sys.stdout.write("Notice: No CliqueSNV fasta files provided. Writing empty output tables.\n")
        with freq_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "treatment", "tool", "haplotype_id", "frequency", "mutations", "n_mutations", "tier"])
        with summary_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "treatment", "n_haplotypes", "n_shared", "max_frequency", "shannon_diversity"])
        fasta_out_path.touch()
        sys.exit(0)

    # Process samples
    sample_records: List[dict] = []
    unique_signatures: Dict[str, str] = {}  # sig -> seq
    signature_totals: Dict[str, float] = {}

    for fpath in sorted(fasta_files):
        sample_id = fpath.stem.replace(".cliquesnv", "").replace(".fasta", "")
        # Filter stub or empty fastas
        entries = parse_fasta_entries(fpath)
        if not entries:
            continue

        meta = metadata.get(sample_id, {})
        treatment = meta.get("treatment", "unknown") if isinstance(meta, dict) else str(meta)
        dpi = meta.get("dpi", "") if isinstance(meta, dict) else ""
        if treatment.strip().lower() == "mock":
            continue
        
        for header, freq, seq in entries:
            sig, tokens = compute_mutation_signature(ref_seq, seq)
            sample_records.append({
                "sample_id": sample_id,
                "treatment": treatment,
                "dpi": dpi,
                "tool": "CliqueSNV",
                "header": header,
                "frequency": freq,
                "mutations": sig,
                "n_mutations": len(tokens),
                "sequence": seq,
            })
            if sig not in unique_signatures:
                unique_signatures[sig] = seq
                signature_totals[sig] = 0.0
            signature_totals[sig] += freq

    if not sample_records:
        sys.stdout.write("Notice: CliqueSNV input files contained no valid sequences. Emitting empty outputs.\n")
        with freq_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "treatment", "dpi", "tool", "haplotype_id", "frequency", "mutations", "n_mutations", "tier"])
        with summary_out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "treatment", "n_haplotypes", "n_shared", "max_frequency", "shannon_diversity"])
        fasta_out_path.touch()
        sys.exit(0)

    # Assign shared HAP_n IDs deterministically sorted by overall dataset prevalence
    sorted_sigs = sorted(unique_signatures.keys(), key=lambda s: (-signature_totals[s], s))
    sig_to_hapid = {sig: f"HAP_{idx+1:04d}" for idx, sig in enumerate(sorted_sigs)}

    # Count sharing across samples
    hap_sample_counts: Dict[str, Set[str]] = {}
    for rec in sample_records:
        hap_id = sig_to_hapid[rec["mutations"]]
        rec["haplotype_id"] = hap_id
        hap_sample_counts.setdefault(hap_id, set()).add(rec["sample_id"])

    # Tier assignment helper
    def get_tier(freq: float) -> str:
        if freq >= args.tier_report:
            return "report"
        elif freq >= args.tier_exploratory:
            return "candidate"
        else:
            return "exploratory"

    # Write haplotype_frequency_by_sample.csv
    freq_rows = []
    for rec in sample_records:
        freq_rows.append({
            "sample_id": rec["sample_id"],
            "treatment": rec["treatment"],
            "dpi": rec["dpi"],
            "tool": rec["tool"],
            "haplotype_id": rec["haplotype_id"],
            "frequency": round(rec["frequency"], 6),
            "mutations": rec["mutations"],
            "n_mutations": rec["n_mutations"],
            "tier": get_tier(rec["frequency"]),
        })

    freq_rows.sort(key=lambda x: (x["sample_id"], -x["frequency"], x["haplotype_id"]))

    fieldnames_freq = ["sample_id", "treatment", "dpi", "tool", "haplotype_id", "frequency", "mutations", "n_mutations", "tier"]
    with freq_out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_freq)
        writer.writeheader()
        writer.writerows(freq_rows)

    # Write haplotype_sequences.fasta
    with fasta_out_path.open("w", encoding="utf-8") as fh:
        for sig in sorted_sigs:
            hap_id = sig_to_hapid[sig]
            seq = unique_signatures[sig]
            fh.write(f">{hap_id}\n{seq}\n")

    # Compute & write haplotype_summary.csv
    by_sample: Dict[str, List[dict]] = {}
    for row in freq_rows:
        by_sample.setdefault(row["sample_id"], []).append(row)

    summary_rows = []
    for sample_id in sorted(by_sample.keys()):
        rows = by_sample[sample_id]
        treatment = rows[0]["treatment"]
        n_hap = len(rows)
        n_shared = sum(1 for r in rows if len(hap_sample_counts[r["haplotype_id"]]) > 1)
        max_freq = max(r["frequency"] for r in rows) if rows else 0.0
        shannon = -sum(r["frequency"] * math.log(r["frequency"]) for r in rows if r["frequency"] > 0)
        summary_rows.append({
            "sample_id": sample_id,
            "treatment": treatment,
            "n_haplotypes": n_hap,
            "n_shared": n_shared,
            "max_frequency": round(max_freq, 6),
            "shannon_diversity": round(shannon, 6),
        })

    fieldnames_sum = ["sample_id", "treatment", "n_haplotypes", "n_shared", "max_frequency", "shannon_diversity"]
    with summary_out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_sum)
        writer.writeheader()
        writer.writerows(summary_rows)

    sys.stdout.write(f"Done. Processed {len(sample_records)} haplotype records across {len(by_sample)} samples.\n")


if __name__ == "__main__":
    main()
