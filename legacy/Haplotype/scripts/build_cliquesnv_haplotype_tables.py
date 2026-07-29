#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable
from dataclasses import dataclass

SAMPLE_PATTERN = re.compile(r"^INH_(?P<dpi>\d+)_DPI_(?P<replicate>R\d+)_.*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build simple CliqueSNV haplotype tables with shared haplotype IDs, "
            "sample frequencies, and mutation lists."
        )
    )
    parser.add_argument(
        "--cliquesnv-dir",
        default="Haplotype/cliquesnv",
        help="Directory containing per-sample CliqueSNV outputs.",
    )
    parser.add_argument(
        "--reference-fasta",
        default="SNPGenie/input/reference/viral_only.fasta",
        help="Reference FASTA used to call haplotype mutations.",
    )
    parser.add_argument(
        "--aa-annotation-csv",
        default="Variant_discovery_pipeline/Analysis_output/VEEV_LoFreq_mutations.csv",
        help="CSV containing DNA_Change to amino-acid annotation mapping.",
    )
    parser.add_argument(
        "--annotation-gff",
        default="Variant_discovery_pipeline/Input/Reference/VEEV_INH_fromGenbank.gff3",
        help="GFF3 annotation used for CDS-based AA fallback when DNA_Change is not in CSV map.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/cliquesnv/Analysis",
        help="Directory where output TSV files are written.",
    )
    parser.add_argument(
        "--threshold-label",
        default="tf_0p01",
        help=(
            "Threshold directory name to include (default: tf_0p01). "
            "Use 'all' to include all threshold folders."
        ),
    )
    parser.add_argument(
        "--include-non-ok",
        action="store_true",
        help="Include rows even if JSON error is not 'none'.",
    )
    return parser.parse_args()


def read_single_fasta_sequence(path: Path) -> tuple[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {path}")

    name = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    raise ValueError(f"Expected single-sequence FASTA, found multiple in {path}")
                name = line[1:].strip()
                continue
            seq_parts.append(line)

    sequence = "".join(seq_parts).upper()
    if not sequence:
        raise ValueError(f"Reference FASTA has no sequence: {path}")

    return name, sequence


def parse_sample_id(sample_id: str) -> tuple[str, str]:
    match = SAMPLE_PATTERN.match(sample_id)
    if not match:
        return "", ""
    return match.group("dpi"), match.group("replicate")


def iter_sample_dirs(cliquesnv_dir: Path) -> Iterable[Path]:
    for sample_dir in sorted(cliquesnv_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        if sample_dir.name in {"logs", "Analysis"}:
            continue
        yield sample_dir


def mutation_tokens(reference: str, haplotype_seq: str) -> list[str]:
    seq = haplotype_seq.upper()
    n = min(len(reference), len(seq))
    tokens: list[str] = []

    for idx in range(n):
        ref_base = reference[idx]
        hap_base = seq[idx]

        if hap_base == "N":
            continue
        if hap_base == ref_base:
            continue

        pos_1b = idx + 1
        if hap_base == "-":
            tokens.append(f"{ref_base}{pos_1b}del")
        else:
            tokens.append(f"{ref_base}{pos_1b}{hap_base}")

    return tokens


SNP_TOKEN_PATTERN = re.compile(r"^(?P<ref>[ACGT])(?P<pos>\d+)(?P<alt>[ACGT])$")
DEL_TOKEN_PATTERN = re.compile(r"^(?P<ref>[ACGT])(?P<pos>\d+)del$")


def token_to_dna_change(token: str) -> str | None:
    snp = SNP_TOKEN_PATTERN.match(token)
    if snp:
        return f"{snp.group('pos')}{snp.group('ref')}>{snp.group('alt')}"
    deletion = DEL_TOKEN_PATTERN.match(token)
    if deletion:
        return f"{deletion.group('pos')}{deletion.group('ref')}>del"
    return None


@dataclass(frozen=True)
class CDSFeature:
    start: int
    end: int
    strand: str
    gene: str


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = value
    return attrs


def load_cds_features(gff_path: Path) -> list[CDSFeature]:
    features: list[CDSFeature] = []
    if not gff_path.exists():
        return features

    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            _, _, feature_type, start, end, _, strand, _, attrs_raw = cols
            if feature_type != "CDS":
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attrs = parse_gff_attributes(attrs_raw)
            gene = attrs.get("gene") or attrs.get("Name") or attrs.get("ID") or "CDS"
            if strand not in {"+", "-"}:
                strand = "+"
            features.append(CDSFeature(start=start_i, end=end_i, strand=strand, gene=gene))
    return features


DNA_CHANGE_SNP_PATTERN = re.compile(r"^(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
DNA_CHANGE_DEL_PATTERN = re.compile(r"^(?P<pos>\d+)(?P<ref>[ACGT])>del$")


def infer_effect_for_dna_change(dna_change: str, reference_seq: str, cds_features: list[CDSFeature]) -> str:
    snp = DNA_CHANGE_SNP_PATTERN.match(dna_change)
    deletion = DNA_CHANGE_DEL_PATTERN.match(dna_change)

    if deletion:
        pos = int(deletion.group("pos"))
        in_cds = any(c.start <= pos <= c.end for c in cds_features)
        return "frameshift_del" if in_cds else "noncoding_del"

    if not snp:
        return "NA"

    pos = int(snp.group("pos"))
    ref = snp.group("ref")
    alt = snp.group("alt")

    if pos < 1 or pos > len(reference_seq):
        return "NA"
    if reference_seq[pos - 1] != ref:
        return "NA"

    effects: list[str] = []
    for cds in cds_features:
        if not (cds.start <= pos <= cds.end):
            continue

        if cds.strand == "+":
            codon_start = cds.start + ((pos - cds.start) // 3) * 3
            codon_end = codon_start + 2
            aa_pos = ((pos - cds.start) // 3) + 1
        else:
            codon_end = cds.end - ((cds.end - pos) // 3) * 3
            codon_start = codon_end - 2
            aa_pos = ((cds.end - pos) // 3) + 1

        if codon_start < 1 or codon_end > len(reference_seq):
            continue

        codon_ref_genomic = reference_seq[codon_start - 1:codon_end]
        if len(codon_ref_genomic) != 3:
            continue

        idx_in_codon = pos - codon_start
        codon_alt_genomic = list(codon_ref_genomic)
        codon_alt_genomic[idx_in_codon] = alt
        codon_alt_genomic = "".join(codon_alt_genomic)

        if cds.strand == "+":
            codon_ref = codon_ref_genomic
            codon_alt = codon_alt_genomic
        else:
            codon_ref = reverse_complement(codon_ref_genomic)
            codon_alt = reverse_complement(codon_alt_genomic)

        aa_ref = CODON_TABLE.get(codon_ref, "X")
        aa_alt = CODON_TABLE.get(codon_alt, "X")
        aa_change = f"{aa_ref}{aa_pos}{aa_alt}"

        if aa_ref == aa_alt:
            effect = f"syn:{aa_change}"
        elif aa_alt == "*":
            effect = f"stop:{aa_change}"
        elif aa_ref == "*":
            effect = f"stop_lost:{aa_change}"
        else:
            effect = f"mis:{aa_change}"

        effects.append(f"{cds.gene}:{effect}")

    if effects:
        return "|".join(sorted(set(effects)))
    return "noncoding"


def load_aa_annotation_map(paths: Iterable[Path]) -> dict[str, str]:

    mapping: dict[str, set[str]] = {}
    for path in paths:
        if not path.exists():
            continue

        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                dna_change = (row.get("DNA_Change") or "").strip()
                if not dna_change:
                    continue

                mutation_type = (row.get("Mutation_Type") or "").strip().lower()
                aa_change = (row.get("Amino_Acid_Change") or "").strip()

                if mutation_type.startswith("syn"):
                    effect = f"syn:{aa_change}" if aa_change else "syn"
                elif mutation_type.startswith("missense"):
                    effect = f"mis:{aa_change}" if aa_change else "mis"
                elif mutation_type:
                    effect = f"{mutation_type}:{aa_change}" if aa_change else mutation_type
                else:
                    effect = aa_change if aa_change else "NA"

                if dna_change not in mapping:
                    mapping[dna_change] = set()
                mapping[dna_change].add(effect)

    collapsed: dict[str, str] = {}
    for dna_change, effects in mapping.items():
        collapsed[dna_change] = "|".join(sorted(effects))
    return collapsed


def aa_summary_for_mutations(
    mutations: str,
    aa_map: dict[str, str],
    reference_seq: str,
    cds_features: list[CDSFeature],
) -> str:
    if mutations == "REF":
        return "REF"

    parts: list[str] = []
    for token in mutations.split(";"):
        token = token.strip()
        if not token:
            continue
        dna_change = token_to_dna_change(token)
        if dna_change is None:
            parts.append("NA")
            continue
        mapped = aa_map.get(dna_change)
        if mapped and mapped != "NA":
            parts.append(mapped)
            continue
        parts.append(infer_effect_for_dna_change(dna_change, reference_seq, cds_features))
    return ";".join(parts) if parts else "NA"


def status_from_error(error_value: object) -> str:
    if isinstance(error_value, str) and error_value in {"none", "", "null"}:
        return "ok"
    return "non_ok"


def load_haplotype_rows(
    cliquesnv_dir: Path,
    reference_seq: str,
    threshold_label: str,
    include_non_ok: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for sample_dir in iter_sample_dirs(cliquesnv_dir):
        dpi, replicate_id = parse_sample_id(sample_dir.name)

        threshold_dirs: list[Path]
        if threshold_label == "all":
            threshold_dirs = sorted([d for d in sample_dir.iterdir() if d.is_dir() and d.name.startswith("tf_")])
        else:
            d = sample_dir / threshold_label
            threshold_dirs = [d] if d.is_dir() else []

        for tf_dir in threshold_dirs:
            json_files = sorted(tf_dir.glob("*.json"))
            if not json_files:
                continue
            json_path = json_files[-1]

            try:
                payload = json.loads(json_path.read_text(encoding="utf-8", errors="replace"))
            except Exception:
                continue

            json_error = payload.get("error", "none")
            status = status_from_error(json_error)
            if not include_non_ok and status != "ok":
                continue

            haplotypes = payload.get("haplotypes", [])
            if not isinstance(haplotypes, list):
                continue

            for idx, hap in enumerate(haplotypes, start=1):
                if not isinstance(hap, dict):
                    continue

                frequency = hap.get("frequency")
                sequence = hap.get("haplotype")
                if not isinstance(frequency, (int, float)):
                    continue
                if not isinstance(sequence, str) or not sequence:
                    continue

                mut_list = mutation_tokens(reference_seq, sequence)
                mut_signature = ";".join(mut_list) if mut_list else "REF"

                rows.append(
                    {
                        "sample_id": sample_dir.name,
                        "dpi": dpi,
                        "replicate_id": replicate_id,
                        "threshold_label": tf_dir.name,
                        "source_haplotype_name": hap.get("name", f"hap_{idx}"),
                        "json_error": json_error,
                        "status": status,
                        "frequency": float(frequency),
                        "frequency_percent": float(frequency) * 100.0,
                        "n_mutations": len(mut_list),
                        "mutations": mut_signature,
                    }
                )

    return rows


def assign_global_ids(rows: list[dict[str, object]]) -> None:
    signature_to_id: dict[str, str] = {}
    next_id = 1

    for row in rows:
        signature = str(row["mutations"])
        if signature not in signature_to_id:
            signature_to_id[signature] = f"HAP_{next_id:04d}"
            next_id += 1
        row["haplotype_id"] = signature_to_id[signature]


def write_frequency_table(rows: list[dict[str, object]], out_path: Path) -> None:
    fields = [
        "sample_id",
        "haplotype_id",
        "frequency_percent",
        "n_mutations",
        "mutations",
        "mutations_aa",
    ]

    sorted_rows = sorted(
        rows,
        key=lambda x: (
            str(x["sample_id"]),
            str(x["threshold_label"]),
            -float(x["frequency"]),
            str(x["haplotype_id"]),
        ),
    )

    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted_rows:
            formatted = dict(row)
            formatted["frequency_percent"] = f"{float(row['frequency_percent']):.3f}"
            writer.writerow({k: formatted.get(k, "") for k in fields})


def write_catalog_table(rows: list[dict[str, object]], out_path: Path) -> None:
    by_hap: dict[str, dict[str, object]] = {}

    for row in rows:
        hap_id = str(row["haplotype_id"])
        entry = by_hap.get(hap_id)
        if entry is None:
            entry = {
                "haplotype_id": hap_id,
                "n_mutations": row["n_mutations"],
                "mutations": row["mutations"],
                "mutations_aa": row["mutations_aa"],
                "sample_ids": set(),
                "max_frequency_percent": 0.0,
            }
            by_hap[hap_id] = entry

        entry["sample_ids"].add(str(row["sample_id"]))
        entry["max_frequency_percent"] = max(
            float(entry["max_frequency_percent"]),
            float(row["frequency_percent"]),
        )

    fields = [
        "haplotype_id",
        "n_mutations",
        "mutations",
        "mutations_aa",
        "n_samples",
        "sample_ids",
        "max_frequency_percent",
    ]

    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for hap_id in sorted(by_hap.keys()):
            entry = by_hap[hap_id]
            sample_ids = sorted(entry["sample_ids"])
            writer.writerow(
                {
                    "haplotype_id": hap_id,
                    "n_mutations": entry["n_mutations"],
                    "mutations": entry["mutations"],
                    "mutations_aa": entry["mutations_aa"],
                    "n_samples": len(sample_ids),
                    "sample_ids": ",".join(sample_ids),
                    "max_frequency_percent": f"{float(entry['max_frequency_percent']):.3f}",
                }
            )


def main() -> None:
    args = parse_args()

    cliquesnv_dir = Path(args.cliquesnv_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _, reference_seq = read_single_fasta_sequence(Path(args.reference_fasta))

    rows = load_haplotype_rows(
        cliquesnv_dir=cliquesnv_dir,
        reference_seq=reference_seq,
        threshold_label=args.threshold_label,
        include_non_ok=args.include_non_ok,
    )

    primary_aa_csv = Path(args.aa_annotation_csv)
    extra_aa_csv = primary_aa_csv.parent / "VEEV_variants_LoFreq.csv"
    aa_paths = [primary_aa_csv]
    if extra_aa_csv != primary_aa_csv:
        aa_paths.append(extra_aa_csv)

    aa_map = load_aa_annotation_map(aa_paths)
    cds_features = load_cds_features(Path(args.annotation_gff))
    for row in rows:
        mutations = str(row.get("mutations", ""))
        row["mutations_aa"] = aa_summary_for_mutations(
            mutations,
            aa_map,
            reference_seq,
            cds_features,
        )

    assign_global_ids(rows)

    for old_path in sorted(out_dir.glob("cliquesnv_haplotype_*.tsv")):
        old_path.unlink(missing_ok=True)
    for old_path in sorted(out_dir.glob("cliquesnv_haplotype_*.csv")):
        old_path.unlink(missing_ok=True)

    freq_out = out_dir / "cliquesnv_haplotype_frequency_table.csv"
    by_sample_out = out_dir / "cliquesnv_haplotype_frequency_by_sample.csv"

    write_catalog_table(rows, freq_out)
    write_frequency_table(rows, by_sample_out)

    print(f"Rows written: {len(rows)}")
    print(f"Wrote: {freq_out}")
    print(f"Wrote: {by_sample_out}")


if __name__ == "__main__":
    main()
