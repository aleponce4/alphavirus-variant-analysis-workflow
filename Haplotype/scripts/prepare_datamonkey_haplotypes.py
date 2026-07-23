#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


STOP_CODONS = {"TAA", "TAG", "TGA"}
IUPAC_ALLOWED = set("ACGTRYSWKMBDHVN?-X")


@dataclass(frozen=True)
class GeneRegion:
    name: str
    start: int
    end: int


@dataclass(frozen=True)
class HaplotypeRecord:
    sample: str
    hap_name: str
    frequency: float
    sequence: str
    source_fasta: str
    tf_label: str

    @property
    def seq_id(self) -> str:
        freq_label = f"{self.frequency:.6f}".rstrip("0").rstrip(".")
        base = f"{self.sample}|{self.hap_name}|fr_{freq_label}|tf_{self.tf_label}"
        sanitized = re.sub(r"[^A-Za-z0-9_|.-]", "_", base)
        return re.sub(r"_+", "_", sanitized).strip("_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Datamonkey per-gene coding alignments from CliqueSNV haplotypes."
    )
    parser.add_argument(
        "--cliquesnv-dir",
        default="Haplotype/cliquesnv",
        help="CliqueSNV output root directory (contains per-sample tf_* folders).",
    )
    parser.add_argument(
        "--gff",
        default="Variant_discovery_pipeline/Input/Reference/VEEV_INH_fromGenbank.gff3",
        help="GFF3 file with CDS coordinates and gene_name attributes.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/consensus/DataMonkeyInput/datamonkeyhaplotypes",
        help="Output directory for Datamonkey haplotype FASTAs.",
    )
    parser.add_argument(
        "--sample-regex",
        default=r"^INH_\d+_DPI_.*$",
        help="Regex used to select sample directories in the CliqueSNV root.",
    )
    parser.add_argument(
        "--tf-label",
        default="auto",
        help=(
            "Threshold label to use, e.g. 0p01 or 0p001. "
            "Use 'auto' to pick the lowest available threshold per sample."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default="INH_AllDPI_Haplotypes",
        help="Prefix used for output FASTA filenames.",
    )
    parser.add_argument(
        "--min-frequency",
        type=float,
        default=0.0,
        help="Exclude haplotypes below this frequency.",
    )
    parser.add_argument(
        "--max-n-fraction",
        type=float,
        default=0.5,
        help="Exclude sequences from a gene if N fraction exceeds this threshold.",
    )
    parser.add_argument(
        "--mask-readthrough-site",
        default="",
        help=(
            "Optional genomic coordinate range to mask as NNN before stop checks, "
            "e.g. 5682-5684 for VEEV nsP3 opal readthrough codon."
        ),
    )
    return parser.parse_args()


def parse_gene_regions(gff_path: Path) -> list[GeneRegion]:
    regions: list[GeneRegion] = []
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            if fields[2] != "CDS":
                continue
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            gene_name_match = re.search(r"(?:^|;)gene_name=([^;]+)", attributes)
            id_match = re.search(r"(?:^|;)ID=cds-([^;]+)", attributes)
            if gene_name_match:
                name = gene_name_match.group(1)
            elif id_match:
                name = id_match.group(1)
            else:
                continue
            regions.append(GeneRegion(name=name, start=start, end=end))
    if not regions:
        raise ValueError(f"No CDS regions parsed from {gff_path}")
    return regions


def parse_site_range(site_text: str) -> tuple[int, int] | None:
    if not site_text:
        return None
    match = re.fullmatch(r"(\d+)-(\d+)", site_text.strip())
    if not match:
        raise ValueError(
            f"Invalid --mask-readthrough-site value '{site_text}'. Expected format START-END"
        )
    start = int(match.group(1))
    end = int(match.group(2))
    if end < start:
        raise ValueError("Readthrough site end must be >= start")
    if end - start + 1 != 3:
        raise ValueError("Readthrough masking currently supports exactly one codon (3 nt)")
    return start, end


def has_stop_codon(coding_nt: str, include_terminal: bool) -> bool:
    codon_count = len(coding_nt) // 3
    last_index = codon_count if include_terminal else max(codon_count - 1, 0)
    for codon_idx in range(last_index):
        codon = coding_nt[codon_idx * 3 : codon_idx * 3 + 3]
        if set(codon) <= {"A", "C", "G", "T"} and codon in STOP_CODONS:
            return True
    return False


def format_fasta_block(seq_id: str, sequence: str, width: int = 80) -> str:
    lines = [f">{seq_id}"]
    for index in range(0, len(sequence), width):
        lines.append(sequence[index : index + width])
    return "\n".join(lines)


def parse_frequency(header: str) -> float:
    match = re.search(r"_fr_([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", header)
    if match:
        return float(match.group(1))
    return 0.0


def parse_fasta_multi_sequence(fasta_path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if header is not None:
                    sequence = "".join(chunks).upper().replace("U", "T")
                    cleaned = "".join(base if base in IUPAC_ALLOWED else "N" for base in sequence)
                    records.append((header, cleaned))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        sequence = "".join(chunks).upper().replace("U", "T")
        cleaned = "".join(base if base in IUPAC_ALLOWED else "N" for base in sequence)
        records.append((header, cleaned))

    return records


def parse_tf_numeric(label: str) -> float:
    return float(label.replace("p", "."))


def choose_tf_label(sample_dir: Path, requested_label: str) -> str:
    tf_dirs = [d for d in sample_dir.iterdir() if d.is_dir() and d.name.startswith("tf_")]
    available = [d.name.replace("tf_", "") for d in tf_dirs if (d / "primary_only.fasta").exists()]
    if not available:
        raise FileNotFoundError(f"No tf_* directories with primary_only.fasta in {sample_dir}")

    if requested_label != "auto":
        if requested_label not in available:
            available_str = ", ".join(sorted(available))
            raise FileNotFoundError(
                f"Requested tf label '{requested_label}' not found in {sample_dir.name}. "
                f"Available: {available_str}"
            )
        return requested_label

    return sorted(available, key=parse_tf_numeric)[0]


def collect_haplotypes(
    cliquesnv_dir: Path,
    sample_regex: str,
    tf_label: str,
    min_frequency: float,
) -> list[HaplotypeRecord]:
    sample_pattern = re.compile(sample_regex)
    sample_dirs = sorted(
        d
        for d in cliquesnv_dir.iterdir()
        if d.is_dir() and sample_pattern.match(d.name)
    )
    if not sample_dirs:
        raise ValueError(
            f"No sample directories found in {cliquesnv_dir} matching regex '{sample_regex}'"
        )

    records: list[HaplotypeRecord] = []
    skipped_samples: list[str] = []
    for sample_dir in sample_dirs:
        try:
            chosen_tf = choose_tf_label(sample_dir, tf_label)
        except FileNotFoundError:
            skipped_samples.append(sample_dir.name)
            continue
        fasta_path = sample_dir / f"tf_{chosen_tf}" / "primary_only.fasta"
        for header, sequence in parse_fasta_multi_sequence(fasta_path):
            frequency = parse_frequency(header)
            if frequency < min_frequency:
                continue
            hap_name = header.split()[0]
            records.append(
                HaplotypeRecord(
                    sample=sample_dir.name,
                    hap_name=hap_name,
                    frequency=frequency,
                    sequence=sequence,
                    source_fasta=str(fasta_path),
                    tf_label=chosen_tf,
                )
            )

    if not records:
        raise ValueError("No haplotypes collected after filters.")
    if skipped_samples:
        skipped_str = ", ".join(sorted(skipped_samples))
        print(f"[info] Skipped samples without usable tf outputs: {skipped_str}")
    return records


def main() -> None:
    args = parse_args()
    cliquesnv_dir = Path(args.cliquesnv_dir)
    gff_path = Path(args.gff)
    out_dir = Path(args.out_dir)
    readthrough_site = parse_site_range(args.mask_readthrough_site)

    out_dir.mkdir(parents=True, exist_ok=True)
    gene_regions = parse_gene_regions(gff_path)

    haplotypes = collect_haplotypes(
        cliquesnv_dir=cliquesnv_dir,
        sample_regex=args.sample_regex,
        tf_label=args.tf_label,
        min_frequency=args.min_frequency,
    )

    summary_lines: list[str] = [
        "gene\tstart\tend\tinput_haplotypes\twritten_haplotypes\ttrimmed_terminal_stop\treadthrough_masked\texcluded"
    ]
    source_lines: list[str] = ["sequence_id\tsample\thaplotype\tfrequency\ttf_label\tsource_fasta"]

    for record in haplotypes:
        source_lines.append(
            "\t".join(
                [
                    record.seq_id,
                    record.sample,
                    record.hap_name,
                    f"{record.frequency:.10f}",
                    record.tf_label,
                    record.source_fasta,
                ]
            )
        )

    for gene in gene_regions:
        slices: dict[str, str] = {}
        exclusions: list[str] = []
        masked_count = 0

        for record in haplotypes:
            full_seq = record.sequence
            if len(full_seq) < gene.end:
                exclusions.append(f"{record.seq_id}:sequence_shorter_than_{gene.end}")
                continue

            coding = full_seq[gene.start - 1 : gene.end]

            if readthrough_site is not None:
                rt_start, rt_end = readthrough_site
                if gene.start <= rt_start and rt_end <= gene.end:
                    relative_start = rt_start - gene.start
                    relative_end = rt_end - gene.start + 1
                    if coding[relative_start:relative_end] in STOP_CODONS:
                        coding = coding[:relative_start] + "NNN" + coding[relative_end:]
                        masked_count += 1

            if len(coding) % 3 != 0:
                exclusions.append(f"{record.seq_id}:length_not_multiple_of_3")
                continue

            n_fraction = coding.count("N") / len(coding) if coding else 1.0
            if n_fraction > args.max_n_fraction:
                exclusions.append(f"{record.seq_id}:N_fraction_{n_fraction:.3f}")
                continue

            if has_stop_codon(coding, include_terminal=False):
                exclusions.append(f"{record.seq_id}:internal_stop")
                continue

            slices[record.seq_id] = coding

        trimmed_terminal_stop = "no"
        if slices:
            terminal_codons = {seq[-3:] for seq in slices.values() if len(seq) >= 3}
            if terminal_codons and all(codon in STOP_CODONS for codon in terminal_codons):
                trimmed_terminal_stop = "yes"
                slices = {seq_id: seq[:-3] for seq_id, seq in slices.items()}

        final_slices: dict[str, str] = {}
        for seq_id, coding in slices.items():
            if has_stop_codon(coding, include_terminal=True):
                exclusions.append(f"{seq_id}:stop_after_trim")
                continue
            final_slices[seq_id] = coding

        outfile = out_dir / f"{args.output_prefix}_{gene.name}.fasta"
        with outfile.open("w", encoding="utf-8") as out_handle:
            records = [
                format_fasta_block(seq_id, final_slices[seq_id])
                for seq_id in sorted(final_slices)
            ]
            out_handle.write("\n".join(records))
            if records:
                out_handle.write("\n")

        summary_lines.append(
            "\t".join(
                [
                    gene.name,
                    str(gene.start),
                    str(gene.end),
                    str(len(haplotypes)),
                    str(len(final_slices)),
                    trimmed_terminal_stop,
                    str(masked_count),
                    ";".join(exclusions) if exclusions else "none",
                ]
            )
        )

    summary_path = out_dir / "datamonkey_haplotype_prep_summary.tsv"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    source_map_path = out_dir / "datamonkey_haplotype_sequence_map.tsv"
    source_map_path.write_text("\n".join(source_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
