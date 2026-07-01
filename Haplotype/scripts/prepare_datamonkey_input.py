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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Datamonkey per-gene coding alignments from consensus FASTAs."
    )
    parser.add_argument(
        "--consensus-dir",
        default="Haplotype/consensus",
        help="Directory containing *.consensus.fasta files.",
    )
    parser.add_argument(
        "--gff",
        default="Variant_discovery_pipeline/Input/Reference/target_reference.gff3",
        help="GFF3 file with CDS coordinates and gene_name attributes.",
    )
    parser.add_argument(
        "--out-dir",
        default="Haplotype/consensus/DataMonkeyInput",
        help="Output directory for Datamonkey input FASTAs.",
    )
    parser.add_argument(
        "--sample-regex",
        default=r"^INH_\d+_DPI_.*\.consensus\.fasta$",
        help="Regex used to select consensus files to include.",
    )
    parser.add_argument(
        "--output-prefix",
        default="INH_AllDPI",
        help="Prefix used for output FASTA filenames.",
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
            "e.g. 5682-5684 for TARGET nsP3 opal readthrough codon."
        ),
    )
    return parser.parse_args()


def parse_fasta_single_sequence(fasta_path: Path) -> tuple[str, str]:
    header = None
    chunks: list[str] = []
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if header is not None:
                    raise ValueError(f"Multiple sequences found in {fasta_path}")
                header = line[1:].strip()
            else:
                chunks.append(line)
    if header is None:
        raise ValueError(f"No FASTA header found in {fasta_path}")
    sequence = "".join(chunks).upper().replace("U", "T")
    cleaned = "".join(base if base in IUPAC_ALLOWED else "N" for base in sequence)
    return header, cleaned


def sanitize_id(raw_id: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_]", "_", raw_id)
    return re.sub(r"_+", "_", sanitized).strip("_")


def parse_gene_regions(gff_path: Path) -> list[GeneRegion]:
    regions: list[GeneRegion] = []
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]
            if feature_type != "CDS":
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


def main() -> None:
    args = parse_args()
    consensus_dir = Path(args.consensus_dir)
    gff_path = Path(args.gff)
    out_dir = Path(args.out_dir)
    sample_pattern = re.compile(args.sample_regex)
    readthrough_site = parse_site_range(args.mask_readthrough_site)

    out_dir.mkdir(parents=True, exist_ok=True)
    gene_regions = parse_gene_regions(gff_path)

    consensus_files = sorted(
        fasta for fasta in consensus_dir.glob("*.consensus.fasta") if sample_pattern.match(fasta.name)
    )
    if not consensus_files:
        raise ValueError(
            f"No consensus FASTAs found in {consensus_dir} matching {args.sample_regex}"
        )

    sample_sequences: dict[str, str] = {}
    sample_source: dict[str, str] = {}
    for fasta_path in consensus_files:
        header, sequence = parse_fasta_single_sequence(fasta_path)
        sequence_id = sanitize_id(header or fasta_path.stem.replace(".consensus", ""))
        if not sequence_id:
            raise ValueError(f"Could not derive sequence ID from {fasta_path}")
        sample_sequences[sequence_id] = sequence
        sample_source[sequence_id] = fasta_path.name

    summary_lines: list[str] = [
        "gene\tstart\tend\tinput_sequences\twritten_sequences\ttrimmed_terminal_stop\treadthrough_masked\texcluded"
    ]

    for gene in gene_regions:
        slices: dict[str, str] = {}
        exclusions: list[str] = []
        masked_count = 0

        for seq_id, full_seq in sample_sequences.items():
            if len(full_seq) < gene.end:
                exclusions.append(f"{seq_id}:sequence_shorter_than_{gene.end}")
                continue
            coding = full_seq[gene.start - 1 : gene.end]

            if readthrough_site is not None:
                rt_start, rt_end = readthrough_site
                if gene.start <= rt_start and rt_end <= gene.end:
                    relative_start = rt_start - gene.start
                    relative_end = rt_end - gene.start + 1
                    if coding[relative_start:relative_end] in STOP_CODONS:
                        coding = (
                            coding[:relative_start]
                            + "NNN"
                            + coding[relative_end:]
                        )
                        masked_count += 1

            if len(coding) % 3 != 0:
                exclusions.append(f"{seq_id}:length_not_multiple_of_3")
                continue

            n_fraction = coding.count("N") / len(coding) if coding else 1.0
            if n_fraction > args.max_n_fraction:
                exclusions.append(f"{seq_id}:N_fraction_{n_fraction:.3f}")
                continue

            if has_stop_codon(coding, include_terminal=False):
                exclusions.append(f"{seq_id}:internal_stop")
                continue

            slices[seq_id] = coding

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
                    str(len(sample_sequences)),
                    str(len(final_slices)),
                    trimmed_terminal_stop,
                    str(masked_count),
                    ";".join(exclusions) if exclusions else "none",
                ]
            )
        )

    summary_path = out_dir / "datamonkey_prep_summary.tsv"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

