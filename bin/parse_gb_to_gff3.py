#!/usr/bin/env python3
import sys
import os
import re
import argparse
from pathlib import Path

__version__ = "1.0.0"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse GenBank (.gb) file and generate GFF3 with mature peptides as individual CDS features."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--gb", required=True, help="Input GenBank (.gb) file")
    parser.add_argument("--gff3", required=True, help="Output GFF3 file")
    return parser.parse_args()

def parse_genbank_peptides(gb_path: Path) -> tuple[str, list[dict]]:
    peptides = []
    locus_name = None
    with gb_path.open("r", encoding="utf-8") as f:
        current_feature = None
        feature_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith("LOCUS"):
                tokens = line.split()
                if len(tokens) >= 2:
                    locus_name = tokens[1]
            elif line.startswith("VERSION"):
                tokens = line.split()
                if len(tokens) >= 2:
                    locus_name = tokens[1]
            elif line.startswith("     mat_peptide") or line.startswith("     CDS") or line.startswith("     source"):
                if current_feature == "mat_peptide":
                    pep = process_peptide_block(feature_lines)
                    if pep:
                        peptides.append(pep)
                current_feature = line.strip().split()[0]
                feature_lines = [line]
            elif line.startswith("      ") and current_feature:
                feature_lines.append(line)
            elif line.startswith("ORIGIN") or line.startswith("//"):
                if current_feature == "mat_peptide":
                    pep = process_peptide_block(feature_lines)
                    if pep:
                        peptides.append(pep)
                current_feature = None
    return locus_name, peptides

def process_peptide_block(lines: list[str]) -> dict | None:
    block_text = " ".join(l.strip() for l in lines)
    loc_match = re.search(r"mat_peptide\s+(\d+)\.\.(\d+)", block_text)
    if not loc_match:
        return None
    start = int(loc_match.group(1))
    end = int(loc_match.group(2))
    product_match = re.search(r'/product="([^"]+)"', block_text)
    if not product_match:
        return None
    product_name = product_match.group(1)
    std_name = product_name.lower().replace("envelope glycoprotein", "").replace("membrane protein", "").replace("protein", "").replace("precursor", "").strip().replace(" ", "_")
    if "nsp1" in std_name: std_name = "nsP1"
    elif "nsp2" in std_name: std_name = "nsP2"
    elif "nsp3" in std_name: std_name = "nsP3"
    elif "nsp4" in std_name: std_name = "nsP4"
    elif "capsid" in std_name: std_name = "Capsid"
    elif "e3" in std_name: std_name = "E3"
    elif "e2" in std_name: std_name = "E2"
    elif "6k" in std_name: std_name = "6K"
    elif "e1" in std_name: std_name = "E1"
    return {"start": start, "end": end, "product": product_name, "name": std_name}

def write_gff3(output_path: Path, seqid: str, peptides: list[dict]) -> None:
    with output_path.open("w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
        for idx, pep in enumerate(peptides, start=1):
            name = pep["name"]
            start = pep["start"]
            end = pep["end"]
            product = pep["product"]
            f.write(f"{seqid}\tGenBank\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:{name};Name={name};biotype=protein_coding\n")
            f.write(f"{seqid}\tGenBank\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:{name};Parent=gene:{name};Name={name};biotype=protein_coding\n")
            f.write(f"{seqid}\tGenBank\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds:{name};Parent=transcript:{name};Name={name};biotype=protein_coding;product={product}\n")
            f.write(f"{seqid}\tGenBank\texon\t{start}\t{end}\t.\t+\t.\tID=exon:{name};Parent=transcript:{name}\n")

def main() -> None:
    args = parse_args()
    gb_path = Path(args.gb)
    gff3_path = Path(args.gff3)
    if not gb_path.exists():
        print(f"GenBank file not found: {gb_path}", file=sys.stderr)
        sys.exit(1)
    seqid, peptides = parse_genbank_peptides(gb_path)
    if not seqid:
        seqid = gb_path.stem
    if not peptides:
        print("ERROR: No mature peptide features found in GenBank file.", file=sys.stderr)
        sys.exit(1)
    write_gff3(gff3_path, seqid, peptides)
    print(f"Successfully wrote GFF3 to: {gff3_path}")

if __name__ == "__main__":
    main()
