#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_attributes(attr_text: str) -> dict:
    parsed = {}
    for item in attr_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        elif " " in item:
            key, value = item.split(" ", 1)
            parsed[key] = value.strip().strip('"')
    return parsed


def pick_gene_name(attrs: dict) -> str:
    for key in ("gene_name", "gene", "Name", "locus_tag", "product", "ID", "Parent"):
        value = attrs.get(key)
        if value:
            return value
    return "unknown_gene"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert GFF3 CDS records to SNPGenie-compatible GTF with gene_id tags."
    )
    parser.add_argument("--gff3", required=True, help="Input GFF3 file")
    parser.add_argument("--gtf", required=True, help="Output GTF file")
    args = parser.parse_args()

    gff3_path = Path(args.gff3)
    gtf_path = Path(args.gtf)

    cds_count = 0
    with gff3_path.open("r", encoding="utf-8") as in_handle, gtf_path.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attributes = cols
            if feature != "CDS":
                continue

            parsed_attrs = parse_attributes(attributes)
            gene_name = pick_gene_name(parsed_attrs)

            if gene_name.startswith("cds-"):
                gene_name = gene_name.replace("cds-", "", 1)
            if gene_name.startswith("rna-"):
                gene_name = gene_name.replace("rna-", "", 1)
            if gene_name.startswith("gene-"):
                gene_name = gene_name.replace("gene-", "", 1)

            gtf_attrs = f'gene_id "{gene_name}"; transcript_id "{gene_name}";'
            out_handle.write(
                "\t".join(
                    [
                        seqid,
                        source if source else "SNPGenie",
                        "CDS",
                        start,
                        end,
                        score if score else ".",
                        strand,
                        phase if phase else ".",
                        gtf_attrs,
                    ]
                )
                + "\n"
            )
            cds_count += 1

    if cds_count == 0:
        raise SystemExit("No CDS records found in input GFF3.")

    print(f"Wrote {cds_count} CDS records to {gtf_path}")


if __name__ == "__main__":
    main()
