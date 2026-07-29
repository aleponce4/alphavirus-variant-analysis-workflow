#!/usr/bin/env python3
import sys
import os
import argparse
from datetime import datetime

__version__ = "1.0.0"

def parse_ivar_tsv(tsv_file):
    """Parse iVar TSV file and extract variant information."""
    variants = []
    with open(tsv_file, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 12:
                variant = {
                    'CHROM': fields[0],
                    'POS': int(fields[1]),
                    'REF': fields[2],
                    'ALT': fields[3],
                    'REF_DP': int(fields[4]) if fields[4].isdigit() else 0,
                    'ALT_DP': int(fields[7]) if fields[7].isdigit() else 0,
                    'ALT_FREQ': float(fields[10]) if fields[10].replace('.', '', 1).isdigit() else 0.0,
                    'TOTAL_DP': int(fields[11]) if fields[11].isdigit() else 0,
                    'PVAL': float(fields[12]) if len(fields) > 12 and fields[12].replace('.', '', 1).replace('e-', '', 1).replace('E-', '', 1).isdigit() else 1.0,
                    'PASS': fields[13] if len(fields) > 13 else 'TRUE'
                }
                variants.append(variant)
    return variants

def write_vcf_header(output_file, reference_file, sample_name):
    """Write VCF header."""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=ivar_variants_to_vcf.py\n")
        f.write(f"##reference={reference_file}\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
        f.write("##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description=\"Alternative allele frequency\">\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write("##FILTER=<ID=FAIL,Description=\"Failed quality filters\">\n")
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

def convert_to_vcf(variants, output_file, reference_file, sample_name):
    """Convert parsed variants to VCF format."""
    write_vcf_header(output_file, reference_file, sample_name)
    with open(output_file, 'a', encoding='utf-8') as f:
        for variant in variants:
            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']
            qual = 60
            filter_field = "PASS" if variant['PASS'] == 'TRUE' else "FAIL"
            total_dp = variant['TOTAL_DP']
            alt_freq = variant['ALT_FREQ']
            info = f"DP={total_dp};AF={alt_freq}"
            format_field = "GT:DP:AD:ALT_FREQ"
            genotype = "1/1"
            ref_dp = variant['REF_DP']
            alt_dp = variant['ALT_DP']
            sample_data = f"{genotype}:{total_dp}:{ref_dp},{alt_dp}:{alt_freq}"
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_data}\n"
            f.write(vcf_line)

def main():
    parser = argparse.ArgumentParser(description="Convert iVar TSV output to VCF format.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input-tsv", required=True, help="Input iVar variants TSV file")
    parser.add_argument("--output-vcf", required=True, help="Output VCF file")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA file")

    args = parser.parse_args()

    if not os.path.exists(args.input_tsv):
        print(f"ERROR: Input TSV file not found: {args.input_tsv}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.reference_fasta):
        print(f"ERROR: Reference file not found: {args.reference_fasta}", file=sys.stderr)
        sys.exit(1)

    sample_name = os.path.splitext(os.path.basename(args.input_tsv))[0]
    variants = parse_ivar_tsv(args.input_tsv)
    convert_to_vcf(variants, args.output_vcf, args.reference_fasta, sample_name)
    print(f"Successfully converted {len(variants)} variants from {args.input_tsv} to {args.output_vcf}")

if __name__ == "__main__":
    main()
