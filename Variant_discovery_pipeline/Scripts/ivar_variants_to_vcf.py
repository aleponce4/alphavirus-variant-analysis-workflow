#!/usr/bin/env python3
"""
Convert iVar TSV files to VCF format for downstream annotation

This script converts iVar variant TSV files to proper VCF format
that can be processed by bcftools csq for annotation.

Usage: python ivar_variants_to_vcf.py <input.tsv> <output.vcf> <reference.fasta>
"""

import sys
import os
from datetime import datetime

def parse_ivar_tsv(tsv_file):
    """Parse iVar TSV file and extract variant information."""
    variants = []
    
    with open(tsv_file, 'r') as f:
        # Skip header line
        header = f.readline().strip().split('\t')
        
        # Expected columns in iVar TSV:
        # REGION, POS, REF, ALT, REF_DP, REF_RV, REF_QUAL, ALT_DP, ALT_RV, ALT_QUAL, ALT_FREQ, TOTAL_DP, PVAL, PASS
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            
            if len(fields) >= 12:  # Minimum required fields
                variant = {
                    'CHROM': fields[0],
                    'POS': int(fields[1]),
                    'REF': fields[2],
                    'ALT': fields[3],
                    'REF_DP': int(fields[4]) if fields[4].isdigit() else 0,
                    'ALT_DP': int(fields[7]) if fields[7].isdigit() else 0,
                    'ALT_FREQ': float(fields[10]) if fields[10].replace('.', '').isdigit() else 0.0,
                    'TOTAL_DP': int(fields[11]) if fields[11].isdigit() else 0,
                    'PVAL': float(fields[12]) if len(fields) > 12 and fields[12].replace('.', '').replace('e-', '').replace('E-', '').isdigit() else 1.0,
                    'PASS': fields[13] if len(fields) > 13 else 'TRUE'
                }
                variants.append(variant)
    
    return variants

def write_vcf_header(output_file, reference_file, sample_name):
    """Write VCF header."""
    with open(output_file, 'w') as f:
        # VCF version
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write(f"##source=ivar_variants_to_vcf.py\n")
        f.write(f"##reference={reference_file}\n")
        
        # INFO fields
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        
        # FORMAT fields
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
        f.write("##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description=\"Alternative allele frequency\">\n")
        
        # FILTER fields
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write("##FILTER=<ID=FAIL,Description=\"Failed quality filters\">\n")
        
        # Contig header for viral reference
        f.write("##contig=<ID=VEEV_INH>\n")
        
        # Column headers
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sample_name))

def convert_to_vcf(variants, output_file, reference_file, sample_name):
    """Convert parsed variants to VCF format."""
    
    # Write header
    write_vcf_header(output_file, reference_file, sample_name)
    
    # Write variants
    with open(output_file, 'a') as f:
        for variant in variants:
            # Basic fields
            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']
            qual = 60  # Default quality score
            
            # Filter
            filter_field = "PASS" if variant['PASS'] == 'TRUE' else "FAIL"
            
            # INFO field
            total_dp = variant['TOTAL_DP']
            alt_freq = variant['ALT_FREQ']
            info = f"DP={total_dp};AF={alt_freq}"
            
            # FORMAT field
            format_field = "GT:DP:AD:ALT_FREQ"
            
            # Sample data - use homozygous alt to avoid phasing issues
            genotype = "1/1"  # Homozygous alternate call
            ref_dp = variant['REF_DP']
            alt_dp = variant['ALT_DP']
            sample_data = f"{genotype}:{total_dp}:{ref_dp},{alt_dp}:{alt_freq}"
            
            # Write VCF line (ensure proper tab separation)
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_data}\n"
            f.write(vcf_line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python ivar_variants_to_vcf.py <input.tsv> <output.vcf> <reference.fasta>")
        print("Example: python ivar_variants_to_vcf.py sample.tsv sample.vcf reference.fasta")
        sys.exit(1)
    
    input_tsv = sys.argv[1]
    output_vcf = sys.argv[2]
    reference = sys.argv[3]
    
    # Check input files exist
    if not os.path.exists(input_tsv):
        print(f"ERROR: Input TSV file not found: {input_tsv}")
        sys.exit(1)
    
    if not os.path.exists(reference):
        print(f"ERROR: Reference file not found: {reference}")
        sys.exit(1)
    
    # Extract sample name from input filename
    sample_name = os.path.splitext(os.path.basename(input_tsv))[0]
    
    print(f"Converting {input_tsv} to VCF format...")
    print(f"Sample name: {sample_name}")
    print(f"Reference: {reference}")
    print(f"Output: {output_vcf}")
    
    # Parse TSV file
    try:
        variants = parse_ivar_tsv(input_tsv)
        print(f"Found {len(variants)} variants")
        
        if len(variants) == 0:
            print("WARNING: No variants found in input file")
        
        # Convert to VCF
        convert_to_vcf(variants, output_vcf, reference, sample_name)
        print(f"VCF file created: {output_vcf}")
        
    except Exception as e:
        print(f"ERROR: Failed to convert TSV to VCF: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()