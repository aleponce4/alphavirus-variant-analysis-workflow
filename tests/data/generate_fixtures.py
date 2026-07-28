#!/usr/bin/env python3
"""
Generate deterministic test fixtures for alphavirus variant analysis pipeline.
Reference: KP282671.1 length 11444 bp.
"""
import sys
import random

def generate_ref():
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    # Build 11444 bp synthetic sequence
    seq = ''.join(random.choice(bases) for _ in range(11444))
    
    # Write FASTA
    with open("tests/data/viral_ref.test.fasta", "w") as f:
        f.write(">KP282671.1 Venezuelan equine encephalitis virus strain TC-83, complete genome\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
            
    # Write GFF3
    with open("tests/data/viral_ref.test.gff3", "w") as f:
        f.write("##gff-version 3\n")
        f.write("##sequence-region KP282671.1 1 11444\n")
        f.write("KP282671.1\tRefSeq\tregion\t1\t11444\t.\t+\t.\tID=region0;Name=KP282671.1\n")
        f.write("KP282671.1\tRefSeq\tgene\t50\t7500\t.\t+\t.\tID=gene0;Name=nonstructural_polyprotein\n")
        f.write("KP282671.1\tRefSeq\tCDS\t50\t7500\t.\t+\t0\tID=cds0;Parent=gene0;Name=nsP1-4\n")

    # Generate SAM for sampleA and sampleB
    for sample, alt_pos in [("sampleA", 100), ("sampleB", 200)]:
        sam_file = f"tests/data/{sample}.sam"
        with open(sam_file, "w") as f:
            f.write("@HD\tVN:1.6\tSO:coordinate\n")
            f.write("@SQ\tSN:KP282671.1\tLN:11444\n")
            f.write(f"@RG\tID:{sample}\tSM:{sample}\n")
            # Generate 100 aligned reads covering region 1..300
            for r in range(100):
                qname = f"read_{r}"
                pos = (r % 5) * 40 + 10
                rseq = seq[pos-1:pos+99]
                qual = "I" * 100
                f.write(f"{qname}\t0\tKP282671.1\t{pos}\t60\t100M\t*\t0\t0\t{rseq}\t{qual}\tRG:Z:{sample}\n")

if __name__ == "__main__":
    generate_ref()
