#!/usr/bin/env python3
"""
Generate deterministic test fixtures for alphavirus variant analysis pipeline.
Reference: KP282671.1 length 11444 bp.
"""
import sys
import random
import gzip

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
        f.write("KP282671.1\tRefSeq\tgene\t50\t7549\t.\t+\t.\tID=gene0;Name=nonstructural_polyprotein;biotype=protein_coding\n")
        f.write("KP282671.1\tRefSeq\tmRNA\t50\t7549\t.\t+\t.\tID=rna0;Parent=gene0;Name=nonstructural_polyprotein;biotype=protein_coding\n")
        f.write("KP282671.1\tRefSeq\tCDS\t50\t7549\t.\t+\t0\tID=cds0;Parent=rna0;Name=nsP1-4;protein_id=YP_009118625.1\n")

    # Write primer BED
    with open("tests/data/primers.test.bed", "w") as f:
        f.write("KP282671.1\t10\t30\tprimer_1_F\t1\t+\n")
        f.write("KP282671.1\t120\t140\tprimer_1_R\t1\t-\n")
        f.write("KP282671.1\t200\t220\tprimer_2_F\t2\t+\n")
        f.write("KP282671.1\t280\t300\tprimer_2_R\t2\t-\n")

    # Generate SAM and FASTQ for sampleA and sampleB
    for sample, alt_pos in [("sampleA", 100), ("sampleB", 200)]:
        sam_file = f"tests/data/{sample}.sam"
        with open(sam_file, "w") as f:
            f.write("@HD\tVN:1.6\tSO:coordinate\n")
            f.write("@SQ\tSN:KP282671.1\tLN:11444\n")
            f.write(f"@RG\tID:{sample}\tSM:{sample}\n")
            for r in range(100):
                qname = f"read_{r}"
                pos = (r % 5) * 40 + 10
                rseq = seq[pos-1:pos+99]
                qual = "I" * 100
                f.write(f"{qname}\t0\tKP282671.1\t{pos}\t60\t100M\t*\t0\t0\t{rseq}\t{qual}\tRG:Z:{sample}\n")

        # FASTQ 1 & 2
        f1_path = f"tests/data/{sample}.test_1.fastq.gz"
        f2_path = f"tests/data/{sample}.test_2.fastq.gz"
        with gzip.open(f1_path, "wt") as f1, gzip.open(f2_path, "wt") as f2:
            for r in range(100):
                qname = f"read_{r}"
                pos1 = (r % 5) * 40 + 10
                pos2 = pos1 + 50
                rseq1 = seq[pos1-1:pos1+99]
                rseq2 = seq[pos2-1:pos2+99]
                qual = "I" * 100
                f1.write(f"@{qname}/1\n{rseq1}\n+\n{qual}\n")
                f2.write(f"@{qname}/2\n{rseq2}\n+\n{qual}\n")

if __name__ == "__main__":
    generate_ref()
