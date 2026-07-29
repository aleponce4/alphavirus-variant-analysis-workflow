# Test Fixtures

This directory contains test datasets used by `nf-test` and the `-profile test` execution profile:

- `viral_ref.test.fasta`: Synthetic VEEV-derived sequence (~11.5 kb) with contig `KP282671.1`.
- `viral_ref.test.gff3`: Matching GFF3 containing gene, mRNA, and CDS annotations.
- `sampleA.test_1.fastq.gz` / `sampleA.test_2.fastq.gz`: Paired-end FASTQ reads for sampleA.
- `sampleB.test_1.fastq.gz` / `sampleB.test_2.fastq.gz`: Paired-end FASTQ reads for sampleB.
- `sampleA.test.bam` / `sampleA.test.bam.bai`: Test BAM with seeded variants.
- `sampleB.test.bam` / `sampleB.test.bam.bai`: Test BAM with seeded variants.
- `primers.test.bed`: Synthetic primer BED file for amplicon protocol testing.

Fixtures are generated deterministically using `tests/data/generate_fixtures.sh`.
