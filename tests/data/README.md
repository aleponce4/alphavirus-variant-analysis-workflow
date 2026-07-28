# Test Fixtures

This directory contains test datasets used by `nf-test` and the `-profile test` execution profile:

- `viral_ref.test.fasta`: Synthetic VEEV-derived sequence (~11.5 kb) with contig `KP282671.1`.
- `viral_ref.test.gff3`: Matching GFF3 containing CDS annotations.
- `sampleA.test.bam` / `sampleA.test.bam.bai`: Test BAM with seeded variants.
- `sampleB.test.bam` / `sampleB.test.bam.bai`: Test BAM with seeded variants.

Fixtures are generated deterministically using `tests/data/generate_fixtures.sh`.
