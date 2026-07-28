#!/usr/bin/env bash
set -euo pipefail

python3 tests/data/generate_fixtures.py

# Convert SAM to sorted BAM and index using docker samtools
CONTAINER_SAMTOOLS="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

for sample in sampleA sampleB; do
    samtools view -bS tests/data/${sample}.sam | samtools sort -o tests/data/${sample}.test.bam -
    samtools index tests/data/${sample}.test.bam
    rm -f tests/data/${sample}.sam
done

echo "Fixtures successfully generated."
