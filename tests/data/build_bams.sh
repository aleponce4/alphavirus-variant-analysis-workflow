#!/usr/bin/env bash
set -euo pipefail

python3 tests/data/generate_fixtures.py

for sample in sampleA sampleB; do
    samtools view -bS "tests/data/${sample}.sam" | samtools sort -o "tests/data/${sample}.test.bam" -
    samtools index "tests/data/${sample}.test.bam"
    rm -f "tests/data/${sample}.sam"
done

echo "BAM fixtures generated successfully."
