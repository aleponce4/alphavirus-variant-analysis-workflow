#!/bin/bash
# Fast viral read count using idxstats (reads .bai index, not the BAM)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-env

BAM_DIR="/mnt/d/RNAseq/Alphavirus/045_048_nftcore/variant_analysis/input/mouse_veev/BAMs"
CONTIG="KP282671.1"

echo "sample	viral_reads	mapped_total"
echo "------	-----------	------------"
for bam in "$BAM_DIR"/*.bam; do
    s=$(basename "$bam" .bam)
    # idxstats outputs: contig  length  mapped  unmapped
    viral=$(samtools idxstats "$bam" 2>/dev/null | grep "^${CONTIG}" | awk '{print $3}')
    total=$(samtools idxstats "$bam" 2>/dev/null | awk '{sum+=$3} END{print sum}')
    echo "${s}	${viral:-0}	${total:-0}"
done
