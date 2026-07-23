#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-env

BAM_DIR="/mnt/d/RNAseq/Alphavirus/045_048_nftcore/variant_analysis/input/mouse_veev/BAMs"
CONTIG="KP282671.1"

# Check a few samples: challenge DPI 3-4 (should have virus)
for s in s157 s158 s173 s174 s357 s373; do
    bam="$BAM_DIR/${s}.bam"
    reads=$(samtools view -c -r "$CONTIG" "$bam" 2>/dev/null || echo 0)
    if [ "$reads" -gt 0 ]; then
        depth_stats=$(samtools depth -a -r "$CONTIG" "$bam" 2>/dev/null | awk '
            BEGIN { min=999999; max=0; sum=0; n=0; a100=0; a1000=0 }
            { d=$3; sum+=d; n++; if(d<min)min=d; if(d>max)max=d; if(d>=100)a100++; if(d>=1000)a1000++ }
            END { printf "mean=%.0f min=%d max=%d >=100x=%.1f%% >=1000x=%.1f%%", sum/n, min, max, a100/n*100, a1000/n*100 }
        ')
        echo "$s: $reads reads | $depth_stats"
    else
        echo "$s: 0 viral reads"
    fi
done
