#!/bin/bash
# Quick coverage check on a single sample
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-env

VARIANT_ROOT="/mnt/d/RNAseq/Alphavirus/045_048_nftcore/variant_analysis"
BAM="$VARIANT_ROOT/input/mouse_veev/BAMs/s314.bam"
CONTIG="KP282671.1"

echo "Checking coverage on: $(basename $BAM)"
echo "Contig: $CONTIG"
echo ""

samtools depth -a -r "$CONTIG" "$BAM" | awk '
BEGIN { min=999999; max=0; sum=0; n=0; a100=0; a1000=0; a5000=0 }
{
    d = $3
    sum += d
    n++
    if (d < min) min = d
    if (d > max) max = d
    if (d >= 100) a100++
    if (d >= 1000) a1000++
    if (d >= 5000) a5000++
}
END {
    if (n == 0) { print "No coverage data"; exit }
    printf "Positions:    %d\n", n
    printf "Mean depth:   %.1f\n", sum/n
    printf "Min depth:    %d\n", min
    printf "Max depth:    %d\n", max
    printf "%% >= 100x:   %.1f\n", a100/n*100
    printf "%% >= 1000x:  %.1f\n", a1000/n*100
    printf "%% >= 5000x:  %.1f\n", a5000/n*100
}
'
