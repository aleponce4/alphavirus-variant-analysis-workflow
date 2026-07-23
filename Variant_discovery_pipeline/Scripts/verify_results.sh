#!/bin/bash
echo "=== LoFreq results ==="
ls -d ~/variant_results/mouse_veev/LoFreq/*/ 2>/dev/null | wc -l
echo "sample dirs"
echo ""
echo "=== VCFs ==="
find ~/variant_results/mouse_veev/LoFreq -name 'variants.filtered.vcf.gz' 2>/dev/null | wc -l
echo "filtered VCFs"
echo ""
echo "=== Annotated ==="
ls ~/variant_results/mouse_veev/Annotated_variants/LoFreq/*.vcf 2>/dev/null | wc -l
echo "annotated VCFs"
echo ""
echo "=== Variant counts ==="
for d in ~/variant_results/mouse_veev/LoFreq/*/; do
    s=$(basename "$d")
    v=$(zcat "$d/variants.filtered.vcf.gz" 2>/dev/null | grep -v '^#' | wc -l)
    echo "$s: $v variants"
done
