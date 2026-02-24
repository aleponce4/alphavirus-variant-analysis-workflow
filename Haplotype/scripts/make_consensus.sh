#!/bin/bash
# ============================================================
# make_consensus.sh
# Build a consensus FASTA for each LoFreq sample.
#
#   1. Filter VCF to AF > 0.5  (major-allele consensus)
#   2. Mask low-depth positions to N  (depth < MIN_DEPTH)
#   3. bcftools consensus → FASTA
#
# Usage:
#   cd Variant_discovery_pipeline
#   ../Haplotype/scripts/make_consensus.sh            # all samples
#   ../Haplotype/scripts/make_consensus.sh INH_3_DPI_R1_A3  # one sample
# ============================================================

# ── Settings ────────────────────────────────────────────────
AF_THRESHOLD=0.01        # only apply variants above this AF
MIN_DEPTH=100           # mask positions with depth below this
# ────────────────────────────────────────────────────────────

# ── Paths (relative to Variant_discovery_pipeline/) ─────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../Variant_discovery_pipeline" && pwd)"
OUTDIR="$(cd "$SCRIPT_DIR/.." && pwd)/consensus"

REF="$PIPELINE_DIR/Input/Reference/viral_only.fasta"   # viral contig only
VIRAL_CONTIG="VEEV_INH"                                # contig name in ref
LOFREQ_DIR="$PIPELINE_DIR/LoFreq"
BAM_DIR="$PIPELINE_DIR/Input/BAMs"

mkdir -p "$OUTDIR"

# ── Conda (activate before set -e, conda can return non-zero) ──
eval "$(conda shell.bash hook)"
conda activate lofreq-env 2>/dev/null \
  || conda activate ivar_env 2>/dev/null \
  || { echo "ERROR: could not activate lofreq-env or ivar_env"; exit 1; }

set -euo pipefail

# Quick tool check
for tool in bcftools samtools bedtools tabix bgzip; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not found"; exit 1; }
done

# ── Decide which samples to process ────────────────────────
if [ $# -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for d in "$LOFREQ_DIR"/*/; do
        SAMPLES+=("$(basename "$d")")
    done
fi

echo "================================================================"
echo "  Consensus builder  (AF > ${AF_THRESHOLD}, mask depth < ${MIN_DEPTH})"
echo "  Samples: ${#SAMPLES[@]}"
echo "  Output : $OUTDIR"
echo "================================================================"
echo ""

# ── Process each sample ────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    echo "── $SAMPLE ──────────────────────────────────"
    VCF="$LOFREQ_DIR/$SAMPLE/variants.filtered.vcf.gz"
    OUT_FASTA="$OUTDIR/${SAMPLE}.consensus.fasta"

    # Prefer the small viral-only BAM from the LoFreq run (much faster)
    VIRAL_BAM="$LOFREQ_DIR/$SAMPLE/viral_only.bam"
    FULL_BAM="$BAM_DIR/${SAMPLE}.bam"
    if [ -f "$VIRAL_BAM" ]; then
        BAM="$VIRAL_BAM"
    else
        BAM="$FULL_BAM"
    fi

    # --- check inputs ---
    if [ ! -f "$VCF" ]; then
        echo "  SKIP: no filtered VCF found"
        echo ""
        continue
    fi
    if [ ! -f "$BAM" ]; then
        echo "  SKIP: no BAM found"
        echo ""
        continue
    fi

    TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    # 1. Filter to major-allele variants (AF > threshold)
    echo "  1) Filtering VCF  (AF > $AF_THRESHOLD) ..."
    RAW_NVAR=$(bcftools view -H "$VCF" | wc -l)
    bcftools view -i "INFO/AF > $AF_THRESHOLD" "$VCF" -Oz -o "$TMPDIR/filtered.vcf.gz"
    tabix -p vcf "$TMPDIR/filtered.vcf.gz"

    NVAR=$(bcftools view -H "$TMPDIR/filtered.vcf.gz" | wc -l)
    AF_DROPPED=$((RAW_NVAR - NVAR))
    echo "     → raw VCF variants: $RAW_NVAR"
    echo "     → $NVAR variants pass AF threshold"
    echo "     → $AF_DROPPED variants removed by AF filter"

    # 2. Mask low-depth positions to N  (viral contig only)
    echo "  2) Building depth mask  (depth < $MIN_DEPTH) ..."
    samtools depth -aa -r "$VIRAL_CONTIG" "$BAM" \
        | awk -v min="$MIN_DEPTH" '$3 < min {print $1"\t"($2-1)"\t"$2}' \
        > "$TMPDIR/lowdepth.bed"

    MASKED_SITES=$(wc -l < "$TMPDIR/lowdepth.bed")
    echo "     → $MASKED_SITES sites will be masked to N"

    if [ "$MASKED_SITES" -gt 0 ]; then
        bedtools maskfasta -fi "$REF" -bed "$TMPDIR/lowdepth.bed" -fo "$TMPDIR/ref_masked.fasta"
        CONSENSUS_REF="$TMPDIR/ref_masked.fasta"
    else
        CONSENSUS_REF="$REF"
    fi

    # 3. Build consensus
    echo "  3) Building consensus ..."
    bcftools consensus -f "$CONSENSUS_REF" "$TMPDIR/filtered.vcf.gz" \
        | sed "s/^>.*/>$SAMPLE/" \
        > "$OUT_FASTA"

    # 4. Sanity-check that filtered SNP variants are reflected in consensus
    echo "  4) Variant application sanity-check ..."
    bcftools query -f '%POS\t%REF\t%ALT\n' "$TMPDIR/filtered.vcf.gz" > "$TMPDIR/filtered_variants.tsv"
    python - "$CONSENSUS_REF" "$OUT_FASTA" "$TMPDIR/filtered_variants.tsv" <<'PY'
import sys

ref_fa, cons_fa, var_tsv = sys.argv[1], sys.argv[2], sys.argv[3]


def read_single_fasta(path):
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).upper()


consensus = read_single_fasta(cons_fa)
masked_ref = read_single_fasta(ref_fa)

total = 0
simple_snp = 0
simple_snp_applied = 0
simple_snp_mismatch = 0
simple_snp_oob = 0
indel_or_complex = 0
multiallelic = 0

with open(var_tsv, "r", encoding="utf-8") as handle:
    for row in handle:
        row = row.strip()
        if not row:
            continue
        pos_str, ref, alt = row.split("\t")
        total += 1

        if "," in alt:
            multiallelic += 1
            continue

        ref = ref.upper()
        alt = alt.upper()
        if len(ref) == 1 and len(alt) == 1:
            simple_snp += 1
            pos = int(pos_str)
            idx = pos - 1
            if idx < 0 or idx >= len(consensus):
                simple_snp_oob += 1
                continue
            if consensus[idx] == alt:
                simple_snp_applied += 1
            else:
                simple_snp_mismatch += 1
        else:
            indel_or_complex += 1

print(
    f"     → filtered variants by type: total={total}, simple_snp={simple_snp}, "
    f"indel_or_complex={indel_or_complex}, multiallelic={multiallelic}"
)
print(
    f"     → simple SNPs applied in consensus: {simple_snp_applied}/{simple_snp} "
    f"(mismatch={simple_snp_mismatch}, out_of_range={simple_snp_oob})"
)

if len(masked_ref) == len(consensus):
    # Ignore masked reference sites (N) when counting base differences.
    base_changes = sum(
        1 for r, c in zip(masked_ref, consensus)
        if r != "N" and c != "N" and r != c
    )
    print(f"     → consensus base changes vs masked reference: {base_changes}")
else:
    print(
        "     → consensus/reference lengths differ "
        f"({len(consensus)} vs {len(masked_ref)}); likely due to indel application"
    )
PY

    echo "  ✓  Saved: $OUT_FASTA"
    echo ""

    rm -rf "$TMPDIR"
done

echo "Done! Consensus files are in: $OUTDIR"
