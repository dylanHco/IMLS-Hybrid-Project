#!/bin/bash
set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate extraction

# Input and output paths
VCF_IN="gatk_output/variants.filtered.vcf.gz"
WORK_DIR="gatk_output"
LOG="$WORK_DIR/filtering_log.txt"

echo "📝 VCF filtering log" > "$LOG"
echo "Input VCF: $VCF_IN" >> "$LOG"
echo "Start time: $(date)" >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 0: Initial stats
# --------------------
bcftools stats "$VCF_IN" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 1: Biallelic SNPs
# --------------------
echo "Step 1: Filtering to biallelic SNPs..." >> "$LOG"
bcftools view -m2 -M2 -v snps "$VCF_IN" -Oz -o "$WORK_DIR/step1.vcf.gz"
tabix -p vcf "$WORK_DIR/step1.vcf.gz"
bcftools stats "$WORK_DIR/step1.vcf.gz" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 2: QUAL ≥ 30 and per-sample DP 5–60
# --------------------
echo "Step 2: Filtering by QUAL and per-sample DP..." >> "$LOG"
bcftools filter -e 'QUAL<30 || FORMAT/DP<5 || FORMAT/DP>60' \
  -Oz -o "$WORK_DIR/step2.vcf.gz" "$WORK_DIR/step1.vcf.gz"
tabix -p vcf "$WORK_DIR/step2.vcf.gz"
bcftools stats "$WORK_DIR/step2.vcf.gz" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 3: Filter by mean site depth
# --------------------
echo "Step 3: Filtering by mean site depth..." >> "$LOG"
vcftools --gzvcf "$WORK_DIR/step2.vcf.gz" --site-mean-depth --out "$WORK_DIR/site_depth"
awk 'NR > 1 && $3 >= 5 && $3 <= 60 {print $1"\t"$2}' "$WORK_DIR/site_depth.ldepth.mean" > "$WORK_DIR/keep_sites_step3.txt"

echo "Sites kept: $(wc -l < "$WORK_DIR/keep_sites_step3.txt")" >> "$LOG"
head "$WORK_DIR/keep_sites_step3.txt" >> "$LOG"
echo "" >> "$LOG"

bcftools view -R "$WORK_DIR/keep_sites_step3.txt" "$WORK_DIR/step2.vcf.gz" -Oz -o "$WORK_DIR/step3.vcf.gz"
tabix -p vcf "$WORK_DIR/step3.vcf.gz"
bcftools stats "$WORK_DIR/step3.vcf.gz" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 4: Remove variants with >40% missing data
# --------------------
echo "Step 4: Filtering variants with >40% missing data..." >> "$LOG"
NUM_SAMPLES=$(bcftools query -l "$WORK_DIR/step3.vcf.gz" | wc -l)

bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$WORK_DIR/step3.vcf.gz" | \
awk -v n=$NUM_SAMPLES '{
  m=0;
  for (i=3; i<=NF; i++) if ($i == "./.") m++;
  if ((NF - 2 - m) >= n * 0.6) print $1"\t"$2
}' > "$WORK_DIR/keep_sites_step4.txt"

if [ ! -s "$WORK_DIR/keep_sites_step4.txt" ]; then
    echo "Error: No variants passed the missingness threshold." >> "$LOG"
    exit 1
fi

bcftools view -T "$WORK_DIR/keep_sites_step4.txt" "$WORK_DIR/step3.vcf.gz" \
  -Oz -o "$WORK_DIR/step4.vcf.gz"
tabix -p vcf "$WORK_DIR/step4.vcf.gz"
bcftools stats "$WORK_DIR/step4.vcf.gz" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# --------------------
# Step 5: Remove samples with >40% missing data
# --------------------
echo "Step 5: Removing samples with >40% missing data (vcftools)..." >> "$LOG"
VCFTOOLS_TMP="$WORK_DIR/vcftools_tmp"
mkdir -p "$VCFTOOLS_TMP"

vcftools --gzvcf "$WORK_DIR/step4.vcf.gz" --missing-indv --out "$VCFTOOLS_TMP/missing"
awk -v threshold=0.4 'NR > 1 && $5 <= threshold { print $1 }' "$VCFTOOLS_TMP/missing.imiss" > "$WORK_DIR/samples_to_keep_step5.txt"

TOTAL_SAMPLES=$(($(wc -l < "$VCFTOOLS_TMP/missing.imiss") - 1))
KEPT_SAMPLES=$(wc -l < "$WORK_DIR/samples_to_keep_step5.txt")
REMOVED_SAMPLES=$((TOTAL_SAMPLES - KEPT_SAMPLES))

echo "Total samples: $TOTAL_SAMPLES" >> "$LOG"
echo "Samples kept: $KEPT_SAMPLES" >> "$LOG"
echo "Samples removed: $REMOVED_SAMPLES" >> "$LOG"

if [ "$KEPT_SAMPLES" -eq 0 ]; then
    echo "Error: No samples left after filtering." >> "$LOG"
    exit 1
fi

bcftools view -S "$WORK_DIR/samples_to_keep_step5.txt" "$WORK_DIR/step4.vcf.gz" \
  -Oz -o "$WORK_DIR/step5.vcf.gz"
tabix -p vcf "$WORK_DIR/step5.vcf.gz"
bcftools stats "$WORK_DIR/step5.vcf.gz" | grep "^SN" | cut -f 2- >> "$LOG"
echo "" >> "$LOG"

# Final rename
cp "$WORK_DIR/step5.vcf.gz" "$WORK_DIR/variants.filtered.final.vcf.gz"
cp "$WORK_DIR/step5.vcf.gz.tbi" "$WORK_DIR/variants.filtered.final.vcf.gz.tbi"

echo "✅ Final VCF: $WORK_DIR/variants.filtered.final.vcf.gz" >> "$LOG"
echo "End time: $(date)" >> "$LOG"
echo "🎉 Filtering complete." >> "$LOG"
