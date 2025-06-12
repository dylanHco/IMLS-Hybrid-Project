#!/bin/bash
set -euo pipefail

# Input VCF
VCF=$1
BASENAME=$(basename "$VCF" .vcf.gz)

# Step 1: Initial site-level filtering
vcftools --gzvcf "$VCF" \
  --remove-indels \
  --min-alleles 2 \
  --max-alleles 2 \
  --max-missing 0.2 \
  --minQ 30 \
  --recode --recode-INFO-all \
  --out ${BASENAME}.step1

# Step 2: Apply depth and allele count filters
vcftools --vcf ${BASENAME}.step1.recode.vcf \
  --mac 2 \
  --minDP 5 \
  --min-meanDP 5 \
  --max-meanDP 60 \
  --maxDP 60 \
  --recode --recode-INFO-all \
  --out ${BASENAME}.step2

# Step 3: Prepare filtered outputs with different missing data thresholds
for MISSING in 0.5 0.6 0.7; do
  LABEL=${MISSING/./}  # Convert 0.5 → 05 for filenames

  # Compute missingness per individual
  vcftools --vcf ${BASENAME}.step2.recode.vcf \
    --missing-indv \
    --out ${BASENAME}.miss${LABEL}

  # Remove individuals with >70% missing data
  awk '$5 > 0.7 {print $1}' ${BASENAME}.miss${LABEL}.imiss > ${BASENAME}.remove70_${LABEL}.txt

  # Apply final missingness filter
  vcftools --vcf ${BASENAME}.step2.recode.vcf \
    --max-missing $MISSING \
    --remove ${BASENAME}.remove70_${LABEL}.txt \
    --recode --recode-INFO-all \
    --out ${BASENAME}.final_miss${LABEL}
done

echo "Filtering complete. Outputs:"
echo " - ${BASENAME}.final_miss05.recode.vcf"
echo " - ${BASENAME}.final_miss06.recode.vcf"
echo " - ${BASENAME}.final_miss07.recode.vcf"
