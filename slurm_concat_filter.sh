#!/bin/bash
#SBATCH --job-name=vcf_concat_filter
#SBATCH --output=logs/vcf_concat_filter.out
#SBATCH --error=logs/vcf_concat_filter.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=genomics
#SBATCH --account=b1042

set -euo pipefail

# Load modules if needed
module load gatk bcftools

# ---- Config ----
REF="/scratch/dcl9541/Oak_genome/primary_scaffolds.fasta"
OUT_DIR="/scratch/dcl9541/Oak_real_deal"
SCAF_LIST="$OUT_DIR/scaffold_list.txt"

RAW_VCF="$OUT_DIR/variants.raw.vcf.gz"
FILTERED_VCF="$OUT_DIR/variants.filtered.vcf.gz"

# ---- Concatenate VCFs ----
echo "Concatenating VCFs..."
VCF_FILES=$(awk -v dir="$OUT_DIR" '{print dir "/variants_raw_" $1 ".vcf.gz"}' $SCAF_LIST)

bcftools concat -Oz -o "$RAW_VCF" $VCF_FILES

# ---- Index the concatenated VCF ----
echo "Indexing concatenated VCF..."
bcftools index "$RAW_VCF"

# ---- Apply Variant Filtering ----
echo "Filtering SNPs..."
gatk VariantFiltration \
  -R "$REF" \
  -V "$RAW_VCF" \
  -O "$FILTERED_VCF" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_snp_filter"

echo "Done: Filtered VCF at $FILTERED_VCF"
