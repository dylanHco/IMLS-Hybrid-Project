#!/bin/bash
#SBATCH --job-name=jointgeno
#SBATCH --output=logs/jointgeno_%j.out # change this
#SBATCH --error=logs/jointgeno_%j.err # change this
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=150G
#SBATCH --partition=genomicslong
#SBATCH --account=b1042

set -euo pipefail
module load gatk samtools

REF="/scratch/dcl9541/Oak_genome/primary_scaffolds.fasta"
OUT_DIR="/scratch/dcl9541/Oak_real_deal"
GVCF_DIR="$OUT_DIR/gvcf"
INTERVALS="$OUT_DIR/intervals.list"

mkdir -p logs

# Create intervals list
cut -f1 "${REF}.fai" > "$INTERVALS"

# List of GVCFs
GVCF_LIST="$OUT_DIR/gvcf_list.args"
ls $GVCF_DIR/*.g.vcf.gz | awk '{print "-V",$1}' > $GVCF_LIST

# GenomicsDBImport
gatk GenomicsDBImport \
  --genomicsdb-workspace-path $OUT_DIR/genomicsdb \
  --batch-size 50 \
  -L "$INTERVALS" \
  $(cat "$GVCF_LIST")

# Genotype
gatk GenotypeGVCFs \
  -R "$REF" \
  -V gendb://$OUT_DIR/genomicsdb \
  -O "$OUT_DIR/variants.raw.vcf.gz"

# Filter
gatk VariantFiltration \
  -R "$REF" \
  -V "$OUT_DIR/variants.raw.vcf.gz" \
  -O "$OUT_DIR/variants.filtered.vcf.gz" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_snp_filter"
