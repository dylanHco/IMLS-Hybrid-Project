#!/bin/bash
#SBATCH --job-name=jointgeno_%A_%a
#SBATCH --output=logs/jointgeno_%A_%a.out
#SBATCH --error=logs/jointgeno_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --partition=genomicslong
#SBATCH --account=b1042
#SBATCH --array=1-<N>   # replace <N> with number of scaffolds

### Create file with scaffold names
#cut -f1 /scratch/dcl9541/Oak_genome/primary_scaffolds.fasta.fai > /scratch/dcl9541/Oak_real_deal/scaffold_list.txt

set -euo pipefail
module load gatk samtools

REF="/scratch/dcl9541/Oak_genome/primary_scaffolds.fasta"
OUT_DIR="/scratch/dcl9541/Oak_real_deal"
GVCF_DIR="$OUT_DIR/gvcf"
SCAF_LIST="$OUT_DIR/scaffold_list.txt"

# Get scaffold name
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SCAF_LIST)

# GVCF list
GVCF_LIST="$OUT_DIR/gvcf_list.args"
ls $GVCF_DIR/*.g.vcf.gz | awk '{print "-V",$1}' > $GVCF_LIST

# Output workspace for this scaffold
WORKDIR="$OUT_DIR/genomicsdb_$SCAF"

gatk GenomicsDBImport \
  --genomicsdb-workspace-path "$WORKDIR" \
  --batch-size 50 \
  -L "$SCAF" \
  $(cat "$GVCF_LIST")

gatk GenotypeGVCFs \
  -R "$REF" \
  -V gendb://"$WORKDIR" \
  -O "$OUT_DIR/variants_raw_$SCAF.vcf.gz"
