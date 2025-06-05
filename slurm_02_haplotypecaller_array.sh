#!/bin/bash
#SBATCH --job-name=hapcall
#SBATCH --output=logs_ALL/hapcall_%A_%a.out
#SBATCH --error=logs_ALL/hapcall_%A_%a.err
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=0-267  # 268 samples
#SBATCH --partition=genomicslong
#SBATCH --account=b1042

set -euo pipefail
module load bwa/0.7.17 samtools gatk

REF="/scratch/dcl9541/Oak_genome/primary_scaffolds.fasta"
READ_DIR="/scratch/dcl9541/Oak_fastq"
OUT_DIR="/scratch/dcl9541/Oak_real_deal"
SAMPLES_FILE="$READ_DIR/samples.txt"
THREADS=8

mkdir -p $OUT_DIR/{bam,gvcf,tmp}

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLES_FILE)
R1="${READ_DIR}/${SAMPLE}_R1_paired.fastq.gz"
R2="${READ_DIR}/${SAMPLE}_R2_paired.fastq.gz"

# Alignment
bwa mem -t $THREADS $REF "$R1" "$R2" > $OUT_DIR/bam/${SAMPLE}.sam

# Add read groups
gatk AddOrReplaceReadGroups \
  -I $OUT_DIR/bam/${SAMPLE}.sam \
  -O $OUT_DIR/bam/${SAMPLE}.rg.bam \
  -RGID $SAMPLE -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $SAMPLE

# Sort and index
samtools sort -@ $THREADS -T $OUT_DIR/tmp/${SAMPLE} \
  -o $OUT_DIR/bam/${SAMPLE}.sorted.bam \
  $OUT_DIR/bam/${SAMPLE}.rg.bam
samtools index $OUT_DIR/bam/${SAMPLE}.sorted.bam

# Mark duplicates
gatk MarkDuplicatesSpark \
  -I $OUT_DIR/bam/${SAMPLE}.sorted.bam \
  -O $OUT_DIR/bam/${SAMPLE}.dedup.bam \
  --tmp-dir $OUT_DIR/tmp \
  --create-output-bam-index true

# Call variants
gatk HaplotypeCaller \
  -R $REF \
  -I $OUT_DIR/bam/${SAMPLE}.dedup.bam \
  -O $OUT_DIR/gvcf/${SAMPLE}.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads $THREADS
