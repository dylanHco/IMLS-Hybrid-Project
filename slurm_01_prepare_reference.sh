#!/bin/bash
#SBATCH --job-name=prepare_ref
#SBATCH --output=logs/prepare_ref_%j.out
#SBATCH --error=logs/prepare_ref_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=genomicslong
#SBATCH --account=b1042

set -euo pipefail
module load bwa/0.7.17 samtools gatk

REF="/scratch/dcl9541/Oak_genome/primary_scaffolds.fasta"

bwa index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fasta}.dict"
