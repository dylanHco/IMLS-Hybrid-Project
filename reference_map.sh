#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate extraction

# --------------------
# Configurable paths
# --------------------
REF="/home/dcl9541/Juglans_genomes/full_assembly/ncbi_dataset/data/GCA_041380795.1/Jhindsii_Rawlins_hap2_chrom.fasta"
READ_DIR="/data/labs/Fant/Cohen/Juglans/test_map"
OUT_DIR="gatk_output"
THREADS=18

# --------------------
# Step 1: Prepare reference (only once)
# --------------------
echo "Indexing reference genome..."
bwa index $REF
samtools faidx $REF
gatk CreateSequenceDictionary -R $REF -O ${REF%.fasta}.dict

# --------------------
# Step 2: Create output folders
# --------------------
mkdir -p $OUT_DIR/bam
mkdir -p $OUT_DIR/gvcf
mkdir -p $OUT_DIR/tmp

# --------------------
# Step 3: Process each sample in parallel (up to 5 at a time)
# --------------------
MAX_JOBS=5
CURRENT_JOBS=0

process_sample() {
    SAMPLE=$1
    R1="$READ_DIR/${SAMPLE}_R1_paired.fastq.gz"
    R2="$READ_DIR/${SAMPLE}_R2_paired.fastq.gz"
    echo "Processing $SAMPLE..."

    # Align
    bwa mem -t 6 $REF $R1 $R2 > $OUT_DIR/tmp/${SAMPLE}.sam

    # Add read group info (required for GATK)
    samtools addreplacerg \
        -r "ID:$SAMPLE" -r "SM:$SAMPLE" -r "PL:ILLUMINA" \
        -r "LB:lib1" -r "PU:unit1" \
        -o $OUT_DIR/tmp/${SAMPLE}.rg.sam \
        $OUT_DIR/tmp/${SAMPLE}.sam
    rm $OUT_DIR/tmp/${SAMPLE}.sam

    # Convert to BAM and sort
    samtools view -bS $OUT_DIR/tmp/${SAMPLE}.rg.sam | \
        samtools sort -@ 6 -o $OUT_DIR/bam/${SAMPLE}.sorted.bam
    rm $OUT_DIR/tmp/${SAMPLE}.rg.sam

    # Index BAM
    samtools index $OUT_DIR/bam/${SAMPLE}.sorted.bam

    # Mark duplicates using GATK
    gatk MarkDuplicatesSpark \
        -I $OUT_DIR/bam/${SAMPLE}.sorted.bam \
        -O $OUT_DIR/bam/${SAMPLE}.dedup.bam \
        --tmp-dir $OUT_DIR/tmp \
        --create-output-bam-index true \
        --verbosity ERROR

    # Call variants to GVCF
    gatk HaplotypeCaller \
        -R $REF \
        -I $OUT_DIR/bam/${SAMPLE}.dedup.bam \
        -O $OUT_DIR/gvcf/${SAMPLE}.g.vcf.gz \
        -ERC GVCF \
        --native-pair-hmm-threads 6
}

for R1 in $READ_DIR/*_R1_paired.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_paired.fastq.gz)

    # Start sample processing in the background
    process_sample "$SAMPLE" &

    # Control number of parallel jobs
    ((CURRENT_JOBS++))
    if [[ $CURRENT_JOBS -ge $MAX_JOBS ]]; then
        wait -n  # Wait for any job to finish
        ((CURRENT_JOBS--))
    fi
done

wait  # Wait for all background jobs to complete
# --------------------
# Step 4: Joint Genotyping
# --------------------
echo "Combining GVCFs and running joint genotyping..."

GVCF_LIST="$OUT_DIR/gvcf/gvcf_list.txt"
ls $OUT_DIR/gvcf/*.g.vcf.gz | awk '{print " -V "$1}' > $GVCF_LIST

gatk CombineGVCFs \
    -R $REF \
    $(cat $GVCF_LIST) \
    -O $OUT_DIR/combined.g.vcf.gz

gatk GenotypeGVCFs \
    -R $REF \
    -V $OUT_DIR/combined.g.vcf.gz \
    -O $OUT_DIR/variants.raw.vcf.gz

# --------------------
# Step 5: Variant Filtering
# --------------------
gatk VariantFiltration \
	-R $REF \
    -V $OUT_DIR/variants.raw.vcf.gz \
    -O $OUT_DIR/variants.filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_snp_filter"

echo "🎉 All done! Output VCF: $OUT_DIR/variants.filtered.vcf.gz"
