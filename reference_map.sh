#!/bin/bash
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate extraction

# --------------------
# Configurable paths
# --------------------
REF="/home/dcl9541/Juglans_genomes/full_assembly/ncbi_dataset/data/GCA_041380795.1/Jhindsii_Rawlins_hap2_chrom.fasta"
READ_DIR="/data/labs/Fant/Cohen/Juglans/after_trim/test_these"
OUT_DIR="test_gatk_output"
THREADS=60              # Total threads available on your machine
SAMPLE_THREADS=6        # Threads per sample
PARALLEL_JOBS=10        # Number of concurrent jobs

# --------------------
# Step 1: Prepare reference (only once)
# --------------------
#echo "📦 Indexing reference genome..."
#bwa index $REF
#samtools faidx $REF
#gatk CreateSequenceDictionary -R $REF -O ${REF%.fasta}.dict

# --------------------
# Step 2: Create output folders
# --------------------
mkdir -p $OUT_DIR/{bam,gvcf,tmp}

# --------------------
# Step 3: Define sample processing function
# --------------------
process_sample() {
    SAMPLE=$1
    R1="${READ_DIR}/${SAMPLE}_R1_paired.fastq.gz"
    R2="${READ_DIR}/${SAMPLE}_R2_paired.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "❌ Missing files for $SAMPLE"
        return 1
    fi

    echo "🔄 Processing $SAMPLE..."

    bwa mem -t $SAMPLE_THREADS $REF "$R1" "$R2" | \
    samtools addreplacerg -r "ID:$SAMPLE" -r "SM:$SAMPLE" -r "PL:ILLUMINA" \
                          -r "LB:lib1" -r "PU:unit1" - | \
    samtools sort -@ $SAMPLE_THREADS -o $OUT_DIR/bam/${SAMPLE}.sorted.bam

    samtools index $OUT_DIR/bam/${SAMPLE}.sorted.bam

    gatk MarkDuplicatesSpark \
        -I $OUT_DIR/bam/${SAMPLE}.sorted.bam \
        -O $OUT_DIR/bam/${SAMPLE}.dedup.bam \
        --tmp-dir $OUT_DIR/tmp \
        --create-output-bam-index true \
        --verbosity ERROR

    gatk HaplotypeCaller \
        -R $REF \
        -I $OUT_DIR/bam/${SAMPLE}.dedup.bam \
        -O $OUT_DIR/gvcf/${SAMPLE}.g.vcf.gz \
        -ERC GVCF \
        --native-pair-hmm-threads $SAMPLE_THREADS
}
export -f process_sample
export READ_DIR OUT_DIR REF SAMPLE_THREADS

# --------------------
# Step 4: Run all samples in parallel
# --------------------
find "$READ_DIR" -name '*_R1_paired.fastq.gz' | while read -r R1; do
    basename "$R1" _R1_paired.fastq.gz
done | sort -u | \
xargs -n 1 -P $PARALLEL_JOBS -I {} bash -c 'process_sample "$@"' _ {}

# --------------------
# Step 5: Joint Genotyping
# --------------------
echo "🧬 Joint genotyping...

GVCF_LIST="$OUT_DIR/gvcf/gvcf_list.args"
ls $OUT_DIR/gvcf/*.g.vcf.gz | awk '{print "-V",$1}' > $GVCF_LIST

gatk GenomicsDBImport \
    --genomicsdb-workspace-path $OUT_DIR/genomicsdb \
    --batch-size 50 \
    $(cat $GVCF_LIST)

gatk GenotypeGVCFs \
    -R $REF \
    -V gendb://$OUT_DIR/genomicsdb \
    -O $OUT_DIR/variants.raw.vcf.gz



# --------------------
# Step 6: Variant Filtering
# --------------------
gatk VariantFiltration \
    -R $REF \
    -V $OUT_DIR/variants.raw.vcf.gz \
    -O $OUT_DIR/variants.filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_snp_filter
    
    echo "🎉 Done! Filtered VCF: $OUT_DIR/variants.filtered.vcf.gz"
