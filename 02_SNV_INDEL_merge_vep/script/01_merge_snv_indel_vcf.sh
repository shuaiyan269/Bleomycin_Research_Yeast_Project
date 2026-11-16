#!/bin/bash
# the freebayes and gatk4 of same sample merge with filtering

source ./../../../env.sh

id="BY" #BY, WE, WA, NA, SA, MA
INPUT_DIR="./../../01_variant_calling/output/02_SNV_INDEL/${id}" ###id
REFERENCE="./../input/Reference_Genomes/${id}.tidy.lite.fa"
OUTPUT_DIR="./../output/merge_snv_indel"
BCFTOOLS="./../../../build/bcftools-1.17/bcftools"

# make output directory 
mkdir -p "$OUTPUT_DIR"

# get all samples name
SAMPLES=$(find "$INPUT_DIR" -mindepth 2 -name "*.freebayes.filtered.vcf*" | 
          sed 's/.*\/\([^/]*\)\/.*\.freebayes.*/\1/' | sort -u)

for SAMPLE in $SAMPLES; do
  echo "====== processing $SAMPLE ======"
  
  # define path
  GATK_VCF=$(find "$INPUT_DIR" -path "*/${SAMPLE}/${SAMPLE}.gatk4.filtered.vcf*")
  FREEBAYES_VCF=$(find "$INPUT_DIR" -path "*/${SAMPLE}/${SAMPLE}.freebayes.filtered.vcf*")
  SAMPLE_OUTDIR="${OUTPUT_DIR}/${id}/${SAMPLE}"
  mkdir -p "$SAMPLE_OUTDIR"
  OUTPUT="${SAMPLE_OUTDIR}/${SAMPLE}_merged.vcf.gz"
  
  # check
  if [[ ! -f "$GATK_VCF" || ! -f "$FREEBAYES_VCF" ]]; then
    echo "Warning: Input files missing for sample $SAMPLE, skipping"
    continue
  fi

  # Normalize VCFs
  $BCFTOOLS norm -m-any -f "$REFERENCE" "$GATK_VCF" -Oz -o gatk_norm.vcf.gz
  $BCFTOOLS norm -m-any -f "$REFERENCE" "$FREEBAYES_VCF" -Oz -o freebayes_norm.vcf.gz
  $BCFTOOLS index gatk_norm.vcf.gz
  $BCFTOOLS index freebayes_norm.vcf.gz

  # Apply filters based on the image criteria
  # GATK4 filters
  $BCFTOOLS view -i '(FS<60 || FS=".") && 
  (MQ>40 || MQ=".") && 
  (MQRankSum>-3 || MQRankSum=".") && 
  (QD>3 || QD=".") && 
  (ReadPosRankSum>-2 || ReadPosRankSum=".") && 
  (SOR<3 || SOR=".")' \
  gatk_norm.vcf.gz -Oz -o gatk_filtered.vcf.gz
  $BCFTOOLS index gatk_filtered.vcf.gz

  # freebayes filters
  $BCFTOOLS view -i '(EPP>3 || EPP=".") && 
  (SAP>3) && 
  (SRP>3 || INFO/RO==0)' \
  freebayes_norm.vcf.gz -Oz -o freebayes_filtered.vcf.gz
  $BCFTOOLS index freebayes_filtered.vcf.gz

  # Extract chrom/pos/ref/alt for filtered variants
  $BCFTOOLS query -f '%CHROM\t%POS\t%REF\t%ALT\n' gatk_filtered.vcf.gz | 
    awk '{split($4,a,","); for(i in a) print $1"\t"$2"\t"$3"\t"a[i]}' | sort > gatk_simple.txt
  $BCFTOOLS query -f '%CHROM\t%POS\t%REF\t%ALT\n' freebayes_filtered.vcf.gz | 
    awk '{split($4,a,","); for(i in a) print $1"\t"$2"\t"$3"\t"a[i]}' | sort > freebayes_simple.txt

  # Find common variants
  comm -12 gatk_simple.txt freebayes_simple.txt > common_simple.txt

  if [ ! -s common_simple.txt ]; then
    echo "[$SAMPLE] empty"
    touch "${SAMPLE_OUTDIR}/merged.empty"
    continue
  fi

  # Extract GATK4 info for common variants
  $BCFTOOLS view gatk_filtered.vcf.gz | 
    awk -F'\t' 'BEGIN {while(getline <"common_simple.txt") {key=$1"\t"$2"\t"$3"\t"$4; a[key]=1}} 
                /^#/ {print; next} 
                {key=$1"\t"$2"\t"$4"\t"$5; if(a[key]) print}' |
    $BCFTOOLS view -Oz -o "$OUTPUT"

  $BCFTOOLS index "$OUTPUT"
  rm gatk_simple.txt freebayes_simple.txt common_simple.txt \
     gatk_norm.vcf.gz* freebayes_norm.vcf.gz* \
     gatk_filtered.vcf.gz* freebayes_filtered.vcf.gz*

  echo "[$SAMPLE] finish → $OUTPUT"
  echo "variation number: $($BCFTOOLS stats "$OUTPUT" | grep "number of records" | awk '{print $NF}')"
  echo "------------------------------------"
done

echo "All samples have been processed! See the results: $OUTPUT_DIR"
