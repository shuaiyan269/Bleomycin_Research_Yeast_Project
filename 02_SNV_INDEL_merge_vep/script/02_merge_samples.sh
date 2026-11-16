#!/bin/bash
# the freebayes and gatk4 of same sample merge with filtering

source ./../../../env.sh

id="NA" #BY, WE, WA, NA, SA, MA
INPUT_DIR="./../output/merge_snv_indel/${id}" ###id
OUTPUT_DIR="./../output/merge_sample"
BCFTOOLS="./../../../build/bcftools-1.17/bcftools"
LIST_FILE="${OUTPUT_DIR}/${id}_vcf_list.txt"
OUT_VCF_GZ="${OUTPUT_DIR}/${id}_snv_indel_merged.vcf.gz"
OUT_VCF="${OUTPUT_DIR}/${id}_snv_indel_merged.vcf"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "====== Start merging samples for background: ${id} ======"

# 1. Generate file list
find "$INPUT_DIR" -type f -name "*.vcf.gz" | sort > "$LIST_FILE"

# Check if any files were found
if [ ! -s "$LIST_FILE" ]; then
    echo "Error: No .vcf.gz files found in $INPUT_DIR"
    exit 1
fi

# 2. Perform merge and index
echo "Merging VCFs..."
$BCFTOOLS merge $(cat "$LIST_FILE") -Oz -o "$OUT_VCF_GZ"
$BCFTOOLS index "$OUT_VCF_GZ"

# 3. Decompress to plain text VCF
echo "Decompressing to plain text VCF..."
$BCFTOOLS view "$OUT_VCF_GZ" -o "$OUT_VCF"

# 4. Clean up list file
rm "$LIST_FILE"

echo "====== Done! Result file: $OUT_VCF ======"
