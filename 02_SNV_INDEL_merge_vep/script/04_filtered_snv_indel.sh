#!/bin/bash

# input
id="MA" #BY, WE, WA, NA, SA, MA
INPUT_VCF="./../output/vep/${id}/all_${id}_snv_indel_vepanno.vcf"
OUTPUT_VCF="./../output/vep/${id}/${id}.final_filtered.vcf"
BCFTOOLS="./../../../build/bcftools-1.17/bcftools"
TEMP_FILE1="temp_filtered.vcf"
TEMP_FILE2="temp_marked.vcf"

# First: filter out mutations on repeats and subtelomeres
echo "Step 1: Filtering repeats and subtelomere regions..."
$BCFTOOLS view \
    -i 'INFO/CSQ !~ "repeats_" && INFO/CSQ !~ "subtelomere_"' \
    "$INPUT_VCF" -o "$TEMP_FILE1"

# Second: filter out the mutaions of total depth < 10 or alt depth < 3
echo "Step 2: Performing sample-level depth filtering..."
awk -F'\t' '
BEGIN {OFS="\t"}
/^#/ {print; next}
{
  split($9, fmt, ":");
  for (i=10; i<=NF; i++) {
    split($i, vals, ":");
    dp_idx = ad_idx = 0;
    for (j=1; j<=length(fmt); j++) {
      if (fmt[j] == "DP") dp_idx = j;
      if (fmt[j] == "AD") ad_idx = j;
    }
    if (dp_idx && ad_idx && vals[dp_idx] != "." && vals[ad_idx] != ".") {
      split(vals[ad_idx], ad, ",");
      # check total depth
      if (vals[dp_idx] + 0 < 10) {
        $i = "./.:.:.:.:.";
        continue;
      }
      # check alt depth
      has_sufficient_alt = 0;
      for (k=2; k<=length(ad); k++) {  # start from second, skip ref count
        if (ad[k] + 0 >= 3) {
          has_sufficient_alt = 1;
          break;
        }
      }
      if (!has_sufficient_alt) {
        $i = "./.:.:.:.:.";
      }
    }
  }
  print
}' "$TEMP_FILE1" > "$TEMP_FILE2"

# Third: filter the mutations if all samples = 0
echo "Step 3: Removing sites where all samples are missing..."
$BCFTOOLS view \
    -i 'COUNT(GT!="mis")>0' \
    "$TEMP_FILE2" -o "$OUTPUT_VCF"

# clean
echo "Cleaning up temporary files..."
rm -f "$TEMP_FILE1" "$TEMP_FILE2"

echo "Filtering completed! Final output: $OUTPUT_VCF"

# statistic information
echo -e "\n=== Filtering Statistics ==="
original_count=$(grep -v "^#" "$INPUT_VCF" | wc -l)
final_count=$(grep -v "^#" "$OUTPUT_VCF" | wc -l)
filtered_count=$((original_count - final_count))

echo "Original variants: $original_count"
echo "Final variants: $final_count"
echo "Filtered variants: $filtered_count"