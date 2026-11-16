# 04_SV_merge_vcf_to_xlsx
## Local Environment (R)

Run `SV_merge_vcf_to_xlsx.R` to filter and merge the SV data, which consists of five steps:

Step 1: Filter out low-quality variants (depth < 5).

Step 2: Filter out variants that are present in the ancestor (D0).

Step 3: Take the union of manta and delly calls for each sample.

Step 4: Merge the data from multiple samples and output a single file named `all_SV_combined_table.xlsx`.

Step 5: Open the BAM files of the ancestor (D0) and the evolved samples locally in IGV to manually validate the variants, excluding false positives and those already present in the ancestor.
