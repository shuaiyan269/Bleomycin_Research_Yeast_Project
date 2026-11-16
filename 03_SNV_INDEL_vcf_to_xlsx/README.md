# 03_SNV_INDEL_vcf_to_xlsx
## Local Environment (R)

Step 1: For each genetic background, run `01_vcf_to_excel.R` to convert the VCF files into CSV tables. The `id` parameters are BY, WE, WA, NA, SA, and MA sequentially.

Step 2: Open the CSV files locally and filter out the variant rows that are present in the ancestor (D0).

Step 3: Run `02_merge_all_backgrounds.R` to transform the filtered files into a long format and merge the six genetic backgrounds together to generate `combined_long.xlsx`.

Step 4: Open the BAM files of the ancestor and samples in IGV locally to manually validate the variants in `combined_long.xlsx`. Filter out false positives and variants present in the ancestor to obtain the final `Table.S3`.
