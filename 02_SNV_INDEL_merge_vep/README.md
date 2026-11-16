# 02 SNV INDEL Merge and VEP Annotation
## Server Environment (Bash)

Step 1: Run `01_merge_snv_indel_vcf.sh`
For each sample, filter the FreeBayes and GATK4 results according to the specified criteria, take the intersection (consensus) of the two, and output the GATK4-based records.
Run this sequentially for each genetic background by modifying the `id` variable to BY, WE, WA, NA, SA, and MA.

Step 2: Run `02_merge_samples.sh`
Merge the VCF files of multiple samples within the same genetic background into a single multi-sample VCF file.
Run this sequentially for each genetic background by modifying the id variable to BY, WE, WA, NA, SA, and MA.

Step 3: Run `03_VEP.sh`
Annotate the merged VCF files using Ensembl VEP.
Run this sequentially for each genetic background by modifying the `id` variable to BY, WE, WA, NA, SA, and MA.

Step 4: Run `04_filtered_snv_indel.sh`
Filter out mutations located in repeat and subtelomere regions, as well as mutations with a total depth < 10 or an alternative allele depth < 3.
Run this sequentially for each genetic background by modifying the `id` variable to BY, WE, WA, NA, SA, and MA.

Final Output
At the end of this step, one final filtered VCF file is generated for each genetic background. These files can then be downloaded to the local machine for downstream analysis.