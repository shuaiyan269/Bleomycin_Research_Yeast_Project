
library(vcfR)
library(tidyverse)   
library(stringr)     
library(magrittr)   
library(tidyr)
library(ggplot2)

id <-"MA"
vcf_path <- paste0("E:/Bleomycin_Research_Yeast_Project/03_SNV_INDEL_vcf_to_xlsx/input/", id, ".final_filtered.vcf")
vcf <- read.vcfR(vcf_path)
vcf_df <- as.data.frame(vcf@fix)
vcf_gt <- as.data.frame(vcf@gt)
vcf_full <- cbind(vcf_df, vcf_gt[, -1]) 

csq_example <- vcf_full$INFO[1]  
csq_value <- str_extract(csq_example, "CSQ=([^;]+)")
csq_parts <- strsplit(sub("CSQ=", "", csq_value), "\\|")[[1]]
for (i in 1:length(csq_parts)) {
  cat(i, ": ", csq_parts[i], "\n")
}

vcf_clean <- vcf_full %>%
  mutate(
    CSQ = sapply(INFO, function(x) {
      csq <- str_extract(x, "CSQ=([^;]+)")
      ifelse(!is.na(csq), sub("CSQ=", "", csq), NA)
    })
  ) %>%
  separate(
    CSQ, 
    into = c(
      "Allele", "Consequence", "Impact", "Gene_ID", "Gene_Name",
      "Feature_type", "Feature", "Biotype", "Exon", "Intron",
      "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
      "Amino_acids", "Codons", "Existing_variation", "Distance",
      "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID",
      "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT",
      "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO", "SIFT",
      "POLYPHEN", "DOMAINS", "GMAF", "AFR_MAF", "AMR_MAF",
      "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF",
      "EXON_INTRON_RANK", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
      "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
      "Annotation_Source", "Annotation_File"
    ),
    sep = "\\|", 
    remove = FALSE,
    fill = "right"
  ) %>%
  select(
    CHROM, POS, REF, ALT, QUAL, Impact, STRAND,
    Allele, Gene = Gene_Name, Gene_ID, Feature, Feature_type, Consequence, 
    cDNA_position, CDS_position, Protein_position, 
    Amino_acids, Codons, Existing_variation,
    Biotype, Distance,
    everything()
  ) %>%
  select(
    CHROM, POS, REF, ALT, QUAL, Impact, STRAND,
    Allele, Gene, Gene_ID, Feature, Feature_type, Consequence, 
    cDNA_position, CDS_position, Protein_position, 
    Amino_acids, Codons, Existing_variation,
    Biotype, Distance,INFO,
    starts_with("D0_"),  
    starts_with("D90_")  
  )

snv_df <- vcf_clean %>%
  filter(grepl("VT=SNP", INFO) | grepl("VT=SNP;", INFO))

indel_df <- vcf_clean %>%
  filter(grepl("VT=INDEL", INFO) | grepl("VT=INDEL;", INFO))

write.csv(snv_df, paste0("E:/Bleomycin_Research_Yeast_Project/03_SNV_INDEL_vcf_to_xlsx/output/", id, "_anno_snv.csv"), row.names = FALSE)
write.csv(indel_df, paste0("E:/Bleomycin_Research_Yeast_Project/03_SNV_INDEL_vcf_to_xlsx/output/", id, "_anno_indel.csv"), row.names = FALSE)


# Local filtering: remove mutations present in the ancestor (D0_XX)


