
library(tidyverse)   
library(stringr)     
library(magrittr)   
library(tidyr)
library(ggplot2)
library(writexl)
library(maftools)
library(org.Sc.sgd.db)
library(AnnotationDbi)


backgrounds <- c("BY", "WE", "WA", "NA", "SA", "MA")
base_dir <- "E:/Bleomycin_Research_Yeast_Project/03_SNV_INDEL_vcf_to_xlsx/output/filter_ancestor/"


all_combined_long <- map_dfr(backgrounds, function(id) {
  
  indel_path <- paste0(base_dir, id, "_anno_indel_filtered.csv")
  snv_path <- paste0(base_dir, id, "_anno_snv_filtered.csv")
  
  indel <- read_csv(indel_path)
  snv<- read_csv(snv_path)
  
  snv$cDNA_position <- as.character(snv$cDNA_position)
  snv$CDS_position <- as.character(snv$CDS_position)
  snv$Protein_position <- as.character(snv$Protein_position)
  indel$cDNA_position <- as.character(indel$cDNA_position)
  indel$CDS_position <- as.character(indel$CDS_position)
  indel$Protein_position <- as.character(indel$Protein_position)
  
  snv$variant_type <- "SNV"
  indel$variant_type <- "InDel"
  
  combined_df <- bind_rows(snv, indel)
  
  sample_names <- colnames(combined_df)[23:(ncol(combined_df) - 1)]
  
  combined_long <- combined_df %>%
    pivot_longer(
    cols = all_of(sample_names),
    names_to = "Sample",
    values_to = "Genotype_Info") %>%
    mutate(
    Genotype = str_split(Genotype_Info, ":", simplify = TRUE)[,1],
    Has_Variant = !(Genotype %in% c("./.", "0/0", "0|0")) & !is.na(Genotype)) %>%
    filter(Has_Variant) %>%
    separate(
      col = Sample,           
      into = c("dummy", "strain", "drug", "replicate"), 
      sep = "_",              
      remove = FALSE          
      ) 
  combined_long <- combined_long %>% dplyr::select(-Has_Variant)
  return(combined_long)
})

all_combined_long$gene_name <- mapIds(org.Sc.sgd.db,
                                 keys = all_combined_long$Gene_ID,
                                 column = "GENENAME",  
                                 keytype = "ORF",
                                 multiVals = "first")  
list_cols <- sapply(all_combined_long, is.list)
names(all_combined_long)[list_cols]
head(all_combined_long[list_cols])
all_combined_long <- all_combined_long %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ";"))))

write_xlsx(all_combined_long,"E:/Bleomycin_Research_Yeast_Project/03_SNV_INDEL_vcf_to_xlsx/output/combined_long.xlsx")

#The BAM files of the ancestor and the samples were visualized in IGV to manually validate each variant. False positives and variants already present in the ancestor were removed, yielding the final Table S3.


