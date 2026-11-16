library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)

combined_data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
combined_data <- combined_data %>%
  filter(Sample_type == "population", Drug == "BLM")

indel_processed <- combined_data %>%
  filter(Variant_type == "INDEL") %>%
  mutate(
    ref_len = nchar(REF),
    alt_len = nchar(ALT),
    
    len_diff = alt_len - ref_len
  ) %>%
  mutate(
    indel_class = ifelse(len_diff > 0, "insertion", "deletion"),
    indel_size = abs(len_diff),
    size_bucket = ifelse(indel_size >= 5, "5plus", paste0(indel_size, "bp")),
    column_title = paste(indel_class, size_bucket, sep = "_")
  )

indel_summary_table <- indel_processed %>%
  group_by(Sample, column_title) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = column_title, values_from = count, values_fill = 0)

write.csv(indel_summary_table, 
          file = "E:/Bleomycin_Research_Yeast_Project/07_other/output/BLM_INDEL_Size_Sample_Summary.csv", 
          row.names = FALSE, 
          quote = FALSE)

background_summary <- indel_summary_table %>%
  mutate(background = sub("^[^_]+_([^_]+)_.*$", "\\1", Sample)) %>%
  dplyr::select(-Sample) %>%
  group_by(background) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA")
background_summary <- background_summary %>%
  mutate(background = factor(background, levels = strain_order)) %>%
  arrange(background)

write.csv(background_summary, 
          file = "E:/Bleomycin_Research_Yeast_Project/07_other/output/BLM_INDEL_Size_Background_Summary.csv", 
          row.names = FALSE, 
          quote = FALSE)

#deletion_1bp(AGCT)
deletion_at_details <- combined_data %>%
  filter(Variant_type == "INDEL") %>%
  mutate(
    ref_len = nchar(REF),
    alt_len = nchar(ALT),
    len_diff = alt_len - ref_len
  ) %>%
  filter(len_diff == -1) %>%
  mutate(
    deleted_base = str_remove(REF, fixed(ALT))
  ) %>%
  
  filter(deleted_base %in% c("A", "T", "C", "G")) %>%
  mutate(deleted_base = toupper(deleted_base))

sample_deletion_at <- deletion_at_details %>%
  group_by(Sample, deleted_base) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(column_title = paste0("deletion_", deleted_base)) %>%
  dplyr::select(-deleted_base) %>%
  pivot_wider(names_from = column_title, values_from = count, values_fill = 0)

sample_deletion_at <- sample_deletion_at %>% 
  dplyr::select(Sample, deletion_A, deletion_T, deletion_C, deletion_G)

write.csv(sample_deletion_at, 
          file = "E:/Bleomycin_Research_Yeast_Project/07_other/output/BLM_Deletion_1bp_ATCG_Count_Sample.csv", 
          row.names = FALSE, quote = FALSE)


background_deletion_at <- sample_deletion_at %>%
  mutate(background = sub("^[^_]+_([^_]+)_.*$", "\\1", Sample)) %>%
  dplyr::select(-Sample) %>%
  group_by(background) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  mutate(background = factor(background, levels = c("BY", "WE", "WA", "NA", "SA", "MA"))) %>%
  arrange(background)

write.csv(background_deletion_at, 
          file = "E:/Bleomycin_Research_Yeast_Project/07_other/output/BLM_Deletion_1bp_ATCG_Count_Background.csv", 
          row.names = FALSE, quote = FALSE)





