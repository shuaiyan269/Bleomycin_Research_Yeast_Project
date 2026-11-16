library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
all_sample_names <- readLines("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/sample_name.txt")
####six substitution type stats####
data <- data %>%
  filter(Sample_type == "population")
substitution_types <- c(
  "C>A", "C>G", "C>T",    
  "T>A", "T>C", "T>G"    
)
all_combinations <- expand_grid(
  Sample = all_sample_names,
  substitution = substitution_types  
)

raw_counts <- data %>%
  filter(Variant_type == "SNV") %>%
  mutate(
    substitution = case_when(
      REF %in% c("C", "T") ~ paste0(REF, ">", ALT),
      REF == "A" ~ paste0("T>", ifelse(ALT == "A", "A", 
                                       ifelse(ALT == "C", "G", 
                                              ifelse(ALT == "G", "C", "A")))),
      REF == "G" ~ paste0("C>", ifelse(ALT == "A", "T", 
                                       ifelse(ALT == "C", "G", 
                                              ifelse(ALT == "G", "C", "A"))))
    )
  ) %>%
  count(Sample, substitution, name = "count")
raw_counts <- raw_counts %>%
  mutate(
    substitution = ifelse(Sample == "D90_NA_YPD_3" & substitution == "C>T,A", 
                          "C>T", 
                          substitution)
  )

snv_stats <- all_combinations %>%
  left_join(raw_counts, by = c("Sample", "substitution")) %>%
  mutate(count = coalesce(count, 0))

snv_stats <- snv_stats %>%
  group_by(Sample) %>%
  mutate(total_snvs = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = ifelse(total_snvs > 0, count / total_snvs, 0))


snv_stats <- snv_stats %>%
  dplyr::select(Sample, substitution, count, total_snvs, proportion) %>%
  arrange(Sample, substitution) %>%
  separate(
    col = Sample,
    into = c("dummy", "strain", "drug", "replicate"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    variant_class = case_when(
      substitution %in% c("C>T", "T>C") ~ "Transition",  
      TRUE ~ "Transversion"                              
    )
  )

write.csv(snv_stats, "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/snv_stats.csv", row.names = FALSE, fileEncoding = "UTF-8")

#####INDEL_stats####
indel_stats <- data %>%
  filter(Variant_type == "INDEL")
indel <- indel_stats %>%
  mutate(
    Indel_Type = case_when(
      nchar(REF) > nchar(ALT) ~ "Deletion",      
      nchar(REF) < nchar(ALT) ~ "Insertion",     
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV", 
      TRUE ~ "Complex"                            
    ),
    Length_Change = abs(nchar(REF) - nchar(ALT))
  )
all_combinations <- expand_grid(
  Sample = unique(snv_stats$Sample),
  Indel_Type = unique(indel$Indel_Type)
)
variant_counts <- indel %>%
  count(Sample, Indel_Type, name = "Count")
sample_stats <- all_combinations %>%
  left_join(variant_counts, by = c("Sample", "Indel_Type")) %>%
  mutate(Count = coalesce(Count, 0)) %>%
  pivot_wider(
    names_from = Indel_Type, 
    values_from = Count, 
    values_fill = 0
  ) %>%
  mutate(Total = rowSums(across(where(is.numeric))))
sample_stats <- sample_stats %>%
  separate(
    col = Sample,          
    into = c("dummy", "strain", "drug", "replicate"), 
    sep = "_",            
    remove = FALSE      
  )
write.csv(sample_stats, "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/indel_stats.csv", row.names = FALSE, fileEncoding = "UTF-8")

####fig3C
combined_data <- read.csv("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/snv_stats.csv", na.strings = "")
strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA") 
target_order_top_to_bottom <- c("C>T", "T>C", "C>A", "C>G", "T>A", "T>G")
substitution_levels <- rev(target_order_top_to_bottom) 
plot_data <- combined_data %>%
  group_by(strain, drug, substitution) %>%
  summarize(count = sum(count, na.rm = TRUE), .groups = "drop") %>% 
  mutate(
    strain = factor(strain, levels = strain_order),
    
    substitution = factor(substitution, levels = substitution_levels)
  ) %>%
  group_by(strain, drug) %>%
  mutate(
    total_count = sum(count),
    percentage = count / total_count,
    label_text = ifelse(percentage > 0.03, sprintf("%.2f", percentage), "")
  ) %>%
  ungroup()
substitution_colors <- c(
  "C>A" = "#1f77b4",    
  "C>G" = "#ff7f0e",    
  "C>T" = "#2ca02c",    
  "T>A" = "#d62728",   
  "T>C" = "#9467bd",    
  "T>G" = "#8c564b"   
)
total_counts <- plot_data %>%
  group_by(strain, drug) %>%
  summarize(total = first(total_count), .groups = "drop")
snv_count <- ggplot(plot_data, aes(x = drug, y = count, fill = substitution)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.75) + 
  scale_fill_manual(
    values = substitution_colors,
    breaks = target_order_top_to_bottom 
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     expand = expansion(mult = c(0.01, 0.15)))+ 
  geom_text(
    aes(label = label_text),
    position = position_fill(vjust = 0.5, reverse = TRUE),
    color = "black",
    size = 3
  ) +
  geom_text(
    data = total_counts,
    aes(x = drug, y = 1.0, label = total),
    inherit.aes = FALSE, color = "black", size = 3.5, vjust = -0.5       
  ) +
  facet_grid(. ~ strain, scales = "free_x", space = "free_x", switch = "x") +
  labs(title = "SNV substitution proportions: YPD vs BLM", x = NULL, y = "Proportion of SNVs", fill = NULL) +
  theme_bw() + 
  theme(
    text = element_text(size = 10, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 9, color = "black"), 
    strip.text = element_text(size = 10, color = "black"),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank(),      
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
print(snv_count)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3C.pdf", snv_count, width = 7, height = 4, dpi = 300, device = cairo_pdf)


#SNV_statistical_test
contingency_table <- plot_data %>%
  filter(drug == "BLM") %>%  
  dplyr::select(strain, substitution, count) %>%
  pivot_wider(names_from = substitution, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "strain")
mat <- as.matrix(contingency_table)
chi_test <- chisq.test(mat)

print(chi_test)


####fig3D
combined_data <- read.csv("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/indel_stats.csv", na.strings = "")
combined_data <- combined_data %>%
  pivot_longer(
    cols = c(Insertion, Deletion),
    names_to = "Indel_Type",
    values_to = "Count"
  )

strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA") 

indel_order <- c("Insertion", "Deletion")

plot_data <- combined_data %>%
  group_by(strain, drug, Indel_Type) %>%
  summarize(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>% 
  mutate(
    strain = factor(strain, levels = strain_order),
    Indel_Type = factor(Indel_Type, levels = indel_order) 
  ) %>%
  group_by(strain, drug) %>%
  mutate(
    total_count = sum(Count),
    percentage = Count / total_count,
    
    label_text = ifelse(percentage > 0.03, sprintf("%.2f", percentage), "")
  ) %>%
  ungroup()

total_counts <- plot_data %>%
  group_by(strain, drug) %>%
  summarize(total = first(total_count), .groups = "drop")

indel_count <- ggplot(plot_data, aes(x = drug, y = Count, fill = Indel_Type)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.75) + 
  
  scale_fill_manual(values = c("Insertion" = "#0571b0", "Deletion" = "#92c5de"),
                    breaks = indel_order) +
  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     expand = expansion(mult = c(0.01, 0.15))) + 
  geom_text(
    aes(label = label_text),
    position = position_fill(vjust = 0.5, reverse = TRUE),
    color = "black", 
    size = 3
  ) +
  
  geom_text(
    data = total_counts,
    aes(x = drug, y = 1.0, label = total),
    inherit.aes = FALSE, color = "black", size = 3.5, vjust = -0.5       
  ) +
  facet_grid(. ~ strain, scales = "free_x", space = "free_x", switch = "x") +
  
  labs(
    title = "Insertion and Deletion in BLM and YPD",
    x = NULL,         
    y = "Proportion of INDELs",
    fill = NULL
  ) +
  theme_bw() + 
  theme(
    text = element_text(size = 10, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 9, color = "black"), 
    strip.text = element_text(size = 10, color = "black"),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank(),      
    legend.position = "right",          
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
print(indel_count)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3D.pdf", indel_count, width = 7, height = 4, dpi = 300, device = cairo_pdf)









