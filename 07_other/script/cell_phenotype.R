library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(emmeans)
library(stringr)
library(patchwork)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/07_other/input/cell_phenotype_stats.xlsx", sheet = 1)
data_clean <- data %>%
  filter(is.na(Uncertain) | Uncertain == 0) 
data_clean <- data_clean %>%
  mutate(Line_Name = sub("_[0-9]+$", "", Name))
stats <- data_clean %>%
  group_by(Line_Name, Drug, Background, Rep) %>%
  summarise(
    Total_Valid = n(),
    Normal_count = sum(Normal == 1, na.rm = TRUE),
    Incomplete_cell_separation_count = sum(`Incomplete cell separation` == 1, na.rm = TRUE),
    Flocculation_count = sum(Flocculation == 1, na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  mutate(
    Normal_ratio = Normal_count / Total_Valid,
    Incomplete_cell_separation_ratio = Incomplete_cell_separation_count / Total_Valid,
    Flocculation_ratio = Flocculation_count / Total_Valid,
    Abnormal_ratio = (Incomplete_cell_separation_count + Flocculation_count) / Total_Valid
  )
write.csv(stats, "E:/Bleomycin_Research_Yeast_Project/07_other/output/cell_phenotype_stats_summary.csv", row.names = FALSE)



heat_data <- stats %>%
  filter(Drug %in% c("D0", "YPD", "BLM"))
heat_data <- heat_data %>%
  mutate(X_Label = ifelse(Drug == "D0", "D0", paste0(Drug, "_", Rep)))

x_levels <- c("D0", paste0("YPD_", 1:8), paste0("BLM_", 1:8))
heat_data$X_Label <- factor(heat_data$X_Label, levels = x_levels)
heat_data$Background <- factor(heat_data$Background, levels = c("BY", "WE", "WA", "NA", "SA", "MA"))
heat_data_long <- heat_data %>%
  dplyr::select(X_Label, Background, Incomplete_cell_separation_ratio, Flocculation_ratio) %>%
  pivot_longer(
    cols = c(Incomplete_cell_separation_ratio, Flocculation_ratio),
    names_to = "Phenotype",
    values_to = "Ratio"
  ) %>%
  mutate(Phenotype = dplyr::recode(Phenotype, 
                                   "Incomplete_cell_separation_ratio" = "Incomplete cell separation",
                                   "Flocculation_ratio" = "Flocculation"))
heat_data_long$Background <- factor(heat_data_long$Background, 
                                    levels = c("MA", "SA", "NA", "WA", "WE", "BY"))
heat_data_long$Phenotype <- factor(heat_data_long$Phenotype, 
                                   levels = c("Incomplete cell separation", "Flocculation"))
p <- ggplot(heat_data_long, aes(x = X_Label, y = Background, fill = Ratio)) +
  
  geom_tile(color = "gray40", linewidth = 0.8) +
  facet_wrap(~ Phenotype, ncol = 1, strip.position = "right") +
  
  scale_fill_gradient2(low = "white", mid = "#EEA78B", high = "#851718", 
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Ratio") +
  
  theme_minimal() +
  labs(
    title = "Phenotypic evolution heatmap",
    x = NULL,
    y = "Genetic background"
  ) +
  coord_fixed(ratio = 1)+
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    strip.background = element_rect(fill = "#E5E5E5", color = NA),
    plot.margin = margin(10, 30, 10, 20, "pt"),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm") 
  )
print(p)

ggsave("E:/Bleomycin_Research_Yeast_Project/07_other/output/cell_phenotype_result.pdf", p, width = 9, height = 7, device = cairo_pdf, dpi = 300)

###significant
ancestor_data <- stats %>%
  filter(is.na(Rep) | Rep == "" | trimws(Rep) == "")
evo_data <- stats %>%
  filter(!is.na(Rep) & Rep != "" & trimws(Rep) != "")


backgrounds <- unique(evo_data$Background)

test_Incomplete_cell_separation <- wilcox.test(Incomplete_cell_separation_ratio ~ Drug, data = evo_data)
test_Incomplete_cell_separation
test_Flocculation <- wilcox.test(Flocculation_ratio ~ Drug, data = evo_data)
test_Flocculation
