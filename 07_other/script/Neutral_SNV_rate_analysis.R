## Neutral SNV rate analysis
#This script estimates the neutral SNV mutation rate for six genetic backgrounds (BY, WE, WA, NA, SA, MA) and compares each background to BY using a linear model followed by Dunnett's multiple-comparison test.

library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)
library(multcomp)
library(car)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
data <- data %>%
  filter(Drug == "BLM", Sample_type == "population")

target <- c("synonymous_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "intergenic_variant")

# Exclude genes with recurrent/high-frequency mutations, which may reflect selection rather than neutral evolution
exclude_genes <- c("PMA1", "PTK2", "GIN4", "SFL1", "ACE2", "SKY1", "WHI2", "ASK10", "TRK1", "HSL1", "IRA1") 

neutral <- data %>%
  filter(Consequence %in% target) %>%
  filter(!Gene_name %in% exclude_genes)

snv_count <- neutral %>%
  filter(Variant_type == "SNV") %>%
  count(Sample, Background, name = "n_snv") 

#Genome sizes (bp) for the six backgrounds
#Genome size data: https://www.ncbi.nlm.nih.gov/bioproject/PRJEB7245/
genome_sizes <- data.frame(
  Background = c("BY", "WE", "WA", "NA", "SA", "MA"),   
  genome_size = c(12157149, 11813188, 11958041, 11910474, 11881147, 11747346)        
)
neutral_rate <- snv_count %>%
  left_join(genome_sizes, by = "Background") %>%   
  mutate(rate = n_snv / genome_size / 225 / 2,
         SNV_scaled = rate * 1e10)

neutral_rate <- neutral_rate %>%
  mutate(
    Background = factor(Background, levels = c("BY", "WE", "WA", "NA", "SA", "MA"))
  ) %>%
  mutate(
    Background = relevel(Background, ref = "BY")
  )
levels(neutral_rate$Background)

leveneTest(SNV_scaled ~ Background, data = neutral_rate)
model_neutral <- lm(SNV_scaled ~ Background, data = neutral_rate)
summary(model_neutral)

dunnett_res <- glht(model_neutral, linfct = mcp(Background = "Dunnett"))
dunnett_summary <- summary(dunnett_res)
print(dunnett_summary)

p_values <- dunnett_summary$test$pvalues
comp_names <- names(p_values)

if (is.null(comp_names)) {
  comp_names <- levels(neutral_rate$Background)[-1]
} else {
  comp_names <- str_remove(comp_names, " - BY")
}

p_data_neutral <- tibble(
  group1 = "BY",
  group2 = comp_names,
  p.value = as.numeric(p_values)
) %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )
print(p_data_neutral)

max_y <- max(neutral_rate$SNV_scaled, na.rm = TRUE)
p_data <- p_data_neutral %>%
  mutate(
    y.position = max_y * 1.05 + (row_number() * max_y * 0.07)
  )
academic_colors <- c(
  "BY" = "#1f77b4", 
  "WE" = "#ff7f0e", 
  "WA" = "#2ca02c", 
  "NA" = "#d62728",  
  "SA" = "#9467bd",
  "MA" = "#17becf"  
)

ggplot(neutral_rate, aes(x = Background, y = SNV_scaled, color = Background)) + 
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6, 
    size = 1.8, 
    shape = 16 
  ) +
  stat_pvalue_manual(
    p_data, 
    label = "label", 
    xmin = "group1", 
    xmax = "group2", 
    y.position = "y.position",
    color = "black", 
    tip.length = 0.02
  ) +
  scale_color_manual(values = academic_colors) +
  labs(
    title = "Neutral SNV rate",
    x = "Genetic Background",
    y = expression(10^{-10}/base/generation), 
    color = "background"
  ) +
  theme_bw() + 
  theme(
    text = element_text(size = 10, color ="black"),
    plot.title = element_text(hjust = 0.5, size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.line = element_blank(),
    strip.text = element_text(size = 10), 
    panel.grid.major.x = element_blank(), 
    legend.position = "none",
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )


ggsave("E:/Bleomycin_Research_Yeast_Project/07_other/output/BLM_neutral_SNV_rate.pdf", width = 6, height = 5, dpi = 300, device = cairo_pdf)

