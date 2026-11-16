library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)

combined_data <- read.csv("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/snv_stats.csv", na.strings = "")

##fig3F_top
ts_tv_summary <- combined_data %>%
  filter(variant_class %in% c("Transition", "Transversion")) %>%
  group_by(Sample, dummy, strain, drug, replicate) %>%
  summarize(
    Ts = sum(count[variant_class == "Transition"], na.rm = TRUE),
    Tv = sum(count[variant_class == "Transversion"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Ts/Tv` = ifelse(Tv > 0, round(Ts / Tv, 2), NA)
  )
write.csv(ts_tv_summary, "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/Sample_Ts_Tv_Summary.csv", row.names = FALSE)

blm_ts_tv <- ts_tv_summary %>%
  filter(drug == "BLM")
strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA")
blm_ts_tv <- blm_ts_tv %>%
  mutate(strain = factor(strain, levels = strain_order))
academic_colors <- c(
  "BY" = "#1f77b4", 
  "WE" = "#ff7f0e", 
  "WA" = "#2ca02c", 
  "NA" = "#d62728",  
  "SA" = "#9467bd",
  "MA" = "#17becf"  
)

ts_tv <- ggplot(blm_ts_tv, aes(x = strain, y = `Ts/Tv`, color = strain)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             aes(fill = strain), size = 1.8, shape = 21, stroke = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  
  labs(
    title = "Ts/Tv Ratio under bleomycin",
    x = NULL,
    y = "Ts/Tv"
  ) +
  scale_fill_manual(values = academic_colors)+
  scale_color_manual(values = academic_colors)+
  theme_bw() +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20, "pt")
  )

print(ts_tv)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3F_top.pdf", ts_tv, width = 4, height = 2.5, dpi = 300, device = cairo_pdf)

#statistical test
analysis_data <- blm_ts_tv %>%
  mutate(strain = factor(strain, levels = c("BY", "WE", "WA", "NA", "SA", "MA"))) %>%
  filter(!is.na(`Ts/Tv`))
global_kw <- kruskal.test(`Ts/Tv` ~ strain, data = analysis_data)
print(global_kw)
pairwise.wilcox.test(
  x = analysis_data$`Ts/Tv`, 
  g = analysis_data$strain, 
  p.adjust.method = "BH",   
  exact = FALSE            
)
#no sig


######fig3G_top
gc_at_summary <- combined_data %>%
  group_by(Sample, dummy, strain, drug, replicate) %>%
  summarize(
    `GC->AT` = sum(count[substitution %in% c("C>T", "C>A")], na.rm = TRUE),
    `AT->GC` = sum(count[substitution %in% c("T>C", "T>G")], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `(GC->AT)/(AT->GC)` = ifelse(`AT->GC` > 0, round(`GC->AT` / `AT->GC`, 2), NA)
  )
write.csv(gc_at_summary, "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/Sample_GC_AT_summary.csv", row.names = FALSE)

blm_gc_at <- gc_at_summary %>%
  filter(drug == "BLM") 
strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA")
blm_gc_at <- blm_gc_at %>%
  mutate(strain = factor(strain, levels = strain_order))
academic_colors <- c(
  "BY" = "#1f77b4",  
  "WE" = "#ff7f0e",  
  "WA" = "#2ca02c",  
  "NA" = "#d62728",  
  "SA" = "#9467bd",  
  "MA" = "#17becf"  
)
gc_at <- ggplot(blm_gc_at, aes(x = strain, y = `(GC->AT)/(AT->GC)`, color = strain)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             aes(fill = strain), size = 1.8, shape = 21, stroke = 0.8) +
  scale_y_continuous(breaks = 0:6, 
                     limits = c(0, 7), 
                     expand = expansion(mult = c(0, 0))) +
  
  labs(
    title = "(GC->AT) / (AT->GC) Ratio under BLM",
    x = NULL,
    y = "(GC->AT) / (AT->GC)"
  ) +
  scale_fill_manual(values = academic_colors)+
  scale_color_manual(values = academic_colors)+
  theme_bw() +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20, "pt")
  )

print(gc_at)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3G_top.pdf", gc_at, width = 4, height = 2.5, dpi = 300, device = cairo_pdf)

#significant
analysis_data <- blm_gc_at %>%
  mutate(strain = factor(strain, levels = c("BY", "WE", "WA", "NA", "SA", "MA"))) %>%
  filter(!is.na(`(GC->AT)/(AT->GC)`))
global_kw <- kruskal.test(`(GC->AT)/(AT->GC)` ~ strain, data = analysis_data)
print(global_kw)
pairwise.wilcox.test(
  x = analysis_data$`(GC->AT)/(AT->GC)`, 
  g = analysis_data$strain, 
  p.adjust.method = "BH",   
  exact = FALSE            
)
#no sig
