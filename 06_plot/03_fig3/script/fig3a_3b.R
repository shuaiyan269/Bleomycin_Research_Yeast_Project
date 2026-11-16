library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)
library(car)
library(multcomp)

#####mutation_rate####
data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/mutation_rate.xlsx", sheet = 1)

####fig3A_SNV_rate
data_cleaned <- data %>%
  separate(
    col = Samples, 
    into = c("Day", "Background", "Treatment", "Line"), 
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("YPD", "BLM")),
    Background = factor(Background, levels = c("BY", "WE", "WA", "NA", "SA", "MA")),
    SNV_scaled = SNV * 1e10  
  )
data_long <- data_cleaned %>%
  pivot_longer(
    cols = c(InDel, SNV), 
    names_to = "Mut_Type", 
    values_to = "Rate"
  )
data_blm <- data_cleaned %>% filter(Treatment == "BLM")
data_blm$Background <- relevel(factor(data_blm$Background), ref = "BY")
leveneTest(SNV_scaled ~ Background, data = data_blm)
model_blm <- lm(SNV_scaled ~ Background, data = data_blm)
dunnett_res <- glht(model_blm, linfct = mcp(Background = "Dunnett"))
dunnett_summary <- summary(dunnett_res)
p_values <- dunnett_summary$test$pvalues
comp_names <- names(p_values)
if(is.null(comp_names)) {
  
  comp_names <- levels(data_blm$Background)[-1] 
} else {
  comp_names <- str_remove(comp_names, " - BY")
}
p_data <- tibble(
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
print(p_data)
max_y <- max(data_cleaned$SNV_scaled, na.rm = TRUE)
p_data <- p_data %>%
  mutate(
    y.position = max_y * 1.05 + (row_number() * max_y * 0.07)
  )

snv_rate <- ggplot(data_cleaned, aes(x = Background, y = SNV_scaled, color = Treatment)) + 
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6, 
    size = 1.8, 
    shape = 16 
  ) +
  scale_color_manual(values = c("YPD" = "gray50", "BLM" = "#7DAEE0")) +
  stat_pvalue_manual(
    p_data, 
    label = "label", 
    xmin = "group1", 
    xmax = "group2", 
    y.position = "y.position",
    color = "black", 
    tip.length = 0.02
  ) +
  
  labs(
    title = "SNV rate",
    x = "Genetic Background",
    y = expression(10^{-10}/base/generation), 
    color = "Treatment"
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
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
print(snv_rate)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3A.pdf", snv_rate, width = 6, height = 4, dpi = 300, device = cairo_pdf)

###fig3B_INDEL_rate
data_cleaned <- data %>%
  separate(
    col = Samples, 
    into = c("Day", "Background", "Treatment", "Line"), 
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("YPD", "BLM")),
    Background = factor(Background, levels = c("BY", "WE", "WA", "NA", "SA", "MA")),
    InDel_scaled = InDel * 1e10 
  )

data_blm <- data_cleaned %>% filter(Treatment == "BLM")
data_blm$Background <- relevel(factor(data_blm$Background), ref = "BY")

print(leveneTest(InDel_scaled ~ Background, data = data_blm))

model_blm <- lm(InDel_scaled ~ Background, data = data_blm)
dunnett_res <- glht(model_blm, linfct = mcp(Background = "Dunnett"))
dunnett_summary <- summary(dunnett_res)

p_values <- dunnett_summary$test$pvalues
comp_names <- names(p_values)

if(is.null(comp_names)) {
  comp_names <- levels(data_blm$Background)[-1] 
} else {
  comp_names <- str_remove(comp_names, " - BY")
}

p_data <- tibble(
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

print(p_data)
max_y <- max(data_cleaned$InDel_scaled, na.rm = TRUE)
p_data <- p_data %>%
  mutate(
    y.position = max_y * 1.05 + (row_number() * max_y * 0.07)
  )


indel_rate <- ggplot(data_cleaned, aes(x = Background, y = InDel_scaled, color = Treatment)) + 
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6, 
    size = 1.8, 
    shape = 16 
  ) +
  scale_color_manual(values = c("YPD" = "gray50", "BLM" = "#7DAEE0")) +
  
  
  stat_pvalue_manual(
    p_data, 
    label = "label", 
    xmin = "group1", 
    xmax = "group2", 
    y.position = "y.position",
    color = "black", 
    tip.length = 0.02
  ) +
  
  labs(
    title = "INDEL rate",
    x = "Genetic Background",
    y = expression(10^{-10}/base/generation), 
    color = "Treatment"
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
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
print(indel_rate)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3B.pdf", indel_rate, width = 6, height = 4, dpi = 300, device = cairo_pdf)





