library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(ggpubr)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/input/data_fig5DE.xlsx", sheet = 1)
df_long <- data %>%
  pivot_longer(
    cols = -time,
    names_to = "sample",
    values_to = "fluorescence"
  ) %>%
  mutate(
    group = case_when(
      str_detect(sample, "Δptk2") ~ "Δptk2",
      str_detect(sample, "Δpma1") ~ "Δpma1", 
      str_detect(sample, "wt") ~ "WT",
      str_detect(sample, "control") ~ "Control",
      TRUE ~ "Other"
    )
  )
df_summary <- df_long %>%
  group_by(time, group) %>%
  summarise(
    mean_fluo = mean(fluorescence, na.rm = TRUE),
    se_fluo = sd(fluorescence, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

##fig5D
fig5D <- ggplot(df_summary, aes(x = time, y = mean_fluo, color = group, group = group)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_fluo - se_fluo, ymax = mean_fluo + se_fluo), 
                width = 0.2, alpha = 0.5) +
  labs(
    x = "Time",
    y = "FI"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    axis.title = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c("Δptk2" = "#69B7CE", 
               "Δpma1" = "#F4CE4B", 
               "WT" = "#C65840", 
               "Control" = "gray70")
  )+
  scale_y_continuous(
    breaks = c(3000, 6000, 9000, 12000, 15000), 
    limits = c(3000, 16000),                       
    expand = c(0, 0)                           
  )
print(fig5D)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/fig5D.pdf", fig5D,  width = 5, height = 3, dpi = 300, device = cairo_pdf)


# Wilcoxon
df_long_filtered <- df_long %>% filter(group != "Control")
stat_test <- df_long_filtered %>%
  group_by(time) %>%
  summarise(
    median_ptk2 = median(fluorescence[group == "Δptk2"], na.rm = TRUE),
    median_pma1 = median(fluorescence[group == "Δpma1"], na.rm = TRUE),
    median_WT = median(fluorescence[group == "WT"], na.rm = TRUE),
    p_ptk2 = wilcox.test(fluorescence[group == "Δptk2"], 
                         fluorescence[group == "WT"],
                         alternative = "greater")$p.value,
    p_pma1 = wilcox.test(fluorescence[group == "Δpma1"], 
                         fluorescence[group == "WT"],
                         alternative = "greater")$p.value,
    
    max_y = max(fluorescence)
  ) %>%
  mutate(
    sig_ptk2 = case_when(
      p_ptk2 < 0.001 ~ "***",
      p_ptk2 < 0.01 ~ "**",
      p_ptk2 < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_pma1 = case_when(
      p_pma1 < 0.001 ~ "***",
      p_pma1 < 0.01 ~ "**",
      p_pma1 < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(stat_test, file = "E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/membrane_potential_statistics.csv", row.names = FALSE)

###boxplot###
###fig5E
fig5E <-ggplot(df_long %>% filter(group != "Control"), 
       aes(x = factor(time), y = fluorescence, fill = group, color = group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.7, outlier.color = "black") +
  scale_fill_manual(
    values = c("Δptk2" = "#69B7CE", 
               "Δpma1" = "#F4CE4B", 
               "WT" = "#C65840")) +
  scale_color_manual(
    values = c("Δptk2" = "#69B7CE", 
               "Δpma1" = "#F4CE4B", 
               "WT" = "#C65840")
  )+
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)) - 0.25, 
                   xend = as.numeric(factor(time)) + 0.25,
                   y = max_y * 1.08, 
                   yend = max_y * 1.08),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)) - 0.25, 
                   xend = as.numeric(factor(time)) - 0.25,
                   y = max_y * 1.08, 
                   yend = max_y * 1.06),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)) + 0.25, 
                   xend = as.numeric(factor(time)) + 0.25,
                   y = max_y * 1.08, 
                   yend = max_y * 1.06),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_text(data = stat_test,
            aes(x = factor(time), y = max_y * 1.1, 
                label = ifelse(sig_ptk2 != "ns", sig_ptk2, "ns")),
            inherit.aes = FALSE,
            color = "#69B7CE", size = 3, fontface = "bold") +
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)), 
                   xend = as.numeric(factor(time)) - 0.25,
                   y = max_y * 1.15, 
                   yend = max_y * 1.15),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)), 
                   xend = as.numeric(factor(time)),
                   y = max_y * 1.15, 
                   yend = max_y * 1.13),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_segment(data = stat_test,
               aes(x = as.numeric(factor(time)) - 0.25, 
                   xend = as.numeric(factor(time)) - 0.25,
                   y = max_y * 1.15, 
                   yend = max_y * 1.13),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.25) +
  geom_text(data = stat_test,
            aes(x = as.numeric(factor(time)) - 0.125, y = max_y * 1.17, 
                label = ifelse(sig_pma1 != "ns", sig_pma1, "ns")),
            inherit.aes = FALSE,
            color = "#F4CE4B", size = 3, fontface = "bold") +
  
  labs(
    x = "Time (min)",
    y = "FI"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )+
  scale_y_continuous(
    breaks = c(3000, 6000, 9000, 12000, 15000, 18000), 
    limits = c(2000, 20000),                       
    expand = c(0, 0)                           
  )

print(fig5E)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/fig5E.pdf", fig5E, width = 12, height = 6, dpi = 300, device = cairo_pdf)





