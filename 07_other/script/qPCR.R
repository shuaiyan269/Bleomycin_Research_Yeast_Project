##
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

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/07_other/input/qPCR.xlsx", sheet = 1)
Data_avg <- data %>%
  group_by(Background, Drug, rep, Target) %>%
  summarise(
    Mean_Cq = mean(Cq, na.rm = TRUE),
    .groups = "drop" 
  )

Data_wide <- Data_avg %>%
  pivot_wider(names_from = Target, values_from = Mean_Cq)

Data_calc <- Data_wide %>%
  mutate(
    dCt_REV1 = REV1 - ACT1,
    dCt_REV3 = REV3 - ACT1,
    RelExp_REV1 = 2^(-dCt_REV1),
    RelExp_REV3 = 2^(-dCt_REV3)
  )

plot_data_REV1 <- Data_calc %>%
  group_by(Background, Drug) %>%
  summarise(
    Mean = mean(RelExp_REV1, na.rm = TRUE),
    SEM = sd(RelExp_REV1, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plot_data_REV3 <- Data_calc %>%
  group_by(Background, Drug) %>%
  summarise(
    Mean = mean(RelExp_REV3, na.rm = TRUE),
    SEM = sd(RelExp_REV3, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

aov_REV1 <- aov(RelExp_REV1 ~ Background * Drug, data = Data_calc)
summary(aov_REV1)
emm <- emmeans(aov_REV1, ~ Background | Drug)
pairwise_res <- pairs(emm)
res_table_REV1 <- as.data.frame(pairwise_res)

print(res_table_REV1)


p1 <- ggplot(Data_calc, aes(x = Background, y = RelExp_REV1, fill = Drug)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), 
               width = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8), width = 0.2) +
  geom_jitter(position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7, color = "black") +
  labs(
    title = "REV1 relative expression",
    x = "Genetic background",
    y = "Relative expression (2^-dCt)",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4, size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )+
  scale_x_discrete(limits = c("BY", "WE", "WA", "NA", "SA", "MA"))+
  scale_fill_manual(values = c("YPD" = "grey80", "BLM" = "#78B6C8")) +
  stat_compare_means(aes(group = Drug), method = "t.test", label = "p.signif", hide.ns = FALSE,
                     label.y = max(Data_calc$RelExp_REV1, na.rm = TRUE) * 1.15) +
  annotate("text", 
           x = 6.5, 
           y = max(Data_calc$RelExp_REV1, na.rm = TRUE) * 1.5, 
           label = "Two-way ANOVA:\nBackground: p < 0.001 ***\nDrug: p < 0.001 ***\nInteraction: p < 0.001 ***", 
           hjust = 1, 
           vjust = 1, 
           size = 4) + 
  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

print(p1)

#####rev3
aov_REV3 <- aov(RelExp_REV3 ~ Background * Drug, data = Data_calc)
summary(aov_REV3)

emm <- emmeans(aov_REV3, ~ Background | Drug)
pairwise_res <- pairs(emm)
res_table_REV3 <- as.data.frame(pairwise_res)
print(res_table_REV3)
p2 <- ggplot(Data_calc, aes(x = Background, y = RelExp_REV3, fill = Drug)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), 
               width = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8), width = 0.2) +
  geom_jitter(position = position_dodge(width = 0.8), size = 1.5, alpha = 0.7, color = "black") +
  labs(
    title = "REV3 relative expression",
    x = "Genetic background",
    y = "Relative expression (2^-dCt)",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4, size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )+
  scale_x_discrete(limits = c("BY", "WE", "WA", "NA", "SA", "MA"))+
  scale_fill_manual(values = c("YPD" = "grey80", "BLM" = "#78B6C8")) +
  
  stat_compare_means(aes(group = Drug), method = "t.test", label = "p.signif", hide.ns = FALSE,
                     label.y = max(Data_calc$RelExp_REV3, na.rm = TRUE) * 1.15)+
  annotate("text", 
           x = 6.5, 
           y = max(Data_calc$RelExp_REV3, na.rm = TRUE) * 1.5, 
           label = "Two-way ANOVA:\nBackground: p < 0.001 ***\nDrug: p = 0.003 **\nInteraction: p = 0.123 ns", 
           hjust = 1, 
           vjust = 1, 
           size = 4)+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))
print(p2)


#####BLM compare to YPD
Data_FC <- Data_calc %>%
  dplyr::select(Background, rep, Drug, RelExp_REV1, RelExp_REV3) %>%
  pivot_wider(names_from = Drug, values_from = c(RelExp_REV1, RelExp_REV3)) %>%
  mutate(
    FC_REV1 = RelExp_REV1_BLM / RelExp_REV1_YPD,
    FC_REV3 = RelExp_REV3_BLM / RelExp_REV3_YPD
  ) %>%
  drop_na(FC_REV1, FC_REV3) 
Data_FC_long <- Data_FC %>%
  pivot_longer(
    cols = c(FC_REV1, FC_REV3), 
    names_to = "Gene", 
    values_to = "FoldChange"
  ) %>%
  mutate(Gene = gsub("FC_", "", Gene))

mutation_rate <- read_excel("E:/Bleomycin_Research_Yeast_Project/07_other/input/qPCR.xlsx", sheet = 2)

Data_FC_colored <- Data_FC_long %>%
  left_join(mutation_rate, by = "Background")

p3 <- ggplot(Data_FC_colored, aes(x = Background, y = FoldChange, fill = `SNV_rate`)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "BLM compared to YPD",
    x = "Genetic background",
    y = "Relative expression level",
    fill = "SNV Rate" 
  ) +
  scale_x_discrete(limits = c("BY", "WE", "WA", "NA", "SA", "MA")) +
  scale_fill_gradientn(
    colors = c("#ffffbf", "#d73027"), 
    name = "SNV Rate",
    breaks = c(2e-9, 6e-9, 10e-9, 14e-9),
    labels = c("2x10-9", "6x10-9", "10x10-9", "14x10-9")
  ) +
  
  facet_wrap(~ Gene) + 
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
p <- (p1+p2)/p3
print(p)
ggsave("E:/Bleomycin_Research_Yeast_Project/07_other/output/qPCR.pdf", p, width = 12, height = 12, dpi = 300, device = cairo_pdf)

##REV1
Data_corr_REV1 <- Data_FC_long %>%
  filter(Gene == "REV1") %>%
  group_by(Background) %>%
  summarise(Mean_FC = mean(FoldChange, na.rm = TRUE)) %>%
  left_join(mutation_rate, by = "Background")


print(Data_corr_REV1)
cor_test_REV1 <- cor.test(Data_corr_REV1$Mean_FC, Data_corr_REV1$SNV_rate, method = "spearman")
print(cor_test_REV1)


##REV3
Data_corr_REV3 <- Data_FC_long %>%
  filter(Gene == "REV3") %>%
  group_by(Background) %>%
  summarise(Mean_FC = mean(FoldChange, na.rm = TRUE)) %>%
  left_join(mutation_rate, by = "Background")


print(Data_corr_REV3)
cor_test_REV3 <- cor.test(Data_corr_REV3$Mean_FC, Data_corr_REV3$SNV_rate, method = "spearman")
print(cor_test_REV3)


