
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(ggpubr)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/input/data_fig5C.xlsx", sheet = 1)
my_comparisons <- list(
  c("Δpma1", "wt"),
  c("Δptk2", "wt")
)
fig5C <- ggplot(data, aes(x = mutation, y = median, color = mutation)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 1) +
  labs(x = "Strains", y = "median FI", subtitle = "Wilcoxon test", color = NULL) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 3,
                     bracket.size = 0.3)+
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = c(
    "wt" = "gray60",
    "Δptk2" = "#71b8ed",
    "Δpma1" = "#ea8379"
  ))+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.subtitle = element_text(hjust = 0, size = 10, face = "italic"),
        plot.margin = margin(20, 20, 20, 20, "pt"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid.major = element_line(color = "gray80", linewidth = 0.3),  
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.2))

print(fig5C)


ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/fig5C.pdf", fig5C, width = 6, height = 6, dpi = 300, device = cairo_pdf)
