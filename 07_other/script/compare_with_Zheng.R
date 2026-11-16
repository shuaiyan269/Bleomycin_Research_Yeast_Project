#Compare mutation rates (SNV / INDEL) between this study and Zheng et al.

library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/07_other/input/compare_with_Zheng.xlsx", sheet = 1)

targets <- c(
  "Zheng et al.", "BY", "WE", "WA", "NA", "SA", "MA"
)

academic_colors <- c(
  "BY" = "#1f77b4",
  "WE" = "#ff7f0e",
  "WA" = "#2ca02c",
  "NA" = "#d62728",
  "SA" = "#9467bd",
  "MA" = "#17becf",
  "Zheng et al." = "gray60"
)

data_long <- data %>%
  mutate(InDel_rate = InDel_rate * 1e10, SNV_rate = SNV_rate * 1e10) %>%
  pivot_longer(cols = c(SNV_rate, InDel_rate), names_to = "Mutation_Type", values_to = "Rate") %>%
  mutate(Mutation_Type = dplyr::recode(Mutation_Type, "SNV_rate" = "SNV rate", "InDel_rate" = "INDEL rate"
  ),
  Mutation_Type = factor(Mutation_Type, levels = c("SNV rate", "INDEL rate")),
  Source = dplyr::recode(Source, "Zheng DQ et al.," = "Zheng et al."),
  Group = if_else(Source == "Zheng et al.", "Zheng et al.", as.character(Background)),
  Group = factor(Group, levels = targets)
  )

kw_result <- data_long %>%
  filter(Group %in% targets) %>%
  group_by(Mutation_Type) %>%
  kruskal_test(Rate ~ Group)

print(kw_result)

#Wilcoxon rank-sum test (ref = Zheng et al.), BH-adjusted p values.
res <- data_long %>%
  filter(Group %in% targets) %>%
  group_by(Mutation_Type) %>%
  wilcox_test(
    Rate ~ Group,
    ref.group = "Zheng et al."
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(res)

res_plot <- res %>%
  mutate(xmin = 1,
         xmax = match(as.character(group2), targets)
  ) %>%
  group_by(Mutation_Type) %>%
  arrange(xmax, .by_group = TRUE) %>%
  mutate(
    y.position = max(
      data_long$Rate[
        data_long$Mutation_Type == first(Mutation_Type) &
          data_long$Group %in% targets
      ],
      na.rm = TRUE
    ) + row_number() * 0.08 *
      diff(range(
        data_long$Rate[
          data_long$Mutation_Type == first(Mutation_Type) &
            data_long$Group %in% targets
        ],
        na.rm = TRUE
      ))
  ) %>%
  ungroup()


compare <- ggplot(data_long %>% filter(Group %in% targets), aes(x = Group, y = Rate, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  facet_wrap(~ Mutation_Type, scales = "free_y") +
  stat_pvalue_manual(
    res_plot,
    label = "p.adj.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3
  ) +
  scale_x_discrete(
    limits = targets,
    drop = FALSE
  ) +
  scale_color_manual(
    values = academic_colors,
    breaks = targets,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = expression(10^{-10}/base/generation)
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 10, color = "black"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 10, color = "black"),
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.8
    ),
    axis.line = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
print(compare)
ggsave("E:/Bleomycin_Research_Yeast_Project/07_other/output/compare_with_Zheng.pdf", compare, width = 8, height = 4.5, dpi = 300, device = cairo_pdf)



