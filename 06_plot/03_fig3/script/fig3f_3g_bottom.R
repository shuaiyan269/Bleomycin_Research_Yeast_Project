library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/journal.pgen.1011012.s005.xlsx", sheet = 1)
data <- data %>%
  filter(evolved_condition == "YPD") %>%
  dplyr::select(strain, replicate, ref, alt) %>%
  mutate(strain_replicate = paste(strain, replicate, sep = "_"))

####fig3f_bottom
YPD_summary <- data %>%
  mutate(
    class = case_when(
      (ref == "A" & alt == "G") | (ref == "G" & alt == "A") ~ "Ts",
      (ref == "C" & alt == "T") | (ref == "T" & alt == "C") ~ "Ts",
      TRUE ~ "Tv"
    )
  ) %>%
  group_by(strain, strain_replicate) %>%
  summarise(
    drug = "YPD",
    Ts = sum(class == "Ts"),
    Tv = sum(class == "Tv"),
    'Ts/Tv' = Ts / Tv,
    .groups = "drop"
  )

BLM_TS_TV <- read.csv("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/Sample_Ts_Tv_Summary.csv", na.strings = "")
BLM_TS_TV <- BLM_TS_TV %>%
  filter(drug == "BLM") %>%
  rename(`Ts/Tv` = Ts.Tv)

blm_vals <- as.numeric(BLM_TS_TV$`Ts/Tv`)
blm_vals <- blm_vals[!is.na(blm_vals)]

ypd_vals <- YPD_summary$`Ts/Tv`
ypd_vals <- ypd_vals[!is.na(ypd_vals)]

cat("YPD n =", length(ypd_vals), " median =", median(ypd_vals), "\n")
cat("BLM n =", length(blm_vals), " median =", median(blm_vals), "\n")

wilcox.test(ypd_vals, blm_vals)

plot_df <- bind_rows(
  data.frame(group = "YPD", value = ypd_vals),
  data.frame(group = "BLM", value = blm_vals)
)
plot_df$group <- factor(plot_df$group,
                        levels = c("YPD", "BLM"),
                        labels = c("YPD\n(Li et al., 2023)", "BLM\n(This study)"))
p_val <- wilcox.test(ypd_vals, blm_vals)$p.value

ts_tv <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7,
               fill = "white") +                       
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +   
  scale_color_manual(values = c("YPD\n(Li et al., 2023)" = "#4C72B0", "BLM\n(This study)" = "#C44E52")) +
  labs(x = NULL, y = "Ts/Tv ratio") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20, "pt")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.x = 2,
                     label.y = max(plot_df$value) * 1.05,
                     color = "black")

print(ts_tv)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3F_bottom.pdf", ts_tv, width = 4, height = 2.5, dpi = 300, device = cairo_pdf)



####fig3G_bottom

YPD_GCAT <- data %>%
  mutate(
    gc_at = case_when(
      ref %in% c("G", "C") & alt %in% c("A", "T") ~ "GC_AT",   # GC → AT
      ref %in% c("A", "T") & alt %in% c("G", "C") ~ "AT_GC",   # AT → GC
      TRUE ~ NA_character_                                       # GC→GC, AT→AT
    )
  ) %>%
  filter(!is.na(gc_at)) %>%
  group_by(strain, strain_replicate) %>%
  summarise(
    drug              = "YPD",
    `GC->AT`          = sum(gc_at == "GC_AT"),
    `AT->GC`          = sum(gc_at == "AT_GC"),
    `(GC->AT)/(AT->GC)` = ifelse(`GC->AT` == 0 | `AT->GC` == 0,
                                 NA_real_,
                                 `GC->AT` / `AT->GC`),
    .groups = "drop"
  )


#this study
BLM_GCAT <- read.csv("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/Sample_GC_AT_Summary.csv", check.names = FALSE, na.strings = "")

BLM_GCAT <- BLM_GCAT %>%
  filter(drug == "BLM")

BLM_GCAT <- BLM_GCAT %>%
  mutate(`(GC->AT)/(AT->GC)` = as.numeric(`(GC->AT)/(AT->GC)`))
# YPD
ypd_vals <- YPD_GCAT$`(GC->AT)/(AT->GC)`
ypd_vals <- ypd_vals[!is.na(ypd_vals)]

# BLM
blm_vals <- BLM_GCAT$`(GC->AT)/(AT->GC)`
blm_vals <- blm_vals[!is.na(blm_vals)]

cat("YPD n =", length(ypd_vals), " median =", median(ypd_vals), "\n")
cat("BLM n =", length(blm_vals), " median =", median(blm_vals), "\n")

wilcox.test(ypd_vals, blm_vals)

plot_df <- bind_rows(
  data.frame(group = "YPD", value = ypd_vals),
  data.frame(group = "BLM", value = blm_vals)
)
plot_df$group <- factor(plot_df$group,
                        levels = c("YPD", "BLM"),
                        labels = c("YPD\n(Li et al., 2023)", "BLM\n(This study)"))
p_val <- wilcox.test(ypd_vals, blm_vals)$p.value

GC_AT <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7,
               fill = "white") +                       
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +   
  scale_color_manual(values = c("YPD\n(Li et al., 2023)" = "#4C72B0", "BLM\n(This study)" = "#C44E52")) +
  labs(x = NULL, y = "GC->AT/AT->GC ratio") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20, "pt")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.x = 2,
                     label.y = max(plot_df$value) * 1.05,
                     color = "black")

print(GC_AT)
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3G_bottom.pdf", GC_AT, width = 4, height = 2.5, dpi = 300, device = cairo_pdf)







