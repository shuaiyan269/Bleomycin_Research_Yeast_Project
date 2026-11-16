##fig1 and fig_s1

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

#input data
data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/input/Table_S1.xlsx", sheet = 1)

relative_inhibition <- data %>%
  mutate(
    Lag_time = as.numeric(Lag_time),
    Generation_time = as.numeric(Generation_time),
    Yield = as.numeric(Yield)
  ) %>%
  dplyr::select(
    Sample, Population, Background, Time, Rep,
    Drug, Lag_time, Generation_time, Yield
  ) %>%
  pivot_wider(
    names_from = Drug,
    values_from = c(Lag_time, Generation_time, Yield),
    names_sep = "_"
  ) %>%
  group_by(Sample, Population, Background, Time) %>%
  summarise(
    Lag_time_BLM = mean(Lag_time_BLM, na.rm = TRUE),
    Lag_time_YPD = mean(Lag_time_YPD, na.rm = TRUE),
    Generation_time_BLM = mean(Generation_time_BLM, na.rm = TRUE),
    Generation_time_YPD = mean(Generation_time_YPD, na.rm = TRUE),
    Yield_BLM = mean(Yield_BLM, na.rm = TRUE),
    Yield_YPD = mean(Yield_YPD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Relative_Lag_time = (Lag_time_BLM - Lag_time_YPD) / Lag_time_YPD,
    Relative_Generation_time = (Generation_time_BLM - Generation_time_YPD) / Generation_time_YPD,
    Relative_Yield = (Yield_BLM - Yield_YPD) / Yield_YPD
  ) %>%
  dplyr::select(
    Sample, Population, Background, Time,
    Relative_Lag_time,
    Relative_Generation_time,
    Relative_Yield
  )

relative_inhibition <- relative_inhibition %>%
  filter(Sample != "D90_MA_BLM_2")

relative_inhibition$Background <- factor(relative_inhibition$Background,
                          levels = c("BY","WE","WA","NA","SA","MA"))

###########fig1c, fig_s1b, and fig_s1c###################################
plot_relative_inhibition <- function(data, value_col, y_label, title) {
  # D90 and D0
  d90_data <- data %>%
    filter(Time == "D90")
  d0_data <- data %>%
    filter(Time == "D0")
  
  # t-test: D90 vs D0 for each Background
  p_values <- data %>%
    group_by(Background) %>%
    summarise(
      p = if(
        length(.data[[value_col]][Time == "D90"]) > 1 &&
        length(.data[[value_col]][Time == "D0"]) == 1
      ) {
        t.test(
          .data[[value_col]][Time == "D90"],
          mu = .data[[value_col]][Time == "D0"]
        )$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  
  # Significance
  p_values$signif <- cut(
    p_values$p,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )
  
  # Y position of significance stars
  star_y <- max(d90_data[[value_col]], na.rm = TRUE) * 1.12
  
  # Plot
  p <- ggplot() +
    
    geom_boxplot(data = d90_data, aes(x = Background, y = .data[[value_col]], fill = "D90"),
      color = "#7DAEE0", width = 0.5, linewidth = 0.6, outlier.shape = NA, alpha = 0.2
    ) +
    geom_jitter(data = d90_data, aes(x = Background, y = .data[[value_col]]),
      width = 0.12, size = 1, alpha = 0.6, color = "black"
    ) +
    geom_segment(data = d0_data, aes(
        x = as.numeric(Background) - 0.25,
        xend = as.numeric(Background) + 0.25,
        y = .data[[value_col]],
        yend = .data[[value_col]],
        color = "D0"
      ),
      linewidth = 1.2
    ) +
    geom_text(data = p_values, aes(
        x = Background, y = star_y,
        label = signif
      ),
      size = 6
    ) +
    scale_fill_manual(
      name = "",
      values = c("D90" = "#7DAEE0")
    ) +
    scale_color_manual(
      name = "",
      values = c("D0" = "gray50")
    ) +
    guides(
      color = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    ) +
    labs(
      title = title,
      x = "Background",
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.line = element_line(linewidth = 0.8),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.key.width = unit(1.2, "cm"),
      plot.margin = margin(10, 15, 10, 10)
    ) +
    expand_limits(y = star_y * 1.05)
  
  return(p)
}
#fig1c
p_generation <- plot_relative_inhibition(
  relative_inhibition,
  "Relative_Generation_time",
  "Generation Time Normalized by YPD",
  "Boxplots of Bleomycin by Generation Time"
)
plot(p_generation)

#fig_s1b
p_lag <- plot_relative_inhibition(
  relative_inhibition,
  "Relative_Lag_time",
  "Lag Time Normalized by YPD",
  "Boxplots of Bleomycin by Lag Time"
)
plot(p_lag)

#fig_s1c
p_yield <- plot_relative_inhibition(
  relative_inhibition,
  "Relative_Yield",
  "Yield Normalized by YPD",
  "Boxplots of Bleomycin by Yield"
)
plot(p_yield)


############fig_1d, fig_s1d, and fig_s1e##########
relative_inhibition$Background <- factor(relative_inhibition$Background,
                                         levels = c("BY","MA","NA","SA","WA","WE"))
d0_rescaled <- relative_inhibition %>%
  filter(Time == "D0") %>%
  dplyr::select(Background, D0_Gen = Relative_Generation_time, D0_Lag = Relative_Lag_time, D0_Yield = Relative_Yield)

d90_rescaled <- relative_inhibition %>%
  filter(Time == "D90")

delta_rescaled <- d90_rescaled %>%
  left_join(d0_rescaled, by = "Background") %>%
  mutate(
    Delta_Gen = Relative_Generation_time - D0_Gen,
    Delta_Lag = Relative_Lag_time - D0_Lag,
    Delta_Yield = Relative_Yield - D0_Yield
  )

#fig_1d
p_gen_rescaled <- ggplot(delta_rescaled, aes(x = D0_Gen, y = Delta_Gen)) +
  geom_point(aes(color = Background), size = 3.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", size = 1, se = TRUE, fill = "gray93") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
    method = "pearson", label.x.npc = "left", label.y.npc = "top"
  ) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Correlation: generation time (YPD normalized)",
    x = "Ancestral relative inhibition (D0 Value)",
    y = expression(Delta ~ "relative inhibition (D90 - D0)"),
    color = "Backgrounds"
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)
  )

print(p_gen_rescaled)

#fig_s1d
p_Lag_rescaled <- ggplot(delta_rescaled, aes(x = D0_Lag, y = Delta_Lag)) +
  geom_point(aes(color = Background), size = 3.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", size = 1, se = TRUE, fill = "gray93") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
    method = "pearson", label.x.npc = "left", label.y.npc = "top"
  ) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Correlation: lag time (YPD normalized)",
    x = "Ancestral relative inhibition (D0 Value)",
    y = expression(Delta ~ "relative inhibition (D90 - D0)"),
    color = "Backgrounds"
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)
  )

print(p_Lag_rescaled)

#fig_s1e
p_yield_rescaled <- ggplot(delta_rescaled, aes(x = D0_Yield, y = Delta_Yield)) +
  geom_point(aes(color = Background), size = 3.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", size = 1, se = TRUE, fill = "gray93") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
    method = "pearson", label.x.npc = "left", label.y.npc = "top"
  ) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Correlation: yield (YPD normalized)",
    x = "Ancestral relative inhibition (D0 Value)",
    y = expression(Delta ~ "relative inhibition (D90 - D0)"),
    color = "Backgrounds"
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)
  )

print(p_yield_rescaled)

##save pdf
#fig1c
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/fig1C.pdf", p_generation, width = 8, height = 4, device = cairo_pdf, dpi = 300)

#fig.s1b
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/figS1B.pdf", p_lag, width = 8, height = 4, device = cairo_pdf, dpi = 300)

#fig.s1c
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/figS1C.pdf", p_yield, width = 8, height = 4, device = cairo_pdf, dpi = 300)

#fig1d
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/fig1D.pdf", p_gen_rescaled, width = 6, height = 4, device = cairo_pdf, dpi = 300)

#figS1d
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/figS1D.pdf", p_Lag_rescaled, width = 6, height = 4, device = cairo_pdf, dpi = 300)

#figS1e
ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/01_fig1_and_fig_s1/output/figS1E.pdf", p_yield_rescaled, width = 6, height = 4, device = cairo_pdf, dpi = 300)








