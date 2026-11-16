
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggpubr)
library(readxl)
library(patchwork)

#input growth data
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

delta_rescaled <- delta_rescaled %>%
  rename(Populations = Population)

##input mutation rate data
data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/mutation_rate.xlsx", sheet = 1)
data_clean <- data %>%
  mutate(
    Sample_Name = Samples 
  ) %>%
  separate(Sample_Name, into = c("Time_Tag", "Background", "Drug", "Populations"), sep = "_") %>%
  mutate(Populations = as.character(Populations)) %>%
  dplyr::select(Background, Populations, Drug, InDel, SNV, mutation_rate)

data_clean <- data_clean %>%
  filter(Drug == "BLM")

mutation_vs_evolution <- delta_rescaled %>%
  mutate(Populations = as.character(Populations)) %>%
  inner_join(data_clean, by = c("Background", "Populations"))

mutation_vs_d0 <- d0_rescaled %>%
  inner_join(data_clean, by = "Background")


create_corr_plot <- function(df, x_var, y_var, title_text, x_label, y_label) {
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = Background), size = 3.5, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", linewidth = 1, se = TRUE, fill = "gray93") +
    stat_cor(
      aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
      method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title_text,
      x = x_label,
      y = y_label,
      color = "Background"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.margin = margin(10, 10, 10, 10, "pt"),
      legend.position = "right",
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(color = "black", size = 10)
    )
}


x_lab_evo <- "Mutation Rate (per base per generation)"
y_lab_mut <- "D90 Mutation Rate"


p1 <- create_corr_plot(mutation_vs_evolution, "mutation_rate", "Delta_Gen", 
                       "Mutation Rate vs Delta Gen", x_lab_evo, expression(Delta ~ "relative generation time"))
p2 <- create_corr_plot(mutation_vs_evolution, "mutation_rate", "Delta_Lag", 
                       "Mutation Rate vs Delta Lag", x_lab_evo, expression(Delta ~ "relative lag time"))
p3 <- create_corr_plot(mutation_vs_evolution, "mutation_rate", "Delta_Yield", 
                       "Mutation Rate vs Delta Yield", x_lab_evo, expression(Delta ~ "relative yield"))


p4 <- create_corr_plot(mutation_vs_d0, "D0_Gen", "mutation_rate", 
                       "Ancestral Gen vs Mutation Rate", "Ancestral relative generation time", y_lab_mut)
p5 <- create_corr_plot(mutation_vs_d0, "D0_Lag", "mutation_rate", 
                       "Ancestral Lag vs Mutation Rate", "Ancestral relative lag time", y_lab_mut)
p6 <- create_corr_plot(mutation_vs_d0, "D0_Yield", "mutation_rate", 
                       "Ancestral Yield vs Mutation Rate", "Ancestral relative yield", y_lab_mut)


final_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold', size = 12))
print(final_plot)
ggsave(
  filename = "E:/Bleomycin_Research_Yeast_Project/07_other/output/correlation_mutation_rate_with_growth_data.pdf", 
  plot = final_plot, 
  width = 20, 
  height = 10, 
  dpi = 300, 
  device = cairo_pdf
)









