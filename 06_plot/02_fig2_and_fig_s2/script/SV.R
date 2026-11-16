library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(patchwork)


##fig2B

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/input/SV_stats.xlsx", sheet = 1)
data_long <- data %>%
  pivot_longer(cols = -Background, names_to = "SV_Type", values_to = "Count") %>%
  mutate(Count = as.numeric(Count)) %>%
  filter(!is.na(Count) & Count > 0) 

inner_ring <- data_long %>%
  group_by(SV_Type) %>%
  summarise(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(
    fraction = Total / sum(Total),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    label_pos = (ymax + ymin) / 2
  )

outer_ring <- data_long %>%
  group_by(SV_Type) %>%
  mutate(fraction_sub = Count / sum(Count)) %>% 
  ungroup() %>%
  left_join(inner_ring %>% dplyr::select(SV_Type, ymin, ymax, Total), by = "SV_Type") %>%
  group_by(SV_Type) %>%
  mutate(
    sub_fraction = (Count / sum(Count)) * (ymax - ymin),
    sub_ymax = ymin + cumsum(sub_fraction),
    sub_ymin = c(ymin[1], head(sub_ymax, n=-1)),
    outer_label_pos = (sub_ymax + sub_ymin) / 2
  ) %>%
  ungroup() %>%
  mutate(
    angle = 90 - (outer_label_pos * 360),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

sv_colors <- c(
  "DEL"     = "#E41A1C",
  "DUP"     = "#377EB8", 
  "INS"     = "#4DAF4A",
  "INV"     = "#984EA3", 
  "TRA" = "#FF7F00"  
)

bg_levels <- unique(outer_ring$Background)
num_bgs <- length(bg_levels)

outer_colors <- c()
for(sv in names(sv_colors)) {
  sv_grad <- colorRampPalette(c("white", sv_colors[sv]))(num_bgs + 2)[3:(num_bgs + 2)]
  names(sv_grad) <- paste(sv, bg_levels, sep = "_")
  outer_colors <- c(outer_colors, sv_grad)
}

outer_ring <- outer_ring %>% 
  mutate(color_key = paste(SV_Type, Background, sep = "_"))

p1 <- ggplot() +
  geom_rect(data = outer_ring, 
            aes(xmin = 3.2, xmax = 4.2, ymin = sub_ymin, ymax = sub_ymax, fill = color_key), 
            color = "white", size = 0.3) +
  geom_rect(data = inner_ring, 
            aes(xmin = 0, xmax = 3.2, ymin = ymin, ymax = ymax, fill = SV_Type), 
            color = "white", size = 0.5) +
  geom_text(data = inner_ring, 
            aes(x = 2, y = label_pos, 
                label = paste0(SV_Type, "\n", Total)), 
            size = 3.8, fontface = "bold", color = "white") +
  geom_text(data = outer_ring, 
            aes(x = 3.7, y = outer_label_pos, 
                label = paste0(Background, "\n", Count)), 
            size = 3.0, fontface = "bold", color = "black") + 
  coord_polar(theta = "y") +
  xlim(0, 4.6) +  
  scale_fill_manual(values = c(sv_colors, outer_colors), guide = "none") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = -1)) +
  labs(title = "BLM_SV")

print(p1)


data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/input/SV_stats.xlsx", sheet = 2)
data_long <- data %>%
  pivot_longer(cols = -Background, names_to = "SV_Type", values_to = "Count") %>%
  mutate(Count = as.numeric(Count)) %>%
  filter(!is.na(Count) & Count > 0) 

inner_ring <- data_long %>%
  group_by(SV_Type) %>%
  summarise(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(
    fraction = Total / sum(Total),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    label_pos = (ymax + ymin) / 2
  )

outer_ring <- data_long %>%
  group_by(SV_Type) %>%
  mutate(fraction_sub = Count / sum(Count)) %>% 
  ungroup() %>%
  left_join(inner_ring %>% dplyr::select(SV_Type, ymin, ymax, Total), by = "SV_Type") %>%
  group_by(SV_Type) %>%
  mutate(
    sub_fraction = (Count / sum(Count)) * (ymax - ymin),
    sub_ymax = ymin + cumsum(sub_fraction),
    sub_ymin = c(ymin[1], head(sub_ymax, n=-1)),
    outer_label_pos = (sub_ymax + sub_ymin) / 2
  ) %>%
  ungroup() %>%
  mutate(
    angle = 90 - (outer_label_pos * 360),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

sv_colors <- c(
  "DEL"     = "#E41A1C",
  "DUP"     = "#377EB8", 
  "INS"     = "#4DAF4A",
  "INV"     = "#984EA3", 
  "TRA" = "#FF7F00"  
)

bg_levels <- unique(outer_ring$Background)
num_bgs <- length(bg_levels)

outer_colors <- c()
for(sv in names(sv_colors)) {
  sv_grad <- colorRampPalette(c("white", sv_colors[sv]))(num_bgs + 2)[3:(num_bgs + 2)]
  names(sv_grad) <- paste(sv, bg_levels, sep = "_")
  outer_colors <- c(outer_colors, sv_grad)
}

outer_ring <- outer_ring %>% 
  mutate(color_key = paste(SV_Type, Background, sep = "_"))

p2 <- ggplot() +
  geom_rect(data = outer_ring, 
            aes(xmin = 3.2, xmax = 4.2, ymin = sub_ymin, ymax = sub_ymax, fill = color_key), 
            color = "white", size = 0.3) +
  geom_rect(data = inner_ring, 
            aes(xmin = 0, xmax = 3.2, ymin = ymin, ymax = ymax, fill = SV_Type), 
            color = "white", size = 0.5) +
  geom_text(data = inner_ring, 
            aes(x = 2, y = label_pos, 
                label = paste0(SV_Type, "\n", Total)), 
            size = 3.8, fontface = "bold", color = "white") +
  geom_text(data = outer_ring, 
            aes(x = 3.7, y = outer_label_pos, 
                label = paste0(Background, "\n", Count)), 
            size = 3.0, fontface = "bold", color = "black") + 
  coord_polar(theta = "y") +
  xlim(0, 4.6) +  
  scale_fill_manual(values = c(sv_colors, outer_colors), guide = "none") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = -1)) +
  labs(title = "YPD_SV")

print(p2)

p <- p1 / p2

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/output/fig2B.pdf", plot = p, width = 4, height = 8, dpi = 300)





