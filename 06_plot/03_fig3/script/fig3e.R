library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)
library(rtracklayer)

data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
data <- data %>%
  filter(Sample_type == "population",
         Drug == "BLM")

data1 <- data %>%
  mutate(Classified_Consequence = case_when(
    grepl("missense_variant", Consequence) ~ "Missense",
    Consequence == "synonymous_variant" ~ "Synonymous",
    Consequence == "stop_gained" ~ "Nonsense",
    Consequence %in% c("frameshift_variant", "inframe_deletion") ~ "Indel in protein coding gene",
    Consequence %in% c("upstream_gene_variant", "downstream_gene_variant", 
                       "intron_variant", "intergenic_variant") ~ "Non-coding",
    grepl("splice", Consequence) ~ "Non-coding",
    TRUE ~ "Other" 
  ))

data1$Classified_Consequence <- factor(
  data1$Classified_Consequence,
  levels = c("Non-coding", "Synonymous", "Missense", "Nonsense", "Indel in protein coding gene", "Other")
)
table(data1$Classified_Consequence)

mutation_summary <- data1 %>%
  group_by(Background, Classified_Consequence) %>%
  summarise(count = n(), .groups = "drop_last") %>% 
  mutate(percentage = (count / sum(count)) * 100) %>% 
  ungroup()

###statistical test
wide_data <- mutation_summary %>%
  dplyr::select(Background, Classified_Consequence, count) %>%  
  pivot_wider(names_from = Classified_Consequence, values_from = count, values_fill = 0)
mat <- as.matrix(wide_data[, -1])
rownames(mat) <- wide_data$Background
chisq_test <- chisq.test(mat)
chisq_test
chisq_test$expected
sum(chisq_test$expected < 1)

set.seed(123)
chisq.test(mat, simulate.p.value = TRUE, B = 100000)

strains <- rownames(mat)
n <- length(strains)
pairwise_results <- data.frame(
  group1 = character(),
  group2 = character(),
  X_squared = numeric(),
  p_value = numeric(),
  p_adjusted = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    sub_mat <- mat[c(i, j), ]
    
    if (sum(sub_mat[1, ]) == 0 || sum(sub_mat[2, ]) == 0) next
    set.seed(123)
    test_result <- chisq.test(sub_mat, simulate.p.value = TRUE, B = 100000)
    pairwise_results <- rbind(pairwise_results, data.frame(
      group1 = strains[i],
      group2 = strains[j],
      X_squared = test_result$statistic,
      p_value = test_result$p.value
    ))
  }
}
pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, method = "BH")
pairwise_results


#null model
background_params <- tribble(
  ~id, ~total_SNV, ~total_indel, ~total_genome_length,
  "BY", 329, 39, 12157149,
  "WE", 210, 40, 11813188, 
  "WA", 59, 9, 11958041,
  "NA", 161, 19, 11910474,
  "SA", 130, 14, 11881147,
  "MA", 27, 3, 11747346
)

base_gff_path <- "E:/Bleomycin_Research_Yeast_Project/02_SNV_INDEL_merge_vep/input/gff/"
expected_all <- map_dfr(1:nrow(background_params), function(i) {
  id <- background_params$id[i]
  total_SNV <- background_params$total_SNV[i]
  total_indel <- background_params$total_indel[i]
  total_genome_length <- background_params$total_genome_length[i]
  path <- paste0(base_gff_path, id, ".gff.gz")
  gff <- import(path) 
  cds <- gff[gff$type == "CDS"]
  cds_reduce <- reduce(cds)
  total_cds_length <- sum(width(cds_reduce))
  r_coding <- total_cds_length / total_genome_length
  r_non_coding <- 1 - r_coding
  expected_coding_SNV <- total_SNV * r_coding
  expected_noncoding_SNV <- total_SNV * r_non_coding
  expected_syn <- expected_coding_SNV * 0.246
  expected_missense <- expected_coding_SNV * 0.696
  expected_nonsense <- expected_coding_SNV * 0.058
  expected_coding_indel <- total_indel * r_coding
  expected_noncoding_indel <- total_indel * r_non_coding
  data.frame(
    strain = id,
    `Non-coding` = expected_noncoding_SNV + expected_noncoding_indel,
    Synonymous = expected_syn,
    Missense = expected_missense,
    Nonsense = expected_nonsense,
    `Indel in protein coding gene` = expected_coding_indel,
    Type = "Expected"
  )
})

print(expected_all)

expected_all <- expected_all %>%
  rename('Non-coding' = Non.coding,
         'Indel in protein coding gene' = Indel.in.protein.coding.gene,
         Background = strain)

wide_data$type <- "Observed"
wide_data <- wide_data %>%
  rename(Type = type)

plot_data <- bind_rows(
  wide_data,
  expected_all
)


#statistical test
set.seed(123)

results <- lapply(unique(plot_data$Background), function(s) {
  
  tmp <- filter(plot_data, Background == s)
  
  obs <- tmp %>%
    filter(Type == "Observed") %>%
    dplyr::select(-Background, -Type) %>%
    unlist(use.names = FALSE)
  
  exp <- tmp %>%
    filter(Type == "Expected") %>%
    dplyr::select(-Background, -Type) %>%
    unlist(use.names = FALSE)
  
  test <- chisq.test(
    x = obs,
    p = exp / sum(exp),
    simulate.p.value = TRUE,
    B = 100000
  )
  
  data.frame(
    Background = s,
    statistic = unname(test$statistic),
    p.value = test$p.value
  )
})

results <- bind_rows(results) %>%
  mutate(
    p.adjusted = p.adjust(p.value, method = "BH")
  )
results

#####plot
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = -c(Background, Type),
    names_to = "Classified_Consequence",
    values_to = "Count"
  ) %>%
  mutate(Count = as.numeric(Count))


plot_percent <- plot_data_long %>%
  group_by(Background, Type) %>%
  mutate(
    total_count = sum(Count, na.rm = TRUE),
    proportion = Count / total_count
  ) %>%
  ungroup()


totals <- plot_percent %>%
  dplyr::select(Background, Type, total_count) %>%
  distinct() %>%
  mutate(total_count = round(total_count, 0))


plot_percent <- plot_percent %>%
  mutate(
    Background = factor(
      Background,
      levels = c("BY", "WE", "WA", "NA", "SA", "MA")
    ),
    Classified_Consequence = factor(
      Classified_Consequence,
      levels = c(
        "Non-coding",
        "Synonymous",
        "Missense",
        "Nonsense",
        "Indel in protein coding gene"
      )
    ),
    Type = factor(
      Type,
      levels = c("Observed", "Expected")
    )
  )

totals <- totals %>%
  mutate(
    Background = factor(
      Background,
      levels = c("BY", "WE", "WA", "NA", "SA", "MA")
    ),
    Type = factor(
      Type,
      levels = c("Observed", "Expected")
    )
  )

my_spectrum_colors <- c(
  "Non-coding"                   = "gray", 
  "Synonymous"                   = "#5BB381", 
  "Missense"                     = "#E69F00", 
  "Nonsense"                     = "#E06666", 
  "Indel in protein coding gene" = "#8E7CC3"  
)
o_e <- ggplot(plot_percent, aes(x = Type, y = proportion, fill = Classified_Consequence)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(proportion >= 0.05,
                       sprintf("%.2f", proportion),
                       "")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  geom_text(
    data = totals,
    aes(x = Type, y = 1.02, label = total_count),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  )+
  facet_wrap(~Background, nrow = 1) +
  scale_fill_manual(values = my_spectrum_colors) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.2),
    labels = sprintf("%.1f", seq(0, 1, 0.2)),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion of variants",
    fill = "Mutation type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20, "pt")
  )
plot(o_e)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3E.pdf", width = 11, height = 6, dpi = 300, device = cairo_pdf)









