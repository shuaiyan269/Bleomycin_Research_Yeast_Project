
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(rtracklayer)

data_cn <-read_tsv("E:/Bleomycin_Research_Yeast_Project/05_aneuploidy_analysis/output/copy_number.txt")
strain_order <- c("BY", "WE", "WA", "NA", "SA", "MA")
chromosome_order <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", 
                      "XI", "XII", "XIII", "XIV", "XV", "XVI")

data_cn_long <- data_cn %>%
  pivot_longer(
    cols = -sample,
    names_to = "chromosome",
    values_to = "value"
  ) %>%
  mutate(
    chromosome = gsub("chr(.*)_ratio", "\\1", chromosome),
    # extract chromosome name from chr()_ratio
    strain = str_split(sample, "_", simplify = TRUE)[, 2],
    Date = str_split(sample, "_", simplify = TRUE)[, 1],
    # extract strain name from sample 
    drugname = case_when(
      str_count(sample, "_") < 2 ~ NA_character_,
      str_count(sample, "_") == 3 ~ str_split(sample, "_", simplify = TRUE)[, 3],
      str_count(sample, "_") == 4 ~ "Bl_clones",  
      TRUE ~ NA_character_
      ),
    # 
    strain = factor(strain, levels = strain_order),
    chromosome = factor(chromosome, levels = chromosome_order)
  ) %>%
  # arrange 
  arrange(strain, drugname) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

view(data_cn_long)


#fig2a
filtered_data1 <- data_cn_long %>%
  filter(drugname == "Bl") 
filtered_data2 <- data_cn_long %>%
  filter(Date == "D0")
filtered_data <- rbind(filtered_data1, filtered_data2) %>%
  arrange(strain, Date, sample) %>% 
  mutate(sample = factor(sample, levels = unique(sample))) 
plot_Bl <- ggplot(data = filtered_data,
       aes(x = sample, y = chromosome, fill = value)) +
  geom_tile(color = "gray60", size = 0.3) +
  geom_text(aes(label = sprintf("%.1f", value)),  
            color = "gray30",                     
            size = 2,                            
            fontface = "bold"                   
  )+
  scale_fill_gradientn(
    colours = c("darkblue", "#87CEFA", "white", "#FB9376", "#96170B"),
    oob = scales::squish,
    values = scales::rescale(c(0, 1.5, 2, 3, 4)),
    limits = c(0, 4)
  ) +
  labs(x = "",
       y = "",
       title = "Copy Number - Bl"  # title 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )+
  coord_fixed()

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/output/fig2A.pdf", plot_Bl, width = 8, height = 6, device = cairo_pdf)


#fig_s2a
filtered_data1 <- data_cn_long %>%
  filter(drugname == "YPD") 
filtered_data2 <- data_cn_long %>%
  filter(Date == "D0")
filtered_data <- rbind(filtered_data1, filtered_data2) %>%
  arrange(strain, Date, sample) %>% 
  mutate(sample = factor(sample, levels = unique(sample))) 
plot_YPD <- ggplot(data = filtered_data,
                  aes(x = sample, y = chromosome, fill = value)) +
  geom_tile(color = "gray60", size = 0.3) +
  geom_text(aes(label = sprintf("%.1f", value)),  
            color = "gray30",                     
            size = 2,                            
            fontface = "bold"                   
  )+
  scale_fill_gradientn(
    colours = c("darkblue", "#87CEFA", "white", "#FB9376", "#96170B"),
    oob = scales::squish,
    values = scales::rescale(c(0, 1.5, 2, 3, 4)),
    limits = c(0, 4)
  ) +
  labs(x = "",
       y = "",
       title = "Copy Number - YPD"  # title 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )+
  coord_fixed()

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/output/figS2A.pdf", plot_YPD, width = 8, height = 6, device = cairo_pdf)

#fig_s2b
filtered_data1 <- data_cn_long %>%
  filter(drugname == "Bl_clones") 
filtered_data2 <- data_cn_long %>%
  filter(Date == "D0")
filtered_data <- rbind(filtered_data1, filtered_data2) %>%
  arrange(strain, Date, sample) %>% 
  mutate(sample = factor(sample, levels = unique(sample))) 
plot_clones <- ggplot(data = filtered_data,
                   aes(x = sample, y = chromosome, fill = value)) +
  geom_tile(color = "gray60", size = 0.3) +
  geom_text(aes(label = sprintf("%.1f", value)),  
            color = "gray30",                     
            size = 2,                            
            fontface = "bold"                   
  )+
  scale_fill_gradientn(
    colours = c("darkblue", "#87CEFA", "white", "#FB9376", "#96170B"),
    oob = scales::squish,
    values = scales::rescale(c(0, 1.5, 2, 3, 4)),
    limits = c(0, 4)
  ) +
  labs(x = "",
       y = "",
       title = "Copy Number - clones"  # title 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )+
  coord_fixed()

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/output/figS2B.pdf", plot_clones, width = 8, height = 6, device = cairo_pdf)



##fig2C and fig2F

input_dir <- "E:/Bleomycin_Research_Yeast_Project/05_aneuploidy_analysis/output/10kb_depth" 
file_list <- list.files(path = input_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

print(basename(file_list))

combined_data <- map_dfr(file_list, function(file) {
  cat("sample:", basename(file), "\n")
  fread(file, col.names = c("sample", "chr", "window", "depth"))  
})
combined_data <- combined_data %>%
  mutate(
    strain = str_split(sample, "_", simplify = TRUE)[, 2]
  )

combined_data <- combined_data %>%
  filter(!chr %in% c("chrMT")) %>%
  mutate(chr = factor(chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", 
                                      "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")))


data_median_depth <-read_tsv("E:/Bleomycin_Research_Yeast_Project/01_variant_calling/output/01_mapping/all_samples.coverage_summary.txt")
data_median_simplified <- data_median_depth %>%
  dplyr::select(sample, median_depth)
combined_data <- combined_data %>%
  left_join(data_median_simplified, by = "sample")

combined_data <- combined_data %>%
  mutate(
    depth_norm = (depth / median_depth) * 2
  )


id <- "MA" #Run sequentially for BY, WE, WA, NA, and MA.


data <- combined_data %>%
  filter(strain == id)
chr_lengths <- data %>%
  group_by(chr) %>%
  summarise(n_windows = n_distinct(window))

data <- data %>%
  left_join(chr_lengths, by = "chr")

##centromere
gff_file <- paste0("E:/Bleomycin_Research_Yeast_Project/02_SNV_INDEL_merge_vep/input/gff/", id, ".gff.gz")  #
gff_data <- import.gff(gff_file)
centromeres <- gff_data[grep("centromere", gff_data$type, ignore.case = TRUE), ]
centromere_coords <- as.data.frame(centromeres) %>%
  dplyr::select(seqnames, start, end, strand, type) %>%
  filter(seqnames %in% unique(data$chr))
colnames(centromere_coords)[1] <- "chr"

unique_samples <- unique(data$sample)  

n_samples <- length(unique_samples)
pdf_height <- 2
pdf_width <- 25 

pdf(
  file = paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/02_fig2_and_fig_s2/output/10kb_depth_", id, ".pdf"),  #
  width = pdf_width,                   
  height = pdf_height,                  
  onefile = TRUE                        
)

for (sample_name in unique_samples) {
  sample_data <- data %>% 
    filter(sample == sample_name) %>%  
    mutate(chr = factor(chr, levels = c("chrI", "chrII", "chrIII", "chrIV", 
                                        "chrV", "chrVI", "chrVII", "chrVIII", 
                                        "chrIX", "chrX", "chrXI", "chrXII", 
                                        "chrXIII", "chrXIV", "chrXV", "chrXVI")))
  sample_data <- sample_data %>%
    mutate(
      chr_parity = case_when(
        chr %in% c("chrI", "chrIII", "chrV", "chrVII", "chrIX", "chrXI", "chrXIII", "chrXV") ~ "odd",
        chr %in% c("chrII", "chrIV", "chrVI", "chrVIII", "chrX", "chrXII", "chrXIV", "chrXVI") ~ "even",
        TRUE ~ "other"
      )
    )
  
  p <- ggplot(sample_data, aes(x = window, y = depth_norm)) +
    geom_point(aes(color = chr_parity), alpha = 0.6, size = 0.8) +  
    scale_color_manual(values = c(odd = "blue", even = "dodgerblue"))+
    geom_hline(yintercept = 2, color = "black", linewidth = 0.2, linetype = "dashed")+
    geom_vline(data = centromere_coords, 
               aes(xintercept = (start+end)/20000),  
               color = "gray40", linetype = "dashed", linewidth = 0.2) +
    scale_y_continuous(
      breaks = c(0, 1, 2, 3, 4, 5, 6),
      limits = c(0, 6) ,
      oob = scales::squish 
    )+
    labs(
      title = sample_name, 
      x = "Genomic position (10 kb windows)",
      y = "Copy number") +
    theme_minimal() +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +  
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5), 
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 8),  
      plot.title = element_text(size = 10, face = "bold"),  
      #plot.subtitle = element_text(size = 7, color = "gray50"), 
      plot.margin = margin(10, 10, 10, 10, "pt"),  
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "none"
    )
  print(p)
}

dev.off()











