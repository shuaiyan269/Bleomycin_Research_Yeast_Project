library(MutationalPatterns)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(readxl)


process_background <- function(bg_name, combined_data, ref_genome, ref_lengths, all_96_types) {
  
  cat("Processing background:", bg_name, "\n")
  filtered_data <- combined_data %>%
    filter(Background == bg_name) %>%
    filter(Variant_type == "SNV") %>%
    filter(nchar(REF) == 1 & nchar(ALT) == 1)
  if (nrow(filtered_data) == 0) {
    return(matrix(0, nrow = length(all_96_types), ncol = 0, dimnames = list(all_96_types, NULL)))
  }
  
  mutations <- filtered_data %>%
    dplyr::select(chr = "CHROM", pos = "POS", ref = "REF", alt = "ALT", sample = "Sample") %>%
    filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(alt))
  granges_mut <- GRanges(
    seqnames = mutations$chr,
    ranges = IRanges(start = mutations$pos, end = mutations$pos),
    ref = mutations$ref,
    alt = mutations$alt,
    sample = mutations$sample
  )
  get_custom_context <- function(granges, ref_lengths, ref_genome) {
    granges_extended <- resize(granges, width = 3, fix = "center")
    chr_names <- as.character(seqnames(granges_extended))
    valid_ranges <- start(granges_extended) >= 1 & end(granges_extended) <= ref_lengths[chr_names]
    valid_ranges[is.na(valid_ranges)] <- FALSE 
    granges_extended <- granges_extended[valid_ranges]
    context_seqs <- getSeq(ref_genome, granges_extended)
    context_chr <- as.character(context_seqs)
    context_full <- rep(NA, length(granges))
    context_full[valid_ranges] <- context_chr
    return(context_full)
  }
  mut_context <- get_custom_context(granges_mut, ref_lengths, ref_genome)
  valid_context_idx <- !is.na(mut_context) & nchar(mut_context) == 3 & substr(mut_context, 2, 2) == granges_mut$ref  
  granges_mut <- granges_mut[valid_context_idx]
  mut_context <- mut_context[valid_context_idx]

  if (length(granges_mut) == 0) {
    return(matrix(0, nrow = length(all_96_types), ncol = 0, dimnames = list(all_96_types, NULL)))
  }
  
  comp_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "N" = "N")
  raw_five  <- substr(mut_context, 1, 1)
  raw_ref   <- granges_mut$ref
  raw_alt   <- granges_mut$alt
  raw_three <- substr(mut_context, 3, 3)
  
  need_convert <- raw_ref %in% c("A", "G")
  
  final_ref   <- raw_ref
  final_alt   <- raw_alt
  final_five  <- raw_five
  final_three <- raw_three
  
  final_ref[need_convert]   <- comp_map[raw_ref[need_convert]]
  final_alt[need_convert]   <- comp_map[raw_alt[need_convert]]
  final_five[need_convert]  <- comp_map[raw_three[need_convert]]  
  final_three[need_convert] <- comp_map[raw_five[need_convert]]   
  
  mut_types <- paste0(final_five, "[", final_ref, ">", final_alt, "]", final_three)
  mut_table <- table(mut_types, granges_mut$sample)
  mut_mat <- as.matrix(mut_table)
  
  complete_mut_mat <- matrix(0, 
                             nrow = length(all_96_types), 
                             ncol = ncol(mut_mat),
                             dimnames = list(all_96_types, colnames(mut_mat)))
  
  common_types <- intersect(rownames(mut_mat), all_96_types)
  complete_mut_mat[common_types, ] <- mut_mat[common_types, ]
  
  return(complete_mut_mat)
}


combined_data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
combined_data <- combined_data %>%
  filter(Sample_type == "population", Drug == "BLM")

create_96_mutation_types <- function() {
  substitutions <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")  
  flanking_bases <- c("A", "C", "G", "T")
  all_types <- c()
  for (five_prime in flanking_bases) {
    for (sub in substitutions) {
      for (three_prime in flanking_bases) {
        all_types <- c(all_types, paste0(five_prime, "[", sub, "]", three_prime))
      }
    }
  }
  return(all_types)
}
all_96_types <- create_96_mutation_types()

backgrounds <- c("BY", "WE", "WA", "NA", "SA", "MA")
mat_list <- list()


ref_dir <- "E:/Bleomycin_Research_Yeast_Project/02_SNV_INDEL_merge_vep/input/Reference_Genomes/"

for (bg in backgrounds) {
  
  ref_path <- paste0(ref_dir, bg, ".tidy.lite.fa")
  ref_genome <- readDNAStringSet(ref_path)
  ref_lengths <- width(ref_genome)
  names(ref_lengths) <- names(ref_genome)
  mat <- process_background(bg, combined_data, ref_genome, ref_lengths, all_96_types)
  if (ncol(mat) > 0) {
    mat_list[[bg]] <- mat
  }
}

combined_mut_mat <- do.call(cbind, mat_list)

print(dim(combined_mut_mat))
print(colSums(combined_mut_mat))




write.csv(combined_mut_mat, 
          file = "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/BLM_combined_96_mutation_matrix.csv",
          row.names = TRUE)


plot_matrix <- combined_mut_mat %>% 
  as.matrix()    
print(head(rownames(plot_matrix)))

p <- plot_96_profile(plot_matrix) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
  )
print(p)


target_order <- c("BY", "WE", "WA", "NA", "SA", "MA")
df_long <- as.data.frame(plot_matrix) %>%
  rownames_to_column(var = "mutation_type") %>%
  pivot_longer(cols = -mutation_type, names_to = "sample", values_to = "count") %>%
  mutate(background = sub("^[^_]+_([^_]+)_.*$", "\\1", sample)) %>%
  mutate(background = factor(background, levels = target_order))
merged_df <- df_long %>%
  group_by(mutation_type, background) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = background, values_from = total_count, names_expand = TRUE)

merged_matrix <- merged_df %>%
  column_to_rownames(var = "mutation_type") %>%
  as.matrix()

merged_matrix <- merged_matrix[all_96_types, , drop = FALSE]

available_order <- intersect(target_order, colnames(merged_matrix))
merged_matrix <- merged_matrix[, available_order, drop = FALSE]

print(colnames(merged_matrix))
write.csv(merged_matrix, 
          file = "E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/BLM_background_combined_96_mutation_matrix.csv",
          row.names = TRUE)
p <- plot_96_profile(merged_matrix) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
  )
print(p)

ggsave("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/output/fig3H.pdf", 
       plot = p, 
       width = 15, 
       height = 8, 
       device = "pdf")










