library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(GenomicRanges)
library(writexl)

dir <- "E:/Bleomycin_Research_Yeast_Project/04_SV_merge_vcf_to_xlsx/output/"

##step1. filter variants with an ALT depth ≥ 5.

base_dir <- "E:/Bleomycin_Research_Yeast_Project/01_variant_calling/output/03_SV_calling"
output_dir <- paste0(dir,"filter_low_depth/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#manta
all_vcf_files <- list.files(path = base_dir, pattern = "\\.vcf$", recursive = TRUE, full.names = TRUE)
manta_files <- grep("\\.manta\\.SV\\.vcf$", all_vcf_files, value = TRUE)

process_and_save_vcf <- function(vcf_path, out_dir, min_altdp = 5) {
  all_lines <- readLines(vcf_path)
  header_lines <- all_lines[grep("^#", all_lines)]
  vcf_cols_line <- header_lines[length(header_lines)]
  vcf_cols <- unlist(strsplit(vcf_cols_line, "\t"))
  sample_col_name <- tail(vcf_cols, 1) 
  message(paste("processing file:", basename(vcf_path), " | checked sample col:", sample_col_name))
  vcf_data <- fread(vcf_path, skip = "#CHROM", data.table = FALSE)
  filtered_data <- vcf_data %>%
    mutate(
      pr_val = str_extract(!!sym(sample_col_name), "(?<=:)[0-9]+,[0-9]+(?=(:|$))"),
      sr_val = str_extract(!!sym(sample_col_name), "[0-9]+,[0-9]+$")
    ) %>%
    separate(pr_val, c("PR_REF", "PR_ALT"), sep = ",", convert = TRUE, fill = "right") %>%
    separate(sr_val, c("SR_REF", "SR_ALT"), sep = ",", convert = TRUE, fill = "right") %>%
    mutate(across(c(PR_ALT, SR_ALT), ~ replace_na(.x, 0))) %>%
    mutate(DP_ALT = PR_ALT + SR_ALT) %>%
    filter(DP_ALT >= min_altdp) %>%
    dplyr::select(-c(PR_REF, PR_ALT, SR_REF, SR_ALT, DP_ALT))
  output_path <- file.path(out_dir, basename(vcf_path))
  writeLines(header_lines, output_path)
  fwrite(filtered_data, output_path, append = TRUE, sep = "\t", col.names = FALSE)
  message(paste("successful save ->", output_path))
}

for (f in manta_files) {
  process_and_save_vcf(vcf_path = f, out_dir = output_dir, min_altdp = 5)
}

#delly
delly_files <- grep("\\.delly\\.SV\\.vcf$", all_vcf_files, value = TRUE)

process_and_save_delly <- function(vcf_path, out_dir, min_altdp = 5) {
  
  all_lines <- readLines(vcf_path)
  header_lines <- all_lines[grep("^#", all_lines)]
  message(paste("processing file:", basename(vcf_path)))
  vcf_data <- fread(vcf_path, skip = "#CHROM", data.table = FALSE)
  filtered_data <- vcf_data %>%
    mutate(
      PE_val = str_extract(INFO, "(?<=PE=)[0-9]+"),
      SR_val = str_extract(INFO, "(?<=SR=)[0-9]+")
    ) %>%
    mutate(
      PE_val = as.integer(coalesce(PE_val, "0")),
      SR_val = as.integer(coalesce(SR_val, "0")),
      RSV = PE_val + SR_val
    ) %>%
    filter(RSV >= min_altdp) %>%
    dplyr::select(-c(PE_val, SR_val, RSV))
  output_path <- file.path(out_dir, basename(vcf_path))
  writeLines(header_lines, output_path)
  fwrite(filtered_data, output_path, append = TRUE, sep = "\t", col.names = FALSE)
  message(paste("successful save ->", output_path))
}

for (f in delly_files) {
  process_and_save_delly(vcf_path = f, out_dir = output_dir, min_altdp = 5)
}

##step2. Filter out variants present in the ancestor.

output_dir <- paste0(dir,"filter_ancestral/")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


base_dir <- paste0(dir,"filter_low_depth/")  
all_files <- list.files(path = base_dir, pattern = "\\.vcf$", full.names = TRUE)

manta_files <- grep("\\.manta\\.SV\\.vcf$", all_files, value = TRUE)
delly_files <- grep("\\.delly\\.SV\\.vcf$", all_files, value = TRUE)

get_vcf_end <- function(info_col, pos_col) {
  end_val <- str_extract(info_col, "(?<=END=)[0-9]+")
  ifelse(is.na(end_val), pos_col + 1, as.integer(end_val))
}

filter_ancestor_vcf <- function(vcf_path, ancestor_path, out_dir, max_gap = 100) {
  
  lines_ev <- readLines(vcf_path)
  header_ev <- lines_ev[grep("^#", lines_ev)]
  data_ev <- fread(vcf_path, skip = "#CHROM", data.table = FALSE)
  
  if (nrow(data_ev) == 0) {
    writeLines(header_ev, file.path(out_dir, basename(vcf_path)))
    return(NULL)
  }
  
  data_an <- fread(ancestor_path, skip = "#CHROM", data.table = FALSE)
  
  if (nrow(data_an) == 0) {
    writeLines(header_ev, file.path(out_dir, basename(vcf_path)))
    fwrite(data_ev, file.path(out_dir, basename(vcf_path)), append = TRUE, sep = "\t", col.names = FALSE)
    return(NULL)
  }
  
  data_ev$END <- get_vcf_end(data_ev$INFO, data_ev$POS)
  data_ev$SVTYPE <- str_extract(data_ev$INFO, "(?<=SVTYPE=)[A-Z]+")
  
  data_an$END <- get_vcf_end(data_an$INFO, data_an$POS)
  data_an$SVTYPE <- str_extract(data_an$INFO, "(?<=SVTYPE=)[A-Z]+")
  
  gr_ev <- GRanges(seqnames = data_ev$`#CHROM`, ranges = IRanges(start = data_ev$POS, end = data_ev$END))
  gr_an <- GRanges(seqnames = data_an$`#CHROM`, ranges = IRanges(start = data_an$POS, end = data_an$END))
  
  overlaps <- findOverlaps(gr_ev, gr_an, maxgap = max_gap)
  
  ev_idx <- queryHits(overlaps)
  an_idx <- subjectHits(overlaps)
  
  true_ancestor_matches <- data_ev$SVTYPE[ev_idx] == data_an$SVTYPE[an_idx]
  ev_rows_to_remove <- unique(ev_idx[true_ancestor_matches])
  
  if (length(ev_rows_to_remove) > 0) {
    filtered_data <- data_ev[-ev_rows_to_remove, ]
  } else {
    filtered_data <- data_ev
  }
  
  clean_data <- filtered_data %>% dplyr::select(-c(END, SVTYPE))
  
  output_path <- file.path(out_dir, basename(vcf_path))
  writeLines(header_ev, output_path)
  if (nrow(clean_data) > 0) {
    fwrite(clean_data, output_path, append = TRUE, sep = "\t", col.names = FALSE)
  }
}

process_multi_backgrounds <- function(files_list, software_name) {
  ancestor_files <- grep("/D0_", files_list, value = TRUE)
  evolved_files <- files_list[!files_list %in% ancestor_files]
  message(paste("\n====== start process", software_name, "file ======"))
  for (f_ev in evolved_files) {
    file_name <- basename(f_ev)
    bg_tag <- unlist(strsplit(file_name, "_"))[2]
    match_pattern <- paste0("D0_", bg_tag)
    f_an <- ancestor_files[grep(match_pattern, basename(ancestor_files))]
    if (length(f_an) == 0) {
      warning(paste("cant find sample", file_name, "ancestral [", match_pattern, "]，skip this sample"))
      next
    }
    message(paste("evolved sample:", file_name, " ---> [successful match ancestral]:", basename(f_an)))
    filter_ancestor_vcf(vcf_path = f_ev, ancestor_path = f_an, out_dir = output_dir, max_gap = 100)
  }
}
process_multi_backgrounds(manta_files, "Manta")
process_multi_backgrounds(delly_files, "Delly")

##step3. merge manta and delly

input_dir <- paste0(dir,"filter_ancestral/")
output_dir <- paste0(dir,"final_vcf/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

all_files <- list.files(path = input_dir, pattern = "\\.vcf$", full.names = TRUE)
manta_files <- grep("\\.manta\\.SV\\.vcf$", all_files, value = TRUE)
delly_files <- grep("\\.delly\\.SV\\.vcf$", all_files, value = TRUE)

get_vcf_end <- function(info_col, pos_col) {
  end_val <- str_extract(info_col, "(?<=END=)[0-9]+")
  ifelse(is.na(end_val), pos_col + 1, as.integer(end_val))
}


get_sample_id <- function(file_path) {
  file_name <- basename(file_path)
  str_extract(file_name, ".*?(?=\\.(manta|delly))")
}

manta_samples <- sapply(manta_files, get_sample_id)
delly_samples <- sapply(delly_files, get_sample_id)

all_samples <- unique(c(manta_samples, delly_samples))

no_variants_list <- c()

message(paste("sample total count:", length(all_samples)))
for (sample in all_samples) {
  
  m_file <- manta_files[manta_samples == sample]
  d_file <- delly_files[delly_samples == sample]
  
  m_data <- if(length(m_file) > 0) fread(m_file, skip = "#CHROM", data.table = FALSE) else data.frame()
  d_data <- if(length(d_file) > 0) fread(d_file, skip = "#CHROM", data.table = FALSE) else data.frame()
  
  if (length(m_file) > 0) {
    vcf_header <- readLines(m_file)[grep("^#", readLines(m_file))]
  } else if (length(d_file) > 0) {
    vcf_header <- readLines(d_file)[grep("^#", readLines(d_file))]
  }
  
  if (nrow(m_data) == 0 && nrow(d_data) == 0) {
    no_variants_list <- c(no_variants_list, sample)
    writeLines(vcf_header, file.path(output_dir, paste0(sample, ".final_union.vcf")))
    message(paste("➔ sample", sample, ": no any new variant [already add list]"))
    next
  }
  
  if (nrow(m_data) == 0) {
    final_union <- d_data %>% mutate(SOURCE = "Delly_only")
    
  } else if (nrow(d_data) == 0) {
    final_union <- m_data %>% mutate(SOURCE = "Manta_only")
    
  } else {
    m_data$END <- get_vcf_end(m_data$INFO, m_data$POS)
    m_data$SVTYPE <- str_extract(m_data$INFO, "(?<=SVTYPE=)[A-Za-z0-9]+")
    m_data$SOURCE <- "Manta_only"  
    
    d_data$END <- get_vcf_end(d_data$INFO, d_data$POS)
    d_data$SVTYPE <- str_extract(d_data$INFO, "(?<=SVTYPE=)[A-Za-z0-9]+")
    d_data$SOURCE <- "Delly_only"  
    
    gr_manta <- GRanges(seqnames = m_data$`#CHROM`, ranges = IRanges(start = m_data$POS, end = m_data$END))
    gr_delly <- GRanges(seqnames = d_data$`#CHROM`, ranges = IRanges(start = d_data$POS, end = d_data$END))
    overlaps <- suppressWarnings(findOverlaps(gr_manta, gr_delly, maxgap = 100))
    
    if (length(overlaps) > 0) {
      m_idx <- queryHits(overlaps)
      d_idx <- subjectHits(overlaps)
      type_match <- m_data$SVTYPE[m_idx] == d_data$SVTYPE[d_idx]
      
      m_shared_rows <- unique(m_idx[type_match])
      d_shared_rows <- unique(d_idx[type_match])
      
      if(length(m_shared_rows) > 0) {
        m_data$SOURCE[m_shared_rows] <- "Both"
      }
      d_unique_data <- if(length(d_shared_rows) > 0) d_data[-d_shared_rows, ] else d_data
    } else {
      d_unique_data <- d_data
    }
    
    final_union <- rbind(m_data, d_unique_data) %>% 
      arrange(`#CHROM`, POS)
    final_union <- final_union %>% dplyr::select(-c(END, SVTYPE))
  }
  
  final_union <- final_union %>%
    mutate(INFO = paste0("METHOD=", SOURCE, ";", INFO)) %>%
    dplyr::select(-SOURCE) 
  
  output_vcf <- file.path(output_dir, paste0(sample, ".final_union.vcf"))
  writeLines(vcf_header, output_vcf)
  fwrite(final_union, output_vcf, append = TRUE, sep = "\t", col.names = FALSE)
  
  message(paste("sample", sample, ": finish，mutation count:", nrow(final_union)))
}

##step4. merge to xlsx
input_dir  <- paste0(dir,"final_vcf/")
output_dir <- dir
output_file <- file.path(output_dir, "all_SV_combined_table.xlsx")

vcf_files <- list.files(path = input_dir, pattern = "\\.final_union\\.vcf$", full.names = TRUE)

message(paste("vcf count:", length(vcf_files)))

all_tables_list <- list()

for (f in vcf_files) {
  
  sample_id <- str_extract(basename(f), ".*?(?=\\.final_union)")
  
  lines <- readLines(f, n = 200) 
  data_lines <- lines[!grepl("^#", lines)]
  
  if (length(data_lines) == 0 && file.info(f)$size < 5000) {
    message(paste("sample", sample_id, "No variants found; skipping the merge."))
    next
  }
  
  vcf_data <- fread(f, skip = "#CHROM", data.table = FALSE)
  
  if (nrow(vcf_data) == 0) next
  
  
  standard_colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE_GENOTYPE")
  colnames(vcf_data) <- standard_colnames
  
  vcf_data <- vcf_data %>%
    mutate(
      SAMPLE_ID = sample_id,  
      METHOD    = str_extract(INFO, "(?<=METHOD=)[A-Za-z_]+"),
      SVTYPE    = str_extract(INFO, "(?<=SVTYPE=)[A-Za-z0-9]+"),
      END       = as.integer(str_extract(INFO, "(?<=END=)[0-9]+"))
    )
  
  
  vcf_data$END <- ifelse(is.na(vcf_data$END), vcf_data$POS + 1, vcf_data$END)
  
  final_vcf_data <- vcf_data %>%
    dplyr::select(
      SAMPLE_ID, CHROM, POS, END, SVTYPE, METHOD, 
      ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE_GENOTYPE
    )
  
  all_tables_list[[sample_id]] <- final_vcf_data
  message(paste("sample", sample_id, "add list"))
}

if (length(all_tables_list) > 0) {
  combined_table <- as.data.frame(rbindlist(all_tables_list))
  write_xlsx(combined_table, path = output_file, format_headers = TRUE)
  cat("\nfinish")
  cat("\nresult save path:", output_file)
  cat("\nmutation counts:", nrow(combined_table), "\n")
} else {
  print(" warning：no data")
}


##step5. The BAM files of the ancestor (D0) and the evolved samples were visualized in IGV locally to manually validate each SV. 
##Structural variants present in the ancestor and false positives were excluded.

