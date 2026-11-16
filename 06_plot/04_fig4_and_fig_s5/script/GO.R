library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)
library(org.Sc.sgd.db)
library(reshape2)
library(enrichplot)
library(clusterProfiler)

#input data
data <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/03_fig3/input/Table_S3.xlsx", sheet = 1)
data <- data %>%
  filter(Sample_type == "population", Drug == "BLM")

drug_high_moderate <- data %>%
  filter(Impact %in% c("HIGH", "MODERATE")) %>%
  group_by(Drug) %>%
  summarise(genes = list(unique(na.omit(Gene_ID))), 
            gene_count = length(unique(na.omit(Gene_ID))))

drug_gene_list <- drug_high_moderate$genes
names(drug_gene_list) <- drug_high_moderate$Drug

drug_go_results <- list()

for (drug_name in names(drug_gene_list)) {
  genes <- drug_gene_list[[drug_name]]
  
  cat("Analyzing", drug_name, "-", length(genes), "genes\n")
  
  # BP
  ego_bp <- enrichGO(
    gene          = genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  # MF
  ego_mf <- enrichGO(
    gene          = genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF", 
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  # CC
  ego_cc <- enrichGO(
    gene          = genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  drug_go_results[[drug_name]] <- list(
    BP = ego_bp,
    MF = ego_mf,
    CC = ego_cc
  )
}
###save result
for (Drug in names(drug_go_results)) {
  write.csv(drug_go_results[[Drug]]$BP@result, 
            paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/", Drug, "_GO_BP.csv"), row.names = FALSE)
  write.csv(drug_go_results[[Drug]]$MF@result,
            paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/", Drug, "_GO_MF.csv"), row.names = FALSE)
  write.csv(drug_go_results[[Drug]]$CC@result,
            paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/", Drug, "_GO_CC.csv"), row.names = FALSE)
}

sig_threshold <- 0.05

for (Drug in names(drug_go_results)) {
  message("processing: ", Drug)  
  
  #BP
  # extract drug
  bp_result <- drug_go_results[[Drug]]$BP@result
  # p.adjust < 0.05
  if (nrow(bp_result) > 0 && any(bp_result$p.adjust < sig_threshold)) {
    sig_bp <- bp_result[bp_result$p.adjust < sig_threshold, ]
    p_bp <- dotplot(
      drug_go_results[[Drug]]$BP, 
      showCategory = min(15, nrow(sig_bp)),  
      orderBy = "x", 
      color = "p.adjust"
    ) +
      ggtitle(paste(Drug, "- GO: Biological Process"))
    
    print(p_bp)
    # save figs4A
    ggsave(
      paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/fig_s5A_", Drug, "_GO_BP_padj", sig_threshold, ".pdf"), 
      plot = p_bp, 
      width = 8, 
      height = 10, 
      device = cairo_pdf
    )
  } else {
    if (nrow(bp_result) == 0) {
      message(paste("  -", Drug, "no BP result"))
    } else {
      message(paste("  -", Drug, "No significantly enriched BP were identified (minimum adjusted P-value = ", round(min(bp_result$p.adjust), 3), "）"))
    }
  }
  
  # MF
  mf_result <- drug_go_results[[Drug]]$MF@result
  if (nrow(mf_result) > 0 && any(mf_result$p.adjust < sig_threshold)) {
    sig_mf <- mf_result[mf_result$p.adjust < sig_threshold, ]
    p_mf <- dotplot(
      drug_go_results[[Drug]]$MF, 
      showCategory = min(15, nrow(sig_mf)), 
      orderBy = "x", 
      color = "p.adjust"
    ) +
      ggtitle(paste(Drug, "- GO: Molecular Function "))
    
    print(p_mf) #figs5B
    ggsave(
      paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/fig_s5B_", Drug, "_GO_MF_padj", sig_threshold, ".pdf"), 
      plot = p_mf, 
      width = 8, 
      height = 10, 
      device = cairo_pdf
    )
  } else {
    if (nrow(mf_result) == 0) {
      message(paste("  -", Drug, "no MF result"))
    } else {
      message(paste("  -", Drug, "No significantly enriched MF were identified (minimum adjusted P-value = ", round(min(mf_result$p.adjust), 3), "）"))
    }
  }
  
  #CC
  cc_result <- drug_go_results[[Drug]]$CC@result
  if (nrow(cc_result) > 0 && any(cc_result$p.adjust < sig_threshold)) {
    sig_cc <- cc_result[cc_result$p.adjust < sig_threshold, ]
    p_cc <- dotplot(
      drug_go_results[[Drug]]$CC, 
      showCategory = min(15, nrow(sig_cc)), 
      orderBy = "x", 
      color = "p.adjust"
    ) +
      ggtitle(paste(Drug, "- GO: Cellular Component"))
    
    print(p_cc) #fig4E
    ggsave(
      paste0("E:/Bleomycin_Research_Yeast_Project/06_plot/04_fig4_and_fig_s5/output/fig_4E_", Drug, "_GO_CC_padj", sig_threshold, ".pdf"), 
      plot = p_cc, 
      width = 8, 
      height = 10, 
      device = cairo_pdf
    )
  } else {
    if (nrow(cc_result) == 0) {
      message(paste("  -", Drug, "no CC result"))
    } else {
      message(paste("  -", Drug, "No significantly enriched MF were identified (minimum adjusted P-value = ", round(min(cc_result$p.adjust), 3), "）"))
    }
  }
}
