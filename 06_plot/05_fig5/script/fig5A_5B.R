library(Biostrings)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)
library(g3viz)
library(trackViewer)
library(GenomicRanges)

#pma1
#fig5A
PMA1 <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/input/data_fig5AB.xlsx", sheet = 2)
sample.gr <- GRanges("PMA1", 
                     IRanges(
                       start = PMA1$SNP, width = 1, names = PMA1$names)
)
features <- GRanges("PMA1", IRanges(c(78, 176, 375), 
                                    width=c(56, 100, 271),
                                    names=c("Cation transporter/ATPase, N-terminus", "P-type ATPase actuator domain", "haloacid dehalogenase-like hydrolase"
                                    )))


features$fill <- c("#FF6B6B", "#4ECDC4", "#344e99")
bg_col <- c(
  "BY" = "#F2B56E",
  "WE" = "#DDA0DD",
  "WA" = "#71BBED",
  "NA" = "#FBE79E",
  "SA" = "#84D8D7",
  "MA" = "#436289"
)


sample.gr$color <- bg_col[PMA1$backgrounds]

features$height <- rep(0.02, 1)
features$featureLayerID <- c("upper_layer")  
xaxis <- c(1, 200, 400, 600, 800, 918)
yaxis <- 6
legend <- list(
  list(
    labels = names(bg_col),  
    fill = bg_col           
  )
)


pdf("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/fig5A.pdf", width = 10, height = 4)
lolliplot(sample.gr, features, 
          ranges = GRanges("PMA1", IRanges(1, 918)),
          xaxis=xaxis,
          yaxis = yaxis,
          legend=legend,
          ylab = "PMA1",
          type = "circle")
dev.off()



###fig5B
PTK2 <- read_excel("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/input/data_fig5AB.xlsx", sheet = 1)
sample.gr <- GRanges("PTK2", 
                     IRanges(
                       start = PTK2$SNP, width = 1, names = PTK2$names))
features <- GRanges("PTK2", IRanges(c(258), 
                                    width=c(304),
                                    names=c("Protein kinase domain"
                                    )))
features$fill <- c("#F7DC6F")
bg_col <- c(
  "BY" = "#F2B56E", 
  "WE" = "#DDA0DD",
  "WA" = "#71BBED"
)
sample.gr$color <- bg_col[PTK2$backgrounds]

features$height <- rep(0.02, 1)
features$featureLayerID <- c("upper_layer")  
xaxis <- c(1, 200, 400, 600, 818)
yaxis <- 6
legend <- list(
  list(
    labels = names(bg_col),  
    fill = bg_col           
  )
)


pdf("E:/Bleomycin_Research_Yeast_Project/06_plot/05_fig5/output/fig5B.pdf", width = 8, height = 4)
lolliplot(sample.gr, features, 
          ranges = GRanges("PTK2", IRanges(1, 818)),
          xaxis=xaxis,
          yaxis = yaxis,
          legend=legend,
          ylab = "PTK2",
          type = "circle")
dev.off()

