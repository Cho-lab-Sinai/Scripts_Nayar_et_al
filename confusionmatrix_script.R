library(ggplot2)
library(ggpubr)
library(heatmaply)
library(viridis)
library(pheatmap)

#humantohuman
confusionhuman<- read.csv("~/Desktop/Extrafish_scRNAseq/Human_fish/Final/confusionhuman.csv", row.names = 1)
datahuman<- as.matrix(confusionhuman)
pheatmap(datahuman, display_numbers = TRUE, color = viridis(60), cluster_rows = FALSE, cluster_cols = FALSE, fontsize_number = 15, fontsize_row = 15, fontsize_col = 15,angle_col = 0) 

#humantofish
confusionprop<- read.csv("~/Desktop/Extrafish_scRNAseq/Human_fish/Final/confusionfishhuman_prop.csv", row.names = 1)
dataprop<- as.matrix(confusionprop)
heatmap.2(dataprop, dendogram = c("none"), Rowv = FALSE, Colv = FALSE, trace = "none", col = "viridis")
pheatmap(dataprop, display_numbers = TRUE, color = viridis(60), cluster_rows = FALSE, cluster_cols = FALSE, fontsize_number = 15, fontsize_row = 15, fontsize_col = 15,angle_col = 0) 