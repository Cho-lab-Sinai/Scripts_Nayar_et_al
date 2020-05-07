library(Seurat)
library(dplyr)
cluster<-c("Col1a1a fibroblasts","Col5a1 fibroblasts","Activated fibroblasts","Mesothelial","Chondrocytes","Smooth muscle")
mydata<-subset(fish.new, idents = cluster)
mydata <- NormalizeData(mydata)
mydata <- FindVariableFeatures(mydata)
mydata <- ScaleData(mydata)
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
ElbowPlot(mydata)
mydata <- FindNeighbors(mydata, dims = 1:7)
mydata <- FindClusters(mydata, resolution = 0.5)
mydata <- RunUMAP(mydata, dims.use = 1:7)

#for overall dimplot
DimPlot(mydata, reduction = "umap")
p1 <- DimPlot(object = mydata, reduction = "umap", group.by = "group")
p2 <- DimPlot(object = mydata, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

#genebycluster
FeaturePlot(mydata, features = c("col5a1", "col1a1a","postnb","mfap4","pcna","il11ra","kdrl"), min.cutoff = "q9")

#find markers
four.markers <- FindConservedMarkers(mydata, ident.1 = 4, grouping.var = "group", 
                                        verbose = FALSE)
head(four.markers)

saveRDS(stromal.recluster, file = "/Users/shikhanayar/Desktop/Graduate School/Cho lab/R files/Fishnew_singlecell/stromal.recluster.rds")

#rename clusters
stromal.recluster <- RenameIdents(mydata, `0` = "pdgfra+ fibroblasts", `1` = "Mesenchymal/lymphatic", `2` = "Chondrocytes", 
                         `3` = "Smooth muscle", `4` = "Smooth muscle", `5` = "col1a1a+ fibroblasts", `6` = "Mesothelial", `7` = "Mesenchymal/lymphatic", `8`= "Mesenchymal/lymphatic", `9` = "Chondrocytes")
DimPlot(stromal.recluster, label = TRUE)

#number of cells per cluster per group
table(Idents(stromal.recluster), stromal.recluster$orig.ident)

