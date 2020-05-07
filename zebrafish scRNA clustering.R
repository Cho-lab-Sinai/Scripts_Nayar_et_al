library(Seurat)
library(dplyr)
library(cowplot)
WT1.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNJCZ1/raw")
WT2.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNZF1/raw")
WT3.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNJC02/raw")
DSS1.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNJCZ2/raw")
DSS2.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNZF2/raw")
DSS3.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNJC03/raw")
Ecoli.data <- Read10X(data.dir = "/Users/shikhanayar/Google Drive/fish/SNZF3/raw")

#set up object 1
ctrl1 <- CreateSeuratObject(counts = WT1.data, project = "ctrl1", min.cells = 5)
ctrl1$group <- "CTRL"
ctrl1 <- subset(ctrl1, subset = nFeature_RNA > 150)
ctrl1 <- NormalizeData(ctrl1, verbose = FALSE)
ctrl1 <- FindVariableFeatures(ctrl1, selection.method = "vst", nfeatures = 2000)

#set up object 2
ctrl2 <- CreateSeuratObject(counts = WT2.data, project = "ctrl2", min.cells = 5)
ctrl2$group <- "CTRL"
ctrl2 <- subset(ctrl2, subset = nFeature_RNA > 150)
ctrl2 <- NormalizeData(ctrl2, verbose = FALSE)
ctrl2 <- FindVariableFeatures(ctrl2, selection.method = "vst", nfeatures = 2000)

#set up object 3
ctrl3 <- CreateSeuratObject(counts = WT3.data, project = "ctrl3", min.cells = 5)
ctrl3$group <- "CTRL"
ctrl3 <- subset(ctrl3, subset = nFeature_RNA > 150)
ctrl3 <- NormalizeData(ctrl3, verbose = FALSE)
ctrl3 <- FindVariableFeatures(ctrl3, selection.method = "vst", nfeatures = 2000)

#set up object 4
stim1 <- CreateSeuratObject(counts = DSS1.data, project = "stim1", min.cells = 5)
stim1$group <- "DSS"
stim1 <- subset(stim1, subset = nFeature_RNA > 150)
stim1 <- NormalizeData(stim1, verbose = FALSE)
stim1 <- FindVariableFeatures(stim1, selection.method = "vst", nfeatures = 2000)

#set up object 5
stim2 <- CreateSeuratObject(counts = DSS2.data, project = "stim2", min.cells = 5)
stim2$group <- "DSS"
stim2 <- subset(stim2, subset = nFeature_RNA > 150)
stim2 <- NormalizeData(stim2, verbose = FALSE)
stim2 <- FindVariableFeatures(stim2, selection.method = "vst", nfeatures = 2000)

#set up object 6
stim3 <- CreateSeuratObject(counts = DSS3.data, project = "stim3", min.cells = 5)
stim3$group <- "DSS"
stim3 <- subset(stim3, subset = nFeature_RNA > 150)
stim3 <- NormalizeData(stim3, verbose = FALSE)
stim3 <- FindVariableFeatures(stim3, selection.method = "vst", nfeatures = 2000)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl1, ctrl2, ctrl3, stim1, stim2, stim3), dims = 1:15)
fish.all <- IntegrateData(anchorset = immune.anchors, dims = 1:15)

#perform integrated analysis
DefaultAssay(fish.all) <- "integrated"
fish.new <- ScaleData(fish.all, verbose = FALSE)
fish.new <- RunPCA(fish.new, npcs = 30, verbose = FALSE)

#tSNE and clustering
fish.new <- RunUMAP(fish.new, reduction = "pca", dims = 1:15)
fish.new <- FindNeighbors(fish.new, reduction = "pca", dims = 1:15)
fish.new <- FindClusters(fish.new, resolution = 0.5)

#visualization
p1 <- DimPlot(fish.new, reduction = "umap", group.by = "group")
p2 <- DimPlot(fish.new, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(fish.new, reduction = "umap", split.by = "group")

#identify conserved cell type markers
DefaultAssay(fish.new) <- "RNA"
mesothelial.markers <- FindConservedMarkers(fish.all, ident.1 = "Mesothelial", grouping.var = "group", verbose = FALSE)
head(ten.markers)

#marker expression
FeaturePlot(fish.new, features = c("postnb","pcna", "il11ra","pcolce2","ackr4", "mfap4", "itgav" 
                                          ), min.cutoff = "q9")
FeaturePlot(fish.new, features = c("col1a2", "il11ra", "mfap4", "pcna"), split.by = "group", max.cutoff = 3, 
            cols = c("lightgrey", "red"))

#rename clusters
fish.new <- RenameIdents(fish.new, `0` = "Exocrine", `1` = "Enterocytes", `2` = "Parietal cells", 
                              `3` = "Goblet cells", `4` = "Activated fibroblasts", `5` = "Smooth muscle", `6` = "Macrophages", `7` = "Col5a1 fibroblasts", `8` = "Chondrocytes", `9` = "Hepatocytes", 
                              `10` = "Enteroendocrine", `11` = "LREs", `12` = "HSPCs", `13` = "Adipocytes", `14` = "Epithelial", `15` = "Col1a1a fibroblasts", `16` = "Pancreatic", `17` = "Epithelial", `18` = "Enterocytes", `19` = "Endocrine", `20` = "Mesothelial", `21` = "RBCs")
DimPlot(fish.new, label = TRUE)

#savefile
saveRDS(fish.new, file = "/Users/shikhanayar/Desktop/Graduate School/Cho lab/R files/fish.new.rds")



#dotplot with splitby
Idents(fish.new) <- factor(Idents(fish.new), levels = c("Intestinal exocrine", "RBCs", 
                                                                      "Hepatocytes", "Chondrocytes", "Epithelial", "Neuroendocrine", "Macrophages?", "Muscle", "Enteroendocrine", "Stromal", "Col1a1a fibro", "Keratinocytes", 
                                                                      "Fibroblasts", "Activated fibro", "Endothelial", "Epithelial", "HSPCs"))
markers.to.plot <- c("il6st", "cebpb", "pdgfra")
DotPlot(fish.new, features = rev(markers.to.plot), cols = c("blue", "red", "green", "yellow"), dot.scale = 8, 
        split.by = "group") + RotatedAxis()

#violin plot
plots <- VlnPlot(fish.new, features = c("il6st", "cebpb", "pdgfra"), split.by = "group", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

#heatmap
top10 <- stromal_recluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(stromal_recluster, features = c("ccl2", "lmo2", "ccl25b","cfd","rbpms2b","cxcl12a","cald1b","acta2", "tagln","acta1b","myl9a","desmb","fhl1a","rna_pdgfra","ptgs2a","ptgs2b","col4a6","col4a5","cxcl18b", "ifitm1","pfn1","slc16a3","angptl5","alcama","bricd5","matn1","gyg1b","col11a1a","krt8","loxl1","col1a1a","tat","msx2b","twist1a", "col1a2", "postnb","ptgdsb.1","prrx1b","mdh1aa","vcanb"))

#differential expression between clusters
col1a1a.cells <- subset(stromal_recluster, idents = "col1a1a+ fibroblasts")
Idents(col1a1a.cells) <- "group"
avg.col1a1a.cells <- log1p(AverageExpression(col1a1a.cells, verbose = FALSE)$RNA)
avg.col1a1a.cells$gene <- rownames(avg.col1a1a.cells)

meso.cells <- subset(stromal_recluster, idents = "Mesothelial cells")
Idents(meso.cells) <- "group"
avg.meso.cells <- log1p(AverageExpression(meso.cells, verbose = FALSE)$RNA)
avg.meso.cells$gene <- rownames(avg.meso.cells)

genes.to.label = c("col1a1a", "col1a2", "cebpb", "mmp9", "pcna", "il11a", "wt1b", "nod2", "il11ra", "pdgfra")
p1 <- ggplot(avg.col1a1a.cells, aes(CTRL, DSS)) + geom_point() + ggtitle("col1a1a+ activated fibroblasts")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.meso.cells, aes(CTRL, DSS)) + geom_point() + ggtitle("Mesothelial cells")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

