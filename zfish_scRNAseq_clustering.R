library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)

#load datasets
S1.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S1raw_feature_bc_matrix/")
S2.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S2raw_feature_bc_matrix/")
S3.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S3raw_feature_bc_matrix/")
S4.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S4raw_feature_bc_matrix/")
S5.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S5raw_feature_bc_matrix/")
S6.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S6raw_feature_bc_matrix/")
S7.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S7raw_feature_bc_matrix/")
S8.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S8raw_feature_bc_matrix/")
S9.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S9raw_feature_bc_matrix/")
S10.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S10raw_feature_bc_matrix/")

#initialize each Seurat object
S1 <- CreateSeuratObject(counts = S1.data, project = "WTctrl1X", min.cells = 3, min.features = 200)
S2 <- CreateSeuratObject(counts = S2.data, project = "nod2ctrl1X", min.cells = 3, min.features = 200)
S3 <- CreateSeuratObject(counts = S3.data, project = "WTDSS1X", min.cells = 3, min.features = 200)
S4 <- CreateSeuratObject(counts = S4.data, project = "nod2DSS1X", min.cells = 3, min.features = 200)
S5 <- CreateSeuratObject(counts = S5.data, project = "WTctrl2X", min.cells = 3, min.features = 200)
S6 <- CreateSeuratObject(counts = S6.data, project = "WTDSS2X", min.cells = 3, min.features = 200)
S7 <- CreateSeuratObject(counts = S7.data, project = "WTDSSBZA2X", min.cells = 3, min.features = 200)
S8 <- CreateSeuratObject(counts = S8.data, project = "nod2ctrl2X", min.cells = 3, min.features = 200)
S9 <- CreateSeuratObject(counts = S9.data, project = "nod2DSS2X", min.cells = 3, min.features = 200)
S10 <- CreateSeuratObject(counts = S10.data, project = "nod2DSSBZA2X", min.cells = 3, min.features = 200)

#set up each object after choosing QC parameters
S1$group <- "WT"
S1 <- subset(S1, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S1 <- NormalizeData(S1, normalization.method = "LogNormalize")
S1 <- FindVariableFeatures(S1, selection.method = "vst", nfeatures = 2000)

S2$group <- "nod2"
S2 <- subset(S2, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S2 <- NormalizeData(S2, normalization.method = "LogNormalize")
S2 <- FindVariableFeatures(S2, selection.method = "vst", nfeatures = 2000)

S3$group <- "WT"
S3 <- subset(S3, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S3 <- NormalizeData(S3, normalization.method = "LogNormalize")
S3 <- FindVariableFeatures(S3, selection.method = "vst", nfeatures = 2000)

S4$group <- "nod2"
S4 <- subset(S4, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S4 <- NormalizeData(S4, normalization.method = "LogNormalize")
S4 <- FindVariableFeatures(S4, selection.method = "vst", nfeatures = 2000)

S5$group <- "WT"
S5 <- subset(S5, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S5 <- NormalizeData(S5, normalization.method = "LogNormalize")
S5 <- FindVariableFeatures(S5, selection.method = "vst", nfeatures = 2000)

S6$group <- "WT"
S6 <- subset(S6, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S6 <- NormalizeData(S6, normalization.method = "LogNormalize")
S6 <- FindVariableFeatures(S6, selection.method = "vst", nfeatures = 2000)

S7$group <- "WT"
S7 <- subset(S7, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S7 <- NormalizeData(S7, normalization.method = "LogNormalize")
S7 <- FindVariableFeatures(S7, selection.method = "vst", nfeatures = 2000)

S8$group <- "nod2"
S8 <- subset(S8, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S8 <- NormalizeData(S8, normalization.method = "LogNormalize")
S8 <- FindVariableFeatures(S8, selection.method = "vst", nfeatures = 2000)

S9$group <- "nod2"
S9 <- subset(S9, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S9 <- NormalizeData(S9, normalization.method = "LogNormalize")
S9 <- FindVariableFeatures(S9, selection.method = "vst", nfeatures = 2000)

S10$group <- "nod2"
S10 <- subset(S10, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S10 <- NormalizeData(S10, normalization.method = "LogNormalize")
S10 <- FindVariableFeatures(S10, selection.method = "vst", nfeatures = 2000)

#perform integration of all samples
immune.anchors <- FindIntegrationAnchors(object.list = list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10), dims = 1:20)
new.fish <- IntegrateData(anchorset = immune.anchors)

#perform integrated analysis
DefaultAssay(new.fish) <- "integrated"
new.fish <- ScaleData(new.fish, verbose = FALSE)
new.fish <- RunPCA(new.fish, verbose = FALSE)
ElbowPlot(new.fish)

#tSNE and clustering
new.fish <- RunUMAP(new.fish, reduction = "pca", dims = 1:18)
new.fish <- FindNeighbors(new.fish, reduction = "pca", dims = 1:18)
new.fish <- FindClusters(new.fish, resolution = 0.8)

#visualization
p1 <- DimPlot(new.fish, reduction = "umap", group.by = "group")
p2 <- DimPlot(new.fish, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(new.fish, reduction = "umap", label = TRUE)


#histograms to plot UMI distribution
Meta<- new.fish@meta.data
hist(Meta$nCount_RNA)

#subset cluster first
sixteencluster<- subset(new.fish, idents = "16")
Metasixteen<- sixteencluster@meta.data
hist(x= Metasixteen$nCount_RNA)

#number of cells per cluster per group
table(Idents(new.fish), new.fish$orig.ident)

#number of cells expressing a certain gene
pdgfra<-subset(new.fish, subset = rna_pdgfra > 0)
table(Idents(pdgfra), pdgfra$orig.ident)
subset(new.fish, subset=rna_pdgfra>0 & rna_cebpb>0)

#rename clusters
new.fish <- RenameIdents(new.fish, `0` = "mix", `1` = "Proliferative myeloid + ILC", `2` = "Smooth muscle", 
                          `3` = "Goblet cells", `4` = "Brush border cells", `5` = "Ionocytes", `6` = "Squamous epithelial", `7` = "mix", `8`= "Enteroendocrine + ILC3", `9` = "pdgfra hi fibro", `10` = "Mesenchymal", `11` = "Corpuscles of stannius", `12` = "Absorptive enterocytes",
                         `13` = "Enteroendocrine", `14` = "Absorptive adipocytes", `15` = "Absorptive adipocytes", `16` = "Stromal-like", `17` = "RBCs", `18` = "HSPCs", `19` = "wt1 hi stromal", `20` = "Epithelial", `21` = "Endothelial", `22` = "LREs", `23` = "Epithelial", `24` = "Pancreatic cells", `25` = "Liver cells", `26` = "Chondrocytes",
                         `27` = "Stromal", `28` = "Activated fibroblasts", `29` = "Respiratory ciliated cells", `30` = "Enteric neuron", `31` = "Neuroendocrine", `32` = "Enteroendocrine", `33` = "Epidermal", `34` = "Epithelial", `35` = "Macrophages", `36` = "Primordial germ cells")
DimPlot(new.fish, label = TRUE, label.size = 3)

#rescale data
DefaultAssay(new.fish)<- "RNA"
scaled<- ScaleData(new.fish, features = rownames(new.fish))
DefaultAssay(scaled)<- "RNA"

#heatmap top 5 markers per cluster
new.fish.markers <- FindAllMarkers(scaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- new.fish.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(scaled, features = top5$gene, size = 2) + NoLegend() + scale_fill_viridis_b()
