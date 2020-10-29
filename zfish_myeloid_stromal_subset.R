library(Seurat)
library(dplyr)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

#recluster stromal and myeloid populations
cluster<-c("2", "10", "9", "28", "19", "27", "30", "1", "21", "35", "16", "26")

#subset, run PCA and cluster
mydata<-subset(new.fish, idents = cluster)
DefaultAssay(mydata)<- "RNA"
mydata <- NormalizeData(mydata)
mydata <- FindVariableFeatures(mydata)
mydata <- ScaleData(mydata)
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
ElbowPlot(mydata)
mydata <- FindNeighbors(mydata, dims = 1:15)
mydata <- FindClusters(mydata, resolution = 0.8)
mydata <- RunUMAP(mydata, dims = 1:15)

#visualize
DimPlot(mydata, label = TRUE)

#feature plots split by group
FeaturePlot(DSS, features = c("wt1a", "cebpb", "mmp9", "tnfa"), shape.by = "group", max.cutoff = 3, 
            cols = c("lightgrey", "red"))
markers.to.plot <- c("rna_wt1a", "rna_il11a", "rna_il6st")
DotPlot(scaled, features = c("rna_twist1a","rna_cebpb", "rna_stat3", "rna_col1a1b","rna_col5a1", "rna_col1a2", "rna_tgfb1a","rna_wt1b", "rna_il11a", "rna_mmp9", "rna_nod2"), cols = c("RdYlBu"), dot.scale = 8, 
        group.by = "treatment") + RotatedAxis()

FeaturePlot(DSS, features = c("wt1b", "il11a", "mmp9","tgfb1a","mfap4", "tnfa"), shape.by = "group", 
            +             cols = c("lightgrey", "red"), order = TRUE, pt.size = 1.8, ncol = 3) & NoLegend() + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 0), axis.line = element_line(size = 0.2))

FeaturePlot(stromal_myeloidanno, features = c("il6st", "il6r", "il11ra","osmr","lifrb", "stat3","mapk3", "twist1a", "col1a1b"), 
            +             cols = c("lightgrey", "blue"), order = TRUE, pt.size = 0.2, ncol = 3) & NoLegend() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 0), axis.line = element_line(size = 0.2))
#number of cells per cluster per group
table(Idents(mydata), mydata$orig.ident)

#savefile
saveRDS(stromal_myeloidanno, file = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Myeloid_stromal_recluster/stromal_myeloidanno.rds")

stromal_myeloidanno <- RenameIdents(stromal_myeloid, `0` = "EMT transitioning", `1` = "wt1+ stromal cells", `2` = "Residential fibroblasts", 
                         `3` = "ILC3s", `4` = "Smooth muscle cells", `5` = "Endothelial cells", `6` = "Smooth muscle cells", `7` = "PGE-producing epidermal", `8`= "mfap4+ residential macrophages", `9` = "Mesothelial", `10` = "pdgfra+ fibroblasts", `11` = "Enteric neuron", `12` = "Glial cells",
                         `13` = "Chondrocytes", `14` = "Myoblasts", `15` = "Osteoblasts", `16` = "ILC1?", `17` = "Activated macrophages", `18` = "Metabolic muscle", `19` = "Striated muscle", `20` = "Mesodermal/glial", `21` = "Platelets/osteoblasts")
DimPlot(stromal_myeloidanno, label = TRUE, label.size = 3)

#number of cells expressing a certain gene
pdgfra<-subset(mydata, subset = rna_pdgfra > 0)
table(Idents(pdgfra), pdgfra$orig.ident)
subset(mydata, subset=rna_pdgfra>0 & rna_cebpb>0)


#heatmap
DefaultAssay(stromal_myeloidanno)<- "RNA"
scaled<- ScaleData(stromal_myeloidanno, features = rownames(stromal_myeloidanno))
DefaultAssay(scaled)<- "RNA"
DoHeatmap(scaled, features = c("wt1a", "col1a1b", "col5a1", "tgfb1a","il11a", "mmp9", "tnfa"), size = 1) + scale_fill_viridis_b() + NoLegend()
#redblue heatmap
DoHeatmap(scaled, features = c("wt1a","cebpb", "stat3", "col1a1b", "col5a1", "col1a2", "tgfb1a","il11a", "mmp9", "nod2"), size = 3) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))


#heatmap top 5 markers per cluster
stromye.markers <- FindAllMarkers(scaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- stromye.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(scaled, features = top5$gene, size = 2) + NoLegend() + scale_fill_viridis_b()

#subset into uninflamed, DSS, DSS+BZA
stromal_myeloidanno$annotation<- stromal_myeloidanno@active.ident
Idents(stromal_myeloidanno)<- "treatment"
Untreated<- subset(stromal_myeloidanno, idents = "Untreated")
DSS<- subset(stromal_myeloidanno, idents = "DSS only")
DSSandBZA<- subset(stromal_myeloidanno, idents = "DSS+BZA")
Idents(Untreated)<- "annotation"
Idents(DSS)<- "annotation"
Idents(DSSandBZA)<- "annotation"

#differential expression in one cluster between conditions
Idents(scaled)<- "cluster"
markers<- FindMarkers(scaled, ident.1 = "DSS only", ident.2 = "DSS+BZA", group.by = "treatment", subset.ident = "wt1+ stromal cells")

#differential expression in all clusters between conditions/genotypes--genes of interest (default wilcoxin)
nod2WTDSSnew<- FindMarkers(DSS, ident.1 = "nod2", ident.2 = "WT", group.by = "group", features = c("wt1b","il11b", "mmp9", "mfap4", "tnfa", "col1a1b", "wt1a", "cebpb", "stat3", "nod2"), logfc.threshold = 0.005, min.pct = 0.005)
DSSvsDSSBZA<- FindMarkers(stromal_myeloidanno, ident.1 = "DSS only", ident.2 = "DSS+BZA", features = c("cebpb", "mmp9", "stat3", "col1a1b"), logfc.threshold = 0.01, min.pct = 0.01, group.by = "treatment")

#differential expression to include negbinom
#first check distribution
MetaDSS<- DSS@meta.data
hist(MetaDSS$nCount_RNA)
nod2WTDSSnew<- FindMarkers(DSS, ident.1 = "nod2", ident.2 = "WT", group.by = "group", features = c("wt1b","il11b", "mmp9", "mfap4", "tnfa", "col1a1b", "wt1a", "cebpb", "stat3", "nod2"), logfc.threshold = 0.005, min.pct = 0.005, slot = "counts", test.use = "negbinom")
