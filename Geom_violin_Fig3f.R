library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(Seurat)

#create dataframe for genes to plot
vln_dfWT = data.frame( mmp9 = WT[["RNA"]]@data["mmp9",], 
                       col1a1 = WT[["RNA"]]@data["col1a1b",], 
                       col5a1 = WT[["RNA"]]@data["col5a1",], 
                       tnfa = WT[["RNA"]]@data["tnfa",], 
                       wt1a = WT[["RNA"]]@data["wt1a",], 
                       wt1b = WT[["RNA"]]@data["wt1b",], 
                       twist1a = WT[["RNA"]]@data["twist1a",], 
                       pdgfra = WT[["RNA"]]@data["pdgfra",], 
                       tgfb1a = WT[["RNA"]]@data["tgfb1a",], 
                       treatment = WT$treatment, cluster = WT$cluster)

vln_dfnod2 = data.frame( mmp9 = nod2[["RNA"]]@data["mmp9",], 
                       col1a1 = nod2[["RNA"]]@data["col1a1b",], 
                       col5a1 = nod2[["RNA"]]@data["col5a1",], 
                       tnfa = nod2[["RNA"]]@data["tnfa",], 
                       wt1a = nod2[["RNA"]]@data["wt1a",], 
                       wt1b = nod2[["RNA"]]@data["wt1b",], 
                       twist1a = nod2[["RNA"]]@data["twist1a",], 
                       pdgfra = nod2[["RNA"]]@data["pdgfra",], 
                       tgfb1a = nod2[["RNA"]]@data["tgfb1a",], 
                       treatment = nod2$treatment, cluster = nod2$cluster)

#statistical comparison definitions 
my_comparisons <- list(c("Untreated", "DSS only"))

#means for each gene
meanscol1a1WT<- aggregate(col1a1 ~ treatment, vln_dfWT, mean)
meanscol5a1WT<- aggregate(col5a1 ~ treatment, vln_dfWT, mean)
meanswt1aWT<- aggregate(wt1a ~ treatment, vln_dfWT, mean)
meanswt1bWT<- aggregate(wt1b ~ treatment, vln_dfWT, mean)
meansmmp9WT<- aggregate(mmp9 ~ treatment, vln_dfWT, mean)
meanstgfb1aWT<- aggregate(tgfb1a ~ treatment, vln_dfWT, mean)
meanspdgfraWT<- aggregate(pdgfra ~ treatment, vln_dfWT, mean)

meanscol1a1nod2<- aggregate(col1a1 ~ treatment, vln_dfnod2, mean)
meanscol5a1nod2<- aggregate(col5a1 ~ treatment, vln_dfnod2, mean)
meanswt1anod2<- aggregate(wt1a ~ treatment, vln_dfnod2, mean)
meanswt1bnod2<- aggregate(wt1b ~ treatment, vln_dfnod2, mean)
meansmmp9nod2<- aggregate(mmp9 ~ treatment, vln_dfnod2, mean)
meanstgfb1anod2<- aggregate(tgfb1a ~ treatment, vln_dfnod2, mean)
meanspdgfranod2<- aggregate(pdgfra ~ treatment, vln_dfnod2, mean)

#only plot untreated and DSS
col1a1nod2<- ggplot(subset(vln_dfnod2, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = col1a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
col5a1nod2<- ggplot(subset(vln_dfnod2, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = col5a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
wt1nod2<- ggplot(subset(vln_dfnod2, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = wt1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 4))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
mmp9nod2<- ggplot(subset(vln_dfnod2, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = mmp9)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 8))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
tgfbnod2<- ggplot(subset(vln_dfnod2, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = tgfb1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())

col1a1WT<- ggplot(subset(vln_dfWT, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = col1a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
col5a1WT<- ggplot(subset(vln_dfWT, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = col5a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
wt1WT<- ggplot(subset(vln_dfWT, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = wt1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
mmp9WT<- ggplot(subset(vln_dfWT, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = mmp9)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())
tgfbWT<- ggplot(subset(vln_dfWT, treatment %in% c("Untreated", "DSS only")), aes(x = treatment, y = tgfb1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + coord_cartesian(ylim = c(-1, 6))+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6)  + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_blank())

ggarrange(col1a1WT, col5a1WT, wt1WT,mmp9WT, tgfbWT,col1a1nod2,col5a1nod2, wt1nod2, mmp9nod2,tgfbnod2, ncol = 5, nrow = 2)
