library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(Seurat)

#create dataframe for genes to plot
vln_dfall = data.frame(stat3 = stromal_myeloidanno[["RNA"]]@data["stat3",], 
                       mmp9 = stromal_myeloidanno[["RNA"]]@data["mmp9",], 
                       cebpb = stromal_myeloidanno[["RNA"]]@data["cebpb",], 
                       col1a1 = stromal_myeloidanno[["RNA"]]@data["col1a1b",], 
                       col5a1 = stromal_myeloidanno[["RNA"]]@data["col5a1",], 
                       tgfb1a = stromal_myeloidanno[["RNA"]]@data["tgfb1a",], 
                       tnfa = stromal_myeloidanno[["RNA"]]@data["tnfa",], 
                       wt1a = stromal_myeloidanno[["RNA"]]@data["wt1a",], 
                       wt1b = stromal_myeloidanno[["RNA"]]@data["wt1b",], 
                       il11a = stromal_myeloidanno[["RNA"]]@data["il11a",], 
                       il6st = stromal_myeloidanno[["RNA"]]@data["il6st",], 
                       pcna = stromal_myeloidanno[["RNA"]]@data["pcna",], 
                       pdgfra = stromal_myeloidanno[["RNA"]]@data["pdgfra",], 
                       tgfb1a = stromal_myeloidanno[["RNA"]]@data["tgfb1a",], 
                       treatment = stromal_myeloidanno$treatment, genotype = stromal_myeloidanno$group, cluster = stromal_myeloidanno$cluster)

#statistical comparison definitions 
my_comparisons <- list( c("Untreated", "DSS only"), c("Untreated", "DSS+BZA"), c("DSS only", "DSS+BZA"))


#means for each gene
meansstat3<- aggregate(stat3 ~ treatment, vln_dfall, mean)
meanscebpb<- aggregate(cebpb ~ treatment, vln_dfall, mean)
meanscol1a1<- aggregate(col1a1 ~ treatment, vln_dfall, mean)
meanscol5a1<- aggregate(col5a1 ~ treatment, vln_dfall, mean)
meanswt1a<- aggregate(wt1a ~ treatment, vln_dfall, mean)
meanswt1b<- aggregate(wt1b ~ treatment, vln_dfall, mean)
meansil11a<- aggregate(il11a ~ treatment, vln_dfall, mean)
meansil6st<- aggregate(il6st ~ treatment, vln_dfall, mean)
meansmmp9<- aggregate(mmp9 ~ treatment, vln_dfall, mean)
meanstgfb1a<- aggregate(tgfb1a ~ treatment, vln_dfall, mean)
meanspdgfra<- aggregate(pdgfra ~ treatment, vln_dfall, mean)
meanspcna<- aggregate(pcna ~ treatment, vln_dfall, mean)


#individual gene plots
stat3plot<- ggplot(vln_dfall, aes(x = treatment, y = stat3)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
cebpbplot<- ggplot(vln_dfall, aes(x = treatment, y = cebpb)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 8)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
col1a1plot<- ggplot(vln_dfall, aes(x = treatment, y = col1a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 8)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
col5a1plot<- ggplot(vln_dfall, aes(x = treatment, y = col5a1)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 8)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
wt1aplot<- ggplot(vln_dfall, aes(x = treatment, y = wt1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
wt1bplot<- ggplot(vln_dfall, aes(x = treatment, y = wt1b)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
il11aplot<- ggplot(vln_dfall, aes(x = treatment, y = il11a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
il6stplot<- ggplot(vln_dfall, aes(x = treatment, y = il6st)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
mmp9plot<- ggplot(vln_dfall, aes(x = treatment, y = mmp9)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 8)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
tgfb1aplot<- ggplot(vln_dfall, aes(x = treatment, y = tgfb1a)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 5)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
pdgfraplot<- ggplot(vln_dfall, aes(x = treatment, y = pdgfra)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
pcnaplot<- ggplot(vln_dfall, aes(x = treatment, y = pcna)) + geom_violin(aes(fill = treatment), trim=FALSE, scale = "count", adjust = 3) + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 6) + coord_cartesian(ylim = c(-1, 6)) + scale_fill_viridis_d() + theme_bw() + geom_point() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))

#make figure
ggarrange(stat3plot, cebpbplot, mmp9plot, col1a1plot, col5a1plot, il11aplot, wt1bplot, tgfb1aplot, il6stplot, pdgfraplot, pcnaplot, ncol = 3, nrow = 4)


