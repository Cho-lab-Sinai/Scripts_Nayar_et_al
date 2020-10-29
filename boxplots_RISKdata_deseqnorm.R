library(ggplot2)
library(ggpubr)

#convert variable dose from numeric to factor variable
nod2_RISKdeseqdeseq$topNOD2sum <- as.factor(nod2_RISKdeseqdeseq$topNOD2sum)
head(nod2_RISKdeseqdeseq$topNOD2sum)

#basic plot for each gene expression based on nod2 status with stats
topNOD2sum <- nod2_RISKdeseq$topNOD2sum
my_comparisons <- list( c("0", "1"), c("1", "2"), c("0", "2") )

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$PDPN)) + geom_boxplot(outlier.shape = NA) 
pdpn <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$PDPN, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
p <- pdpn+scale_fill_brewer(palette="Reds")
pdpn <- p + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
pdpn <- pdpn + stat_compare_means(comparisons = my_comparisons) + coord_cartesian(ylim = c(0, 2000))
pdpn 

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CHI3L1)) + geom_boxplot(outlier.shape = NA) 
chi3l1 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CHI3L1, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
c <- chi3l1+scale_fill_brewer(palette="Reds")
chi3l1 <- c+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
chi3l1 <- chi3l1 + stat_compare_means(comparisons = my_comparisons) + coord_cartesian(ylim = c(0, 15000))
chi3l1

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MPEG1)) + geom_boxplot()
mpeg1 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MPEG1, fill=topNOD2sum)) + geom_boxplot()
mp <- mpeg1+scale_fill_brewer(palette="Reds")
mpeg1 <- mp+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mpeg1 <- mpeg1 + stat_compare_means(comparisons = my_comparisons)
mpeg1

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MMP3)) + geom_boxplot(outlier.shape = NA) 
mmp3 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MMP3, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
mm <- mmp3+scale_fill_brewer(palette="Reds")
mmp3 <- mm+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mmp3 <- mmp3 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 5e+04))
mmp3

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MFAP4)) + geom_boxplot(outlier.shape = NA) 
mfap4 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MFAP4, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
mf <- mfap4+scale_fill_brewer(palette="Reds")
mfap4 <- mf+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mfap4 <- mfap4 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 4000))
mfap4

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$ITGB3)) + geom_boxplot(outlier.shape = NA) 
itgb3 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$ITGB3, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
it3 <- itgb3+scale_fill_brewer(palette="Reds")
itgb3 <- it3+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
itgb3 <- itgb3 + stat_compare_means(comparisons = my_comparisons)
itgb3

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$COL1A1)) + geom_boxplot(outlier.shape = NA) 
col1a1 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$COL1A1, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
col1 <- col1a1+scale_fill_brewer(palette="Reds")
col1a1 <- col1+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
col1a1 <- col1a1 + stat_compare_means(comparisons = my_comparisons) + coord_cartesian(ylim = c(0, 50000))
col1a1

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL11)) + geom_boxplot(outlier.shape = NA) 
il11 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL11, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
il <- il11+scale_fill_brewer(palette="Reds")
il11 <- il+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il11 <- il11 + stat_compare_means(comparisons = my_comparisons) + coord_cartesian(ylim = c(0, 2000))
il11

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$ITGAV)) + geom_boxplot(outlier.shape = NA) 
itgav <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$ITGAV, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
i5 <- itgav+scale_fill_brewer(palette="Reds")
itgav <- i5+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
itgav <- itgav + stat_compare_means(comparisons = my_comparisons)
itgav

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MMP9)) + geom_boxplot(outlier.shape = NA) 
mmp9 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$MMP9, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
m9 <- mmp9+scale_fill_brewer(palette="Reds")
mmp9 <- m9+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mmp9 <- mmp9 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 10000))
mmp9

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$PDGFRA)) + geom_boxplot(outlier.shape = NA)
pdgfra <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$PDGFRA, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA)
pa <- pdgfra+scale_fill_brewer(palette="Reds")
pdgfra <- pa+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
pdgfra <- pdgfra + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 6000))
pdgfra

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CXCL13)) + geom_boxplot(outlier.shape = NA) 
cxcl13 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CXCL13, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
cx <- cxcl13+scale_fill_brewer(palette="Reds")
cxcl13 <- cx+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
cxcl13 <- cxcl13 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 6000))
cxcl13

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CD14)) + geom_boxplot(outlier.shape = NA) 
cd14 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CD14, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
c14 <- cd14+scale_fill_brewer(palette="Reds")
cd14 <- c14+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
cd14 <- cd14 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 5000))
cd14

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL11RA)) + geom_boxplot()
il11ra <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL11RA, fill=topNOD2sum)) + geom_boxplot()
ilr <- il11ra+scale_fill_brewer(palette="Reds")
il11ra <- ilr+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il11ra <- il11ra + stat_compare_means(comparisons = my_comparisons)
il11ra

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$COL5A1)) + geom_boxplot(outlier.shape = NA) 
col5a1 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$COL5A1, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
c5 <- col5a1+scale_fill_brewer(palette="Reds")
col5a1 <- c5+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
col5a1 <- col5a1 + stat_compare_means(comparisons = my_comparisons)
col5a1

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6)) + geom_boxplot(outlier.shape = NA) 
il6 <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
il6 <- il6+scale_fill_brewer(palette="Reds")
il6 <- il6+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il6 <- il6 + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 1000))
il6

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6R)) + geom_boxplot(outlier.shape = NA) 
il6r <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6R, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
il6r <- il6r+scale_fill_brewer(palette="Reds")
il6r <- il6r+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il6r <- il6r + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 2000))
il6r

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6ST)) + geom_boxplot(outlier.shape = NA) 
il6st <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$IL6ST, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
il6st <- il6st+scale_fill_brewer(palette="Reds")
il6st <- il6st+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il6st <- il6st + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 4000))
il6st

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$OSMR)) + geom_boxplot(outlier.shape = NA) 
osmr <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$OSMR, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
osmr <- osmr+scale_fill_brewer(palette="Reds")
osmr <- osmr+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
osmr <- osmr + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 1500))
osmr

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$OSM)) + geom_boxplot(outlier.shape = NA) 
osm <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$OSM, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
osm <- osm+scale_fill_brewer(palette="Reds")
osm <- osm+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
osm <- osm + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 1000))
osm

ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CEBPB)) + geom_boxplot(outlier.shape = NA) 
cebpb <- ggplot(nod2_RISKdeseq, aes(x=topNOD2sum, y=nod2_RISKdeseq$CEBPB, fill=topNOD2sum)) + geom_boxplot(outlier.shape = NA) 
cebpb <- cebpb+scale_fill_brewer(palette="Reds")
cebpb <- cebpb+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
cebpb <- cebpb + stat_compare_means(comparisons = my_comparisons)+ coord_cartesian(ylim = c(0, 3000))
cebpb

#make figure
figure <- ggarrange(cd14,pdgfra, pdpn, chi3l1, cxcl13 , mmp3, mmp9, col1a1, il11, il6, osm, il6st, labels= c("CD14", "PDGFRA", "PDPN", "CHI3L1", "CXCL13", "MMP3", "MMP9", "COL1A1", "IL11", "IL6", "OSM", "IL6ST"),font.label = list(size = 8, color = "black"),vjust = 2,hjust = -2,ncol = 4, nrow = 3)
figure
