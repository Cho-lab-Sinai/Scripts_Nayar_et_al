library(ggplot2)
library(ggpubr)

#convert variable dose from numeric to factor variable
nod2_RISK$NOD2_status <- as.factor(nod2_RISK$NOD2_status)
head(nod2_RISK$NOD2_status)

#basic plot for each gene expression based on nod2 status with stats
NOD2_status <- nod2_RISK$NOD2_status
my_comparisons <- list( c("0", "1"), c("1", "2"), c("0", "2") )

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$PDPN)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 40))
pdpn <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$PDPN, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 40))
p <- pdpn+scale_fill_brewer(palette="Reds")
pdpn <- p + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
pdpn <- pdpn + stat_compare_means(comparisons = my_comparisons)
pdpn 

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CHI3L1)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 500))
chi3l1 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CHI3L1, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 500))
c <- chi3l1+scale_fill_brewer(palette="Reds")
chi3l1 <- c+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
chi3l1 <- chi3l1 + stat_compare_means(comparisons = my_comparisons)
chi3l1

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MPEG1)) + geom_boxplot()
mpeg1 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MPEG1, fill=NOD2_status)) + geom_boxplot()
mp <- mpeg1+scale_fill_brewer(palette="Reds")
mpeg1 <- mp+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mpeg1 <- mpeg1 + stat_compare_means(comparisons = my_comparisons)
mpeg1

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MMP3)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1200))
mmp3 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MMP3, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1200))
mm <- mmp3+scale_fill_brewer(palette="Reds")
mmp3 <- mm+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mmp3 <- mmp3 + stat_compare_means(comparisons = my_comparisons)
mmp3

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MFAP4)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 200))
mfap4 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MFAP4, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 200))
mf <- mfap4+scale_fill_brewer(palette="Reds")
mfap4 <- mf+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mfap4 <- mfap4 + stat_compare_means(comparisons = my_comparisons)
mfap4

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$ITGB3)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 6))
itgb3 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$ITGB3, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 6))
it3 <- itgb3+scale_fill_brewer(palette="Reds")
itgb3 <- it3+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
itgb3 <- itgb3 + stat_compare_means(comparisons = my_comparisons)
itgb3

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$COL1A1)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 400))
col1a1 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$COL1A1, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 400))
col1 <- col1a1+scale_fill_brewer(palette="Reds")
col1a1 <- col1+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
col1a1 <- col1a1 + stat_compare_means(comparisons = my_comparisons)
col1a1

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$IL11)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 50))
il11 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$IL11, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 50))
il <- il11+scale_fill_brewer(palette="Reds")
il11 <- il+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il11 <- il11 + stat_compare_means(comparisons = my_comparisons)
il11

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$ITGAV)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 30))
itgav <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$ITGAV, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 30))
i5 <- itgav+scale_fill_brewer(palette="Reds")
itgav <- i5+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
itgav <- itgav + stat_compare_means(comparisons = my_comparisons)
itgav

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MMP9)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 30))
mmp9 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$MMP9, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 30))
m9 <- mmp9+scale_fill_brewer(palette="Reds")
mmp9 <- m9+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
mmp9 <- mmp9 + stat_compare_means(comparisons = my_comparisons)
mmp9

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$PDGFRA)) + geom_boxplot()
pdgfra <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$PDGFRA, fill=NOD2_status)) + geom_boxplot()
pa <- pdgfra+scale_fill_brewer(palette="Reds")
pdgfra <- pa+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
pdgfra <- pdgfra + stat_compare_means(comparisons = my_comparisons)
pdgfra

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CXCL13)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 300))
cxcl13 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CXCL13, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 300))
cx <- cxcl13+scale_fill_brewer(palette="Reds")
cxcl13 <- cx+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
cxcl13 <- cxcl13 + stat_compare_means(comparisons = my_comparisons)
cxcl13

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CD14)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 200))
cd14 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$CD14, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 200))
c14 <- cd14+scale_fill_brewer(palette="Reds")
cd14 <- c14+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
cd14 <- cd14 + stat_compare_means(comparisons = my_comparisons)
cd14

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$TNF)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 20))
tnf <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$TNF, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 20))
t <- tnf+scale_fill_brewer(palette="Reds")
tnf <- t+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
tnf <- tnf + stat_compare_means(comparisons = my_comparisons)
tnf

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$IL11RA)) + geom_boxplot()
il11ra <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$IL11RA, fill=NOD2_status)) + geom_boxplot()
ilr <- il11ra+scale_fill_brewer(palette="Reds")
il11ra <- ilr+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
il11ra <- il11ra + stat_compare_means(comparisons = my_comparisons)
il11ra

ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$COL5A1)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 40))
col5a1 <- ggplot(nod2_RISK, aes(x=NOD2_status, y=nod2_RISK$COL5A1, fill=NOD2_status)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 40))
c5 <- col5a1+scale_fill_brewer(palette="Reds")
col5a1 <- c5+theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
col5a1 <- col5a1 + stat_compare_means(comparisons = my_comparisons)
col5a1


#make figure
figure <- ggarrange(cd14, il11, tnf, pdpn, chi3l1, mmp3, mmp9, col1a1, col5a1, mfap4, itgav, itgb3, labels= c("CD14", "IL11", "TNF", "PDPN", "CHI3L1", "MMP3", "MMP9", "COL1A1", "COL5A1", "MFAP4", "ITGAV", "ITGB3"), ncol = 4, nrow = 3)
figure
