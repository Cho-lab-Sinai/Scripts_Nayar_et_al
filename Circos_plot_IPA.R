
#######################
library(dplyr)
library(BioCircos)
tracks=NULL

##Load file with all the columns needed, first columns as chromosome names (or genes names)
n1 <- read.csv("~/Desktop/Graduate School/Cho lab/R files/circosileal_infmac.csv")
n1<-n1  %>% arrange(desc(lengtht,lengthL)) %>% group_by(avg_logFC)
n1$chromosome<-as.character(n1$chromosome)

#make circumferences
genomeChr = as.list(n1$chromosome)
lengthChr =n1$length
names(lengthChr) <- genomeChr

#WT1
genomeChr1 = as.list(n1$chromosome)
lengthChr1 = n1$length2T
names(lengthChr1) <- genomeChr1


#WT1
#STAT3
genomeChr2 = as.list(n1$chromosome)
lengthChr2 = n1$length2L
names(lengthChr2) <- genomeChr2



# Define boxes positions

boxChromosomes = rep(genomeChr,lengthChr)
boxChromosomes1 = rep(genomeChr1,lengthChr1)
boxChromosomes2 = rep(genomeChr2,lengthChr2)


# Define values for two heatmap tracks
boxVal1 = rep(n1$avg_logFC,lengthChr)
#WT1
val2 = rep(n1$lengtht,lengthChr1)
#STAT3
val3 = rep(n1$lengthL,lengthChr2)



#define tracks

tracks = BioCircosHeatmapTrack("heatmap1", boxChromosomes,1,500,
                               boxVal1, minRadius = 0.9, maxRadius = 1.15,color = c("red", "green"),range=c(-1,2.7))
tracks = tracks+BioCircosHeatmapTrack("heatmap2", boxChromosomes1,1,500,values = val2,range=c(0,1),
                                      minRadius = 0.85, maxRadius = 0.7,color = c("white","#14C6CC")) #14c6cc
tracks = tracks+BioCircosHeatmapTrack("heatmap3", boxChromosomes2,1,500,values = val3,range=c(0,1),
                                      minRadius = 0.65, maxRadius = 0.5,color = c("white","#D8BFD8"))

#Put them all together
BioCircos(tracks, genome = as.list(lengthChr),keep_vec_names=TRUE,chrPad = 0.02,displayGenomeBorder = TRUE,genomeBorderSize = 0.5, genomeTicksLen = 0, genomeTicksTextSize = 0, genomeTicksScale = 1e+8,genomeLabelTextSize = "11pt", genomeLabelDy = 50,genomeLabelOrientation = 90)
