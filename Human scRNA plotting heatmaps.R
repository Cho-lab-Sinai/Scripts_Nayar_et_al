#Effi rds file data extraction
#################################################################
for(i in 1:20){
  colnames(f)[i]=annot[i]
}
#inflamed change path
ldm = load_scDissector_data(clustering_data_path="/.......SCD", model_name="combined_081518", sample_names=c(69, 122, 128, 138, 158, 181, 187, 190, 193, 196, 209))

#uninflamed change path
ldm = load_scDissector_data(clustering_data_path="....../SCD", model_name="combined_081518", sample_names=c(68, 123, 129, 135, 159, 180, 186, 189, 192, 195, 208))


inflamedMat<-ldm$model$models
annot<-ldm$clustAnnots
mat<-as.data.frame(inflamedMat)
mat$gene<-rownames(mat)
f<-filter(mat,gene%in%list)
rownames(f)<-f$gene
f<-f[,-51]
f2<-f2[,c("ResidentMacrophages","Inf. Macrophages","Pericytes","Smooth muscle cells" , "Fibroblasts","Activated fibroblasts")]



#blood
ldm = load_scDissector_data(clustering_data_path="/Users/mamtagiri/Downloads/SCD", model_name="pbmc_111718_20", sample_names=c(67,114,119,126,127,134,157,179,185,188,191,194))
inflamedMat<-ldm$model$models
annot<-ldm$clustAnnots
mat<-as.data.frame(inflamedMat)
mat$gene<-rownames(mat)
f<-filter(mat,gene%in%list)
rownames(f)<-f$gene
f<-f[,-21]
f$ClassicalMonocytes = (f[,2] + f[,19]) / 2
f$ResidentMacrophages = (f[,15] + f[,41]) / 2
list<-c("CD14","NOD2","MPEG1","TNF","MMP9","CCL2","IL6ST","CEBPB","ITGAV","ITGB3","COL5A1","PDGFRA","MFAP4","COL1A1","PDPN","CHI3L1","MMP3","IL11","CXCL13","WT1")
f2<-f2[,c("ClassicalMonocytes","Intermediate monocytes","Non classical monocytes","DC (cDC+pDC)", "Platelets")]
f2<-f2[c("CD14","NOD2","MPEG1","TNF","MMP9","CCL2","IL6ST","CEBPB","ITGAV","ITGB3","COL5A1","PDGFRA","MFAP4","COL1A1","PDPN","CHI3L1","MMP3","IL11","CXCL13","WT1"),]


#########################################
#heatmap
#########################################
colgrad_rel_file=paste("colors_brewer_RdBu.txt",sep="")
colgrad_rel<<-read.table(colgrad_rel_file,stringsAsFactors=F)[,1]

datafile<-read.csv("RelativeInflamed.csv",row.names = 1)

heatmaply(t(datafile),Rowv=F,Colv=F, scale = "none",limits=c(-2,4),main_title = "",colors  =colgrad_rel,cex.genes = .6,cex.clusters = .5,line.genes = .1)

