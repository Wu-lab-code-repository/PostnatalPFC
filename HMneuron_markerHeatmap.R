library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(batchelor)
library(Matrix)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)

exnMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_ExN_subcluster.rds')
innMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_InN_subcluster.rds')

DefaultAssay(exnMERGE)<-'RNA'
DefaultAssay(innMERGE)<-'RNA'
HMNmerge<-merge(exnMERGE,innMERGE)

DefaultAssay(HMNcca)<-'RNA'
HMNcca <- NormalizeData(HMNcca, normalization.method = "LogNormalize", scale.factor = 10000)

dend1=readRDS(file="/home/gpfs/home/wulab15/species_merge/HMNeuron_dendrogram_0705.rds")
NeuronOrder=c(labels(dend1))
NeuronOrder


allplotgene=c("SLC17A7","GAD1","CUX2","PDGFD","COL5A2","RORB","ETV1","SULF2","OPRK1","TLE4","PVALB","SST","VIP","ID2")
HMNmerge@meta.data$subcluster=factor(HMNmerge@meta.data$subcluster,levels= NeuronOrder)

Idents(HMNmerge)<-'subcluster'
Avgexp=AverageExpression(HMNmerge,assays = "RNA",features = c(allplotgene),use.scale = F)
dfmatrix=as.matrix(Avgexp$RNA)
df<-data.frame(Avgexp$RNA/rowMaxs(dfmatrix))
df['genename']<-rownames(df)
df
posorder=c()
for(i in 1:length(allplotgene)){
  posi=which(rownames(test)%in%allplotgene[i])
  posorder=c(posorder,posi)
}
newdf=rbind(df[c(posorder),])

allcolor=read.csv(file="/home/gpfs/home/wulab15/species_merge/HMNeuronclusterscolor_0705.csv")
allcolor
cluster_color1=(allcolor[c(1:94),])$clusterColor
cluster_color1
names(cluster_color1)=(allcolor[c(1:94),])$clusterName

my_sample_col=clusterdf
ann_colors=ann_colors
newdf$genename<-NULL
ae=as.matrix(newdf) 

colnames(ae)<-levels(HMNmerge@meta.data$subcluster)
celltypeColor1<-c(rep("#7294D4",43),rep("#0B775E",51))
names(celltypeColor1)<-c(colnames(ae))
B <-Heatmap(ae,
            col =c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026"))),
            show_row_names = T,cluster_rows =F, cluster_columns = F,show_column_names = T,
            top_annotation = HeatmapAnnotation(Cluster= colnames(ae),Celltype= colnames(ae),
                                               col = list(Cluster = cluster_color1,Celltype=celltypeColor1),show_annotation_name = T), heatmap_width = unit(35, "cm"), heatmap_height = unit(15, "cm"))
B

pdf("/home/gpfs/home/wulab15/Figures/HMneuron_markerheatmap.pdf",height = 18,width = 22)   
B
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/HMneuron_markerheatmap.tiff",height = 1200,width = 1500)
B
dev.off()