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

exnMERGE<-readRDS(file='/~/PFCdata/HMintegrate_ExN_subcluster.rds')
innMERGE<-readRDS(file='/~/PFCdata/HMintegrate_InN_subcluster.rds')

HMNmerge<-merge(exnMERGE,innMERGE)

dend1=readRDS(file="/~/PFCdata/HMNeuron_dendrogram_0705.rds")
NeuronOrder=c(labels(dend1))

allplotgene=c("SLC17A7","GAD1","CUX2","PDGFD","COL5A2","RORB","ETV1","SULF2","OPRK1","TLE4","PVALB","SST","VIP","ID2")
HMNmerge@meta.data$subcluster=factor(HMNmerge@meta.data$subcluster,levels= NeuronOrder)

Idents(HMNmerge)<-'subcluster'
Avgexp=AverageExpression(HMNmerge,assays = "RNA",features = c(allplotgene),use.scale = F)
dfmatrix=as.matrix(Avgexp$RNA)
df<-data.frame(Avgexp$RNA/rowMaxs(dfmatrix))
df['genename']<-rownames(df)

allcolor=read.csv(file="/~/PFCdata/HMNeuronclusterscolor_0705.csv")
allcolor
cluster_color1=(allcolor[c(1:94),])$clusterColor
cluster_color1
names(cluster_color1)=(allcolor[c(1:94),])$clusterName

my_sample_col=clusterdf
ann_colors=ann_colors
df$genename<-NULL
ae=as.matrix(df) 

colnames(ae)<-levels(HMNmerge@meta.data$subcluster)
celltypeColor1<-c(rep("#7294D4",43),rep("#0B775E",51))
names(celltypeColor1)<-c(colnames(ae))
p <-Heatmap(ae,
            col =c(colorRampPalette(colors = c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026"))),
            show_row_names = T,cluster_rows =F, cluster_columns = F,show_column_names = T,
            top_annotation = HeatmapAnnotation(Cluster= colnames(ae),Celltype= colnames(ae),
                                               col = list(Cluster = cluster_color1,Celltype=celltypeColor1),show_annotation_name = T), heatmap_width = unit(35, "cm"), heatmap_height = unit(15, "cm"))
p

