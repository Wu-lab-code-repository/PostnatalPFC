library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)
library(patchwork)
library(cowplot)
library(sampling)


################HM ExN
MacacaALL<-readRDS(file="/home/gpfs/home/wulab15/MacacaPFC/MacacaPFC_withsubcluster_0103.rds")
HumanALL<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/HumanPFC_withsubcluster_0103.rds")
length(MacacaALL@meta.data$celltype)
length(HumanALL@meta.data$celltype)
Idents(HumanALL)<-'celltype'
Idents(MacacaALL)<-'celltype'
HExN<-subset(HumanALL,idents=c("Excitatory Neuron"))
MExN<-subset(MacacaALL,idents=c("Excitatory Neuron"))
MExN@meta.data['species']<-'Macaca'
HExN@meta.data['species']<-'Human'

MExN@meta.data$sample_age=paste0('M_',MExN@meta.data$sample_age)
MExN@meta.data$sample_Macaca=paste0('M_',MExN@meta.data$sample_Macaca)

HExN@meta.data$sample_age=paste0('H_',HExN@meta.data$sample_age)
HExN@meta.data$sample_Human=paste0('H_',HExN@meta.data$sample_Human)

HExN@meta.data['sample']=HExN@meta.data$sample_Human
HExN@meta.data$sample_Human=NULL
MExN@meta.data['sample']=MExN@meta.data$sample_Macaca
MExN@meta.data$sample_Macaca=NULL

MExN@meta.data$time=NULL

#downsampling MExN data
A=50000/length(MExN@meta.data$orig.ident) 
cluster_inf <- MExN@meta.data[,c("sample","subcluster")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subcluster),]
k2 <- as.numeric(round(A*table(cluster_inf$subcluster)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subcluster'), size=k2, description=FALSE)
table(sampling_cluster.inf$subcluster)
pos <- rownames(sampling_cluster.inf)
pos <- as.integer(pos)
sampling_cluster.inf <- cluster_inf[pos,]
write.table(sampling_cluster.inf, file = "/home/gpfs/home/wulab15/MacacaPFC/species_merge/MExN_sampling_cluster_inf_0209.csv",sep = ",")
sampling_cluster.inf <- read.csv('/home/gpfs/home/wulab15/MacacaPFC/species_merge/MExN_sampling_cluster_inf_0209.csv')

table(sampling_cluster.inf$subcluster)
(1/A)*table(sampling_cluster.inf$sample) 
table(MExN@meta.data$sample)
subset_cellnames <- rownames(sampling_cluster.inf)
MExN <- SubsetData(MExN,cells = subset_cellnames)
table(MExN@meta.data$sample)
length(MExN@meta.data$sample)

#downsampling HExN data
B=50000/length(HExN@meta.data$orig.ident) 
cluster_inf <- HExN@meta.data[,c("sample","subcluster")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subcluster),]
k2 <- as.numeric(round(B*table(cluster_inf$subcluster)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subcluster'), size=k2, description=FALSE)
table(sampling_cluster.inf$subcluster)
pos <- rownames(sampling_cluster.inf)
pos <- as.integer(pos)
sampling_cluster.inf <- cluster_inf[pos,]
write.table(sampling_cluster.inf, file = "/home/gpfs/home/wulab15/MacacaPFC/species_merge/HExN_sampling_cluster_inf_0209.csv",sep = ",")
sampling_cluster.inf <- read.csv('/home/gpfs/home/wulab15/MacacaPFC/species_merge/HExN_sampling_cluster_inf_0209.csv')

table(sampling_cluster.inf$subcluster)
(1/B)*table(sampling_cluster.inf$sample) 
table(HExN@meta.data$sample)
subset_cellnames <- rownames(sampling_cluster.inf)
HExN <- SubsetData(HExN,cells = subset_cellnames)
table(HExN@meta.data$sample)
length(HExN@meta.data$sample)

objMERGE <- merge(x=HExN,y=MExN)
objMERGE@meta.data$seurat_clusters=NULL

head(objMERGE@meta.data)

DefaultAssay(objMERGE)<-'RNA'
splitlist <- SplitObject(objMERGE,split.by='species')
for (i in 1:length(splitlist)) {
  splitlist[[i]] <- NormalizeData(splitlist[[i]])
  splitlist[[i]] <- FindVariableFeatures(splitlist[[i]],selection.method = 'vst',nfeatures = 3000)
}

splitanchors <- FindIntegrationAnchors(object.list = splitlist,anchor.features = 3000)
newMERGE <- IntegrateData(anchorset = splitanchors)

DefaultAssay(newMERGE) <- 'integrated'
newMERGE <- ScaleData(newMERGE,vars.to.regress = c("percent.mt",'nCount_RNA')) 
newMERGE <- RunPCA(newMERGE, npcs = 100)
ElbowPlot(newMERGE,ndims = 100, reduction = "pca")

batch_info <- as.factor(newMERGE@meta.data$sample)
table(batch_info)
pca_coord <- newMERGE@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
newMERGE@reductions$pca@cell.embeddings <- pca_corrected$corrected

newMERGE <- FindNeighbors(newMERGE, dims = 1:25)
newMERGE <- FindClusters(newMERGE, resolution = 2)
newMERGE <- RunUMAP(newMERGE,dims = 1:25, reduction='pca')
DimPlot(newMERGE,label=TRUE)

Idents(newMERGE)<-'species'
DimPlot(newMERGE,cols=c("#3D578E","plum3"))+facet_wrap(~ident)
DimPlot(newMERGE,cols=c("#3D578E","plum3"))
saveRDS(newMERGE,file='/home/gpfs/home/wulab15/species_merge/HMintegrate_ExN_subcluster.rds')

################HM InN

MacacaALL<-readRDS(file="/home/gpfs/home/wulab15/MacacaPFC/MacacaPFC_withsubcluster_0103.rds")
HumanALL<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/HumanPFC_withsubcluster_0103.rds")
length(MacacaALL@meta.data$celltype)
length(HumanALL@meta.data$celltype)
Idents(HumanALL)<-'celltype'
Idents(MacacaALL)<-'celltype'
HInN<-subset(HumanALL,idents=c("Inhibitory Neuron"))
MInN<-subset(MacacaALL,idents=c("Inhibitory Neuron"))
MInN@meta.data['species']<-'Macaca'
HInN@meta.data['species']<-'Human'

MInN@meta.data$sample_age=paste0('M_',MInN@meta.data$sample_age)
MInN@meta.data$sample_Macaca=paste0('M_',MInN@meta.data$sample_Macaca)

HInN@meta.data$sample_age=paste0('H_',HInN@meta.data$sample_age)
HInN@meta.data$sample_Human=paste0('H_',HInN@meta.data$sample_Human)

HInN@meta.data['sample']=HInN@meta.data$sample_Human
HInN@meta.data$sample_Human=NULL
MInN@meta.data['sample']=MInN@meta.data$sample_Macaca
MInN@meta.data$sample_Macaca=NULL

MInN@meta.data$time=NULL

#downsampling MInN data
A=20000/length(MInN@meta.data$orig.ident) 
cluster_inf <- MInN@meta.data[,c("sample","subcluster")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subcluster),]
k2 <- as.numeric(round(A*table(cluster_inf$subcluster)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subcluster'), size=k2, description=FALSE)
table(sampling_cluster.inf$subcluster)
pos <- rownames(sampling_cluster.inf)
pos <- as.integer(pos)
sampling_cluster.inf <- cluster_inf[pos,]
write.table(sampling_cluster.inf, file = "/home/gpfs/home/wulab15/MacacaPFC/species_merge/MInN_sampling_cluster_inf_0209.csv",sep = ",")
sampling_cluster.inf <- read.csv('/home/gpfs/home/wulab15/MacacaPFC/species_merge/MInN_sampling_cluster_inf_0209.csv')

table(sampling_cluster.inf$subcluster)
(1/A)*table(sampling_cluster.inf$sample) 
table(MInN@meta.data$sample)
subset_cellnames <- rownames(sampling_cluster.inf)
MInN <- SubsetData(MInN,cells = subset_cellnames)
table(MInN@meta.data$sample)
length(MInN@meta.data$sample)

#downsampling HInN data
B=20000/length(HInN@meta.data$orig.ident) 
cluster_inf <- HInN@meta.data[,c("sample","subcluster")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subcluster),]
k2 <- as.numeric(round(B*table(cluster_inf$subcluster)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subcluster'), size=k2, description=FALSE)
table(sampling_cluster.inf$subcluster)
pos <- rownames(sampling_cluster.inf)
pos <- as.integer(pos)
sampling_cluster.inf <- cluster_inf[pos,]
write.table(sampling_cluster.inf, file = "/home/gpfs/home/wulab15/MacacaPFC/species_merge/HInN_sampling_cluster_inf_0209.csv",sep = ",")
sampling_cluster.inf <- read.csv('/home/gpfs/home/wulab15/MacacaPFC/species_merge/HInN_sampling_cluster_inf_0209.csv')

table(sampling_cluster.inf$subcluster)
(1/B)*table(sampling_cluster.inf$sample) 
table(HInN@meta.data$sample)
subset_cellnames <- rownames(sampling_cluster.inf)
HInN <- SubsetData(HInN,cells = subset_cellnames)
table(HInN@meta.data$sample)
length(HInN@meta.data$sample)

objMERGE <- merge(x=HInN,y=MInN)
objMERGE@meta.data$seurat_clusters=NULL

head(objMERGE@meta.data)

DefaultAssay(objMERGE)<-'RNA'
splitlist <- SplitObject(objMERGE,split.by='species')
for (i in 1:length(splitlist)) {
  splitlist[[i]] <- NormalizeData(splitlist[[i]])
  splitlist[[i]] <- FindVariableFeatures(splitlist[[i]],selection.method = 'vst',nfeatures = 3000)
}

splitanchors <- FindIntegrationAnchors(object.list = splitlist,anchor.features = 3000)
newMERGE <- IntegrateData(anchorset = splitanchors)

DefaultAssay(newMERGE) <- 'integrated'
newMERGE <- ScaleData(newMERGE,vars.to.regress = c("percent.mt",'nCount_RNA')) 
newMERGE <- RunPCA(newMERGE, npcs = 100)
ElbowPlot(newMERGE,ndims = 100, reduction = "pca")

batch_info <- as.factor(newMERGE@meta.data$sample)
table(batch_info)
pca_coord <- newMERGE@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
newMERGE@reductions$pca@cell.embeddings <- pca_corrected$corrected

newMERGE <- FindNeighbors(newMERGE, dims = 1:25)
newMERGE <- FindClusters(newMERGE, resolution = 2)
newMERGE <- RunUMAP(newMERGE,dims = 1:25, reduction='pca')
DimPlot(newMERGE,label=TRUE)

Idents(newMERGE)<-'species'
DimPlot(newMERGE,cols=c("#3D578E","plum3"))+facet_wrap(~ident)
DimPlot(newMERGE,cols=c("#3D578E","plum3"))
saveRDS(newMERGE,file='/home/gpfs/home/wulab15/species_merge/HMintegrate_InN_subcluster.rds')
