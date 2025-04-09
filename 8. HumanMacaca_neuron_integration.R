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
MacacaALL<-readRDS(file="/~/PFCdata/MacacaPFC_withsubcluster_0103.rds")
HumanALL<-readRDS(file="/~/PFCdata/HumanPFC_withsubcluster_0103.rds")
Idents(HumanALL)<-'celltype'
Idents(MacacaALL)<-'celltype'
HExN<-subset(HumanALL,idents=c("Excitatory Neuron"))
MExN<-subset(MacacaALL,idents=c("Excitatory Neuron"))

#downsampling MExN data
dscellnum<-50000
dspct=dscellnum/length(MExN@meta.data$orig.ident) 
cluster_inf <- MExN@meta.data[,c("sample","subtype")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subtype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$subtype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subtype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
MExN<- subset(MExN,cells = subset_cellnames)

#downsampling HExN data
dspct=dscellnum/length(HExN@meta.data$orig.ident) 
cluster_inf <- HExN@meta.data[,c("sample","subtype")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subtype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$subtype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subtype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
HExN<- subset(HExN,cells = subset_cellnames)

objMERGE <- merge(x=HExN,y=MExN)

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
newMERGE <- RunUMAP(newMERGE,dims = 1:25, reduction='pca')
newMERGE <- FindNeighbors(newMERGE, dims = 1:25)
newMERGE <- FindClusters(newMERGE, resolution = 2)

DimPlot(newMERGE,label=TRUE)
Idents(newMERGE)<-'species'
DimPlot(newMERGE,cols=c("#3D578E","plum3"))

################HM InN

MacacaALL<-readRDS(file="/~/PFCdata/MacacaPFC_withsubcluster_0103.rds")
HumanALL<-readRDS(file="/~/PFCdata/HumanPFC_withsubcluster_0103.rds")
length(MacacaALL@meta.data$celltype)
length(HumanALL@meta.data$celltype)
Idents(HumanALL)<-'celltype'
Idents(MacacaALL)<-'celltype'
HInN<-subset(HumanALL,idents=c("Inhibitory Neuron"))
MInN<-subset(MacacaALL,idents=c("Inhibitory Neuron"))

#downsampling MInN data
dscellnum<-20000
dspct=dscellnum/length(MInN@meta.data$orig.ident) 
cluster_inf <- MInN@meta.data[,c("sample","subtype")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subtype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$subtype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subtype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
MInN<- subset(MInN,cells = subset_cellnames)

#downsampling HInN data
dspct=dscellnum/length(HInN@meta.data$orig.ident) 
cluster_inf <- HInN@meta.data[,c("sample","subtype")]
class(cluster_inf)
cluster_inf <- cluster_inf[order(cluster_inf$subtype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$subtype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('subtype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
HInN<- subset(HInN,cells = subset_cellnames)

objMERGE <- merge(x=HInN,y=MInN)

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
newMERGE <- RunUMAP(newMERGE,dims = 1:25, reduction='pca')
newMERGE <- FindNeighbors(newMERGE, dims = 1:25)
newMERGE <- FindClusters(newMERGE, resolution = 2)

DimPlot(newMERGE,label=TRUE)
Idents(newMERGE)<-'species'
DimPlot(newMERGE,cols=c("#3D578E","plum3"))
