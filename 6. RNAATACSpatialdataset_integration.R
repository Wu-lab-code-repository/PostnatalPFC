library(Seurat)
library(SeuratDisk)
library(Signac)
library(batchelor)
library(ggplot2)
library(sampling)

RNAdata<-readRDS(file='/~/PFCdata/RNAHuman_data.rds')
ATACdata<-readRDS(file='/~/PFCdata/ATACHuman_data.rds')
Xeniumdata<-readRDS(file='/~/PFCdata/XeniumHuman_data.rds')

######downsampledata
dscellnum<-10000
dspct=dscellnum/length(RNAdata@meta.data$orig.ident)
cluster_inf <- RNAdata@meta.data[,c("sample","celltype")]
cluster_inf <- cluster_inf[order(cluster_inf$celltype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$celltype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('celltype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
RNAdata <- subset(RNAdata,cells = subset_cellnames)

dspct=dscellnum/length(ATACdata@meta.data$orig.ident)
cluster_inf <- ATACdata@meta.data[,c("sample","celltype")]
cluster_inf <- cluster_inf[order(cluster_inf$celltype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$celltype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('celltype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
ATACdata <- subset(ATACdata,cells = subset_cellnames)

dspct=dscellnum/length(Xeniumdata@meta.data$orig.ident)
cluster_inf <- Xeniumdata@meta.data[,c("sample","celltype")]
cluster_inf <- cluster_inf[order(cluster_inf$celltype),]
k2 <- as.numeric(round(dspct*table(cluster_inf$celltype)))
sampling_cluster.inf <- strata(cluster_inf , stratanames=c('celltype'), size=k2, description=FALSE)
pos <- as.integer(rownames(sampling_cluster.inf))
sampling_cluster.inf <- cluster_inf[pos,]
subset_cellnames <- rownames(sampling_cluster.inf)
Xeniumdata <- subset(Xeniumdata,cells = subset_cellnames)

#########integration
allmerge<-merge(RNAdata,c(ATACdata,Xeniumdata))
genes<-c(intersect(c(VariableFeatures(RNAdata)),c(rownames(ATACdata@assays$RNA))))
splitlist <- SplitObject(allmerge,split.by='tech')
splitanchors <- FindIntegrationAnchors(object.list = splitlist,anchor.features =genes,k.anchor =30 )
newMERGE <- IntegrateData(anchorset = splitanchors)

DefaultAssay(newMERGE) <- 'integrated'
newMERGE <- ScaleData(newMERGE)
newMERGE <- RunPCA(newMERGE, npcs = 100)
ElbowPlot(newMERGE,ndims = 100, reduction = "pca")

newMERGE <- FindNeighbors(newMERGE, dims = 1:20)
newMERGE <- FindClusters(newMERGE, resolution = 1)
newMERGE <- RunUMAP(newMERGE,dims = 1:20, reduction='pca')

Idents(newMERGE)<-'tech'
DimPlot(newMERGE,label=TRUE)
Idents(newMERGE)<-'celltype'
DimPlot(newMERGE,label=TRUE)



