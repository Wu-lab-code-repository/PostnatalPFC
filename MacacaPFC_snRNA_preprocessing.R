
##########Macaca_snRNA_postnatal
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

MacacaRNA_postnatal.data <- Read10X(data_dir <- "/home/gpfs/wulab15/MacacaPFC/snRNA/postnatal/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/MacacaPFC/snRNA/postnatal/raw/")
samples<-c('PFC_P0_1','PFC_P0_2','PFC_P0_3',
           'PFC_M6_1','PFC_M6_2',
           'PFC_Y1_1','PFC_Y1_2',
           'PFC_Y2_1','PFC_Y2_2',
           'PFC_Y4_1','PFC_Y4_2')
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,colnames(MacacaRNA_postnatal.data))
  print(samples[i])
  print(length(posi))
  PFCi <- MacacaRNA_postnatal.data[,posi]
  file_name =paste0("MacacaRNA_postnatal_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}



ibrary(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

MacacaRNA_postnatal.data <- Read10X(data_dir <- "/home/gpfs/wulab15/MacacaPFC/snRNA/postnatal/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/MacacaPFC/snRNA/postnatal/raw/")
samples<-c('PFC_P0_1','PFC_P0_2','PFC_P0_3',
           'PFC_M6_1','PFC_M6_2',
           'PFC_Y1_1','PFC_Y1_2',
           'PFC_Y2_1','PFC_Y2_2',
           'PFC_Y4_1','PFC_Y4_2')
cell_doublets<-c()
for(i in 1:length(samples)){
  findi<-paste0('-',i,'$')
  posi<-grep(findi,colnames(MacacaRNA_postnatal.data))
  print(samples[i])
  print(length(posi))
  PFCi <- MacacaRNA_postnatal.data[,posi]
  file_name =paste0("MacacaRNA_postnatal_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_doublets <- colnames(PFCi)[posi]
  cell_doublets<-c(cell_doublets,PFCi_doublets)
}
length(cell_doublets)
TOTALpos <- which(colnames(MacacaRNA_postnatal.data)%in%cell_doublets)
length(TOTALpos)
write.table(as.data.frame(cell_doublets),sep='\t',quote=F,file='MacacaRNA_postnatal_cell_doublet.txt')

MacacaRNA_postnatal_NEW <- MacacaRNA_postnatal.data[,-TOTALpos]
dim(MacacaRNA_postnatal_NEW) 

Macaca_postnatal <- CreateSeuratObject(counts = MacacaRNA_postnatal_NEW, project = "Macaca_postnatal", min.cells = 10)
sample_name <- rownames(Macaca_postnatal@meta.data)
Macaca_postnatal@meta.data['sample_name']<-sample_name
Macaca_postnatal@meta.data['sample_Macaca']<-NA
Macaca_postnatal@meta.data['sample_age']<-NA
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,Macaca_postnatal@meta.data$sample_name)
  Macaca_postnatal@meta.data$sample_Macaca[posi]<-as.character(samples[i])
  Macaca_postnatal@meta.data$sample_age[posi]<-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))}
table(Macaca_postnatal@meta.data$sample_Macaca)
table(Macaca_postnatal@meta.data$sample_age)

head(Macaca_postnatal@meta.data)
tail(Macaca_postnatal@meta.data)
mtgenelist <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
mtgenelist <- intersect(mtgenelist,rownames(Macaca_postnatal))
Macaca_postnatal[["percent.mt"]] <- PercentageFeatureSet(Macaca_postnatal, features = mtgenelist)

VlnPlot(Macaca_postnatal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

Macaca_postnatal$sample_Macaca <- factor(x=Macaca_postnatal$sample_Macaca,levels= c(samples))
Idents(Macaca_postnatal) <- 'sample_Macaca'
VlnPlot(Macaca_postnatal, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(Macaca_postnatal, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(Macaca_postnatal, features = c("percent.mt"),pt.size=0.1)

Macaca_postnatal <- subset(Macaca_postnatal, subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 5)
table(Macaca_postnatal@meta.data$sample_Macaca)

Macaca_postnatal <- NormalizeData(Macaca_postnatal, normalization.method = "LogNormalize", scale.factor = 10000)
Macaca_postnatal <- FindVariableFeatures(Macaca_postnatal, selection.method = "vst", nfeatures = 3000)

Macaca_postnatal <- ScaleData(Macaca_postnatal, vars.to.regress = c("percent.mt","nCount_RNA"))
Macaca_postnatal <- RunPCA(Macaca_postnatal, features = VariableFeatures(object = Macaca_postnatal),npcs = 100)
ElbowPlot(Macaca_postnatal, ndims = 100, reduction = "pca")
batch_info <- as.factor(Macaca_postnatal@meta.data$sample_Macaca)
table(batch_info)
pca_coord <- Macaca_postnatal@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
Macaca_postnatal@reductions$pca@cell.embeddings <- pca_corrected$corrected

Macaca_postnatal <- FindNeighbors(Macaca_postnatal, dims = 1:50)
Macaca_postnatal <- FindClusters(Macaca_postnatal, resolution = 5)
Macaca_postnatal <- RunUMAP(Macaca_postnatal, dims = 1:50)
DimPlot(Macaca_postnatal, reduction = "umap",label=TRUE)
saveRDS(Macaca_postnatal, file = "/home/gpfs/wulab15/MacacaPFC/snRNA/Macaca_postnatal_afterQC.rds")


##########Macaca_snRNA_embryo
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

MacacaRNA_embryo.data <- Read10X(data_dir <- "/home/gpfs/wulab15/MacacaPFC/snRNA/embryo/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/MacacaPFC/snRNA/embryo/raw/")
samples<-c('PFC_E90_1','PFC_E110_1')
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,colnames(MacacaRNA_embryo.data))
  print(samples[i])
  print(length(posi))
  PFCi <- MacacaRNA_embryo.data[,posi]
  file_name =paste0("MacacaRNA_embryo_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}


ibrary(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

MacacaRNA_embryo.data <- Read10X(data_dir <- "/home/gpfs/wulab15/MacacaPFC/snRNA/embryo/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/MacacaPFC/snRNA/embryo/raw/")
samples<-c('PFC_E90_1','PFC_E110_1')
cell_doublets<-c()
for(i in 1:length(samples)){
  findi<-paste0('-',i,'$')
  posi<-grep(findi,colnames(MacacaRNA_embryo.data))
  print(samples[i])
  print(length(posi))
  PFCi <- MacacaRNA_embryo.data[,posi]
  file_name =paste0("MacacaRNA_embryo_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_doublets <- colnames(PFCi)[posi]
  cell_doublets<-c(cell_doublets,PFCi_doublets)
}
length(cell_doublets)
TOTALpos <- which(colnames(MacacaRNA_embryo.data)%in%cell_doublets)
length(TOTALpos)
write.table(as.data.frame(cell_doublets),sep='\t',quote=F,file='MacacaRNA_embryo_cell_doublet.txt')

MacacaRNA_embryo_NEW <- MacacaRNA_embryo.data[,-TOTALpos]
dim(MacacaRNA_embryo_NEW) 

Macaca_embryo <- CreateSeuratObject(counts = MacacaRNA_embryo_NEW, project = "Macaca_embryo", min.cells = 2)
sample_name <- rownames(Macaca_embryo@meta.data)
Macaca_embryo@meta.data['sample_name']<-sample_name
Macaca_embryo@meta.data['sample_Macaca']<-NA
Macaca_embryo@meta.data['sample_age']<-NA
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,Macaca_embryo@meta.data$sample_name)
  Macaca_embryo@meta.data$sample_Macaca[posi]<-as.character(samples[i])
  Macaca_embryo@meta.data$sample_age[posi]<-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))}
table(Macaca_embryo@meta.data$sample_Macaca)
table(Macaca_embryo@meta.data$sample_age)

head(Macaca_embryo@meta.data)
tail(Macaca_embryo@meta.data)
mtgenelist=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
mtgenelist=intersect(mtgenelist,rownames(Macaca_embryo))
Macaca_embryo[["percent.mt"]] <- PercentageFeatureSet(Macaca_embryo, features = mtgenelist)

VlnPlot(Macaca_embryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
dev.off()
length(which(Macaca_embryo@meta.data$percent.mt > 5))  
length(which(Macaca_embryo@meta.data$nFeature_RNA > 7000))  
length(which(Macaca_embryo@meta.data$nFeature_RNA < 800))  

Macaca_embryo$sample_Macaca=factor(x=Macaca_embryo$sample_Macaca,levels= c(samples))
Idents(Macaca_embryo) <-'sample_Macaca'
VlnPlot(Macaca_embryo, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(Macaca_embryo, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(Macaca_embryo, features = c("percent.mt"),pt.size=0.1)

Macaca_embryo <- subset(Macaca_embryo, subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 5)
table(Macaca_embryo@meta.data$sample_Macaca)

Macaca_embryo <- NormalizeData(Macaca_embryo, normalization.method = "LogNormalize", scale.factor = 10000)
Macaca_embryo <- FindVariableFeatures(Macaca_embryo, selection.method = "vst", nfeatures = 3000)

Macaca_embryo <- ScaleData(Macaca_embryo, vars.to.regress = c("percent.mt","nCount_RNA"))
Macaca_embryo <- RunPCA(Macaca_embryo, features = VariableFeatures(object = Macaca_embryo),npcs = 100)
ElbowPlot(Macaca_embryo, ndims = 100, reduction = "pca")
batch_info <- as.factor(Macaca_embryo@meta.data$sample_Macaca)
table(batch_info)
pca_coord <- Macaca_embryo@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
Macaca_embryo@reductions$pca@cell.embeddings <- pca_corrected$corrected

Macaca_embryo <- FindNeighbors(Macaca_embryo, dims = 1:50)
Macaca_embryo <- FindClusters(Macaca_embryo, resolution = 5)
Macaca_embryo <- RunUMAP(Macaca_embryo, dims = 1:50)
DimPlot(Macaca_embryo, reduction = "umap",label=TRUE)
saveRDS(Macaca_embryo, file = "/home/gpfs/wulab15/MacacaPFC/snRNA/Macaca_embryo_afterQC.rds")




###########mergeMacacaRNA EP
MacacaE<-readRDS(file="/home/gpfs/wulab15/MacacaPFC/snRNA/Macaca_embryo_afterQC.rds")
MacacaP<-readRDS(file="/home/gpfs/wulab15/MacacaPFC/snRNA/Macaca_postnatal_afterQC.rds")
head(MacacaE@meta.data)
head(MacacaP@meta.data)

MacacaE@meta.data['time']<-'embryo'
MacacaP@meta.data['time']<-'postnatal'

combineddata<-merge(MacacaE,MacacaP)
splitlist <- SplitObject(combineddata,split.by='time')
for (i in 1:length(splitlist)) {
  splitlist[[i]] <- NormalizeData(splitlist[[i]])
  splitlist[[i]] <- FindVariableFeatures(splitlist[[i]],selection.method = 'vst',nfeatures = 3000)
}

splitanchors <- FindIntegrationAnchors(object.list = splitlist,anchor.features = 3000)
dataMERGE <- IntegrateData(anchorset = splitanchors)

DefaultAssay(dataMERGE) <- 'integrated'
dataMERGE <- ScaleData(dataMERGE)
dataMERGE <- RunPCA(dataMERGE, npcs = 100)
ElbowPlot(dataMERGE,ndims = 100, reduction = "pca")

batch_info <- as.factor(dataMERGE@meta.data$sample_Macaca)
table(batch_info)
pca_coord <- dataMERGE@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
dataMERGE@reductions$pca@cell.embeddings <- pca_corrected$corrected

dataMERGE <- FindNeighbors(dataMERGE, dims = 1:50)
dataMERGE <- FindClusters(dataMERGE, resolution = 5)
dataMERGE <- RunUMAP(dataMERGE,dims = 1:50, reduction='pca')
DimPlot(dataMERGE,label=TRUE)

Idents(dataMERGE)<-'time'
DimPlot(dataMERGE)

saveRDS(dataMERGE, file = "/home/gpfs/wulab15/MacacaPFC/snRNA/MacacaRNA_EPmerge.rds")


###############subtypeUMAP
MacacaEP<-readRDS(file="/home/gpfs/home/wulab15/MacacaPFC/202305/combinedMacacaDatadim50.rds")
MacacaEP@meta.data$subtype<-factor(MacacaEP@meta.data$subtype,levels=c('AST_FB', 'AST_PP','OPC','NOL','MOL','Microglia','Neural Progenitor Cell',
                                                                       'InN_PV', 'InN_SST','InN_VIP','InN_ID2',
                                                                       'ExN_L2_IT','ExN_L2/3_IT','ExN_L3/4_IT','ExN_L4_IT',
                                                                       'ExN_L4/5_IT', 'ExN_L5_IT', 'ExN_L5_PT','ExN_L5/6_NP', 'ExN_L5/6_IT', 
                                                                       'ExN_L6_IT','ExN_L6_CT', 'ExN_L6b',
                                                                       'Endothelial Cell'))
Idents(MacacaEP)<-'subtype'
DimPlot(MacacaEP,label = T,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312","grey",
                                  "#044D2B","palegreen4","olivedrab","forestgreen",
                                  "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                  "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                  "#9C964A"))

Idents(MacacaEP)<-'subtype'
p<-DimPlot(MacacaEP,label = F,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312","grey",
                                     "#044D2B","palegreen4","olivedrab","forestgreen",
                                     "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                     "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                     "#9C964A"),raster=F)
pdf("/home/gpfs/home/wulab15/Figures/UMAP/MacacaAllEP_Subtype_UMAP.pdf",height = 4.5,width = 8)
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/UMAP/MacacaAllEP_Subtype_UMAP.tiff",height = 450,width = 700)
p
dev.off()









