##########Human_snRNA_embryo
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/cellranger_outs/embryo/")
samples<-c('PFC_GW16_1','PFC_GW22_1','PFC_GW26_1')
for(i in 1:length(samples)){
  matrixpath <-paste0("/home/gpfs/home/wulab15/HumanPFC/snRNA/cellranger_outs/embryo/",samples[i],"_outs/filtered_feature_bc_matrix/")
  PFCi <- Read10X(data_dir <- matrixpath)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/raw/")
  file_name =paste0("HumanRNA_embryo_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}

library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/raw/")
samples<-c('PFC_GW16_1','PFC_GW22_1','PFC_GW26_1')
for(i in 1:length(samples)){
  matrixpath <-paste0("/home/gpfs/home/wulab15/HumanPFC/snRNA/cellranger_outs/embryo/",samples[i],"_outs/filtered_feature_bc_matrix/")
  PFCi <- Read10X(data_dir <- matrixpath)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/raw/")
  file_name =paste0("HumanRNA_embryo_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_noDoublets.data <- PFCi[,-posi]
  PFCsamplei <- CreateSeuratObject(counts = PFCi_noDoublets.data, project = "HumanPFC", min.cells = 1)
  PFCsamplei@meta.data['sample_Human'] <-as.character(samples[i])
  PFCsamplei@meta.data['sample_age'] <-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/rds/")
  rdsfilename<-paste0("HumanRNA_embryo_",samples[i],".rds")
  saveRDS(PFCsamplei,file= rdsfilename)
}

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/rds/")
Hdata1<-readRDS(file = 'HumanRNA_embryo_PFC_GW16_1.rds')
for(i in 2:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/embryo/rds/")
  rdsfilename<-paste0("HumanRNA_embryo_",samples[i],".rds")
  PFCsamplei <-readRDS(file=rdsfilename)
  Hdata1<-merge(Hdata1,PFCsamplei)}

Human_embryo <- CreateSeuratObject(counts = Hdata1, project = "Human_embryo", min.cells = 10)
dim(Human_embryo) 
head(Human_embryo@meta.data)
Human_embryo@meta.data['sample_name']<-rownames(Human_embryo@meta.data)
Human_embryo[["percent.mt"]] <- PercentageFeatureSet(Human_embryo, pattern ="^MT-" )
VlnPlot(Human_embryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

Human_embryo$sample_Human=factor(x=Human_embryo$sample_Human,levels= c(samples))
Idents(Human_embryo) <-'sample_Human'
VlnPlot(Human_embryo, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(Human_embryo, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(Human_embryo, features = c("percent.mt"),pt.size=0.1)

Human_embryo <- subset(Human_embryo, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 5)
table(Human_embryo@meta.data$sample_Human)

Human_embryo <- NormalizeData(Human_embryo, normalization.method = "LogNormalize", scale.factor = 10000)
Human_embryo <- FindVariableFeatures(Human_embryo, selection.method = "vst", nfeatures = 3000)
Human_embryo <- ScaleData(Human_embryo, vars.to.regress = c("percent.mt","nCount_RNA"))
Human_embryo <- RunPCA(Human_embryo, features = VariableFeatures(object = Human_embryo),npcs = 100)
ElbowPlot(Human_embryo, ndims = 100, reduction = "pca")
batch_info <- as.factor(Human_embryo@meta.data$sample_Human)
table(batch_info)
pca_coord <- Human_embryo@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
Human_embryo@reductions$pca@cell.embeddings <- pca_corrected$corrected

Human_embryo <- FindNeighbors(Human_embryo, dims = 1:50)
Human_embryo <- FindClusters(Human_embryo, resolution = 5)
Human_embryo <- RunUMAP(HerringData, dims = 1:50)
DimPlot(Human_embryo, reduction = "umap",label=TRUE)
saveRDS(Human_embryo, file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/Human_embryo_afterQC.rds")


############Human_snRNA_postnatal
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

HumanRNA_postnatal.data <- Read10X(data_dir <- "/home/gpfs/wulab15/HumanPFC/snRNA/postnatal/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/postnatal/raw/")
samples<-c('PFC_Y1_1','PFC_Y3_1','PFC_Y4_1',
                             'PFC_Y4M11_1','PFC_Y6_1',
                             'PFC_Y8_1','PFC_Y9_1',
                             'PFC_Y10_1','PFC_Y12_1',
                             'PFC_Y13_1','PFC_Y14_1','PFC_Y15_1','PFC_Y16_1')
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,colnames(HumanRNA_postnatal.data))
  print(samples[i])
  print(length(posi))
  PFCi <- HumanRNA_postnatal.data[,posi]
  file_name =paste0("HumanRNA_postnatal_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}


library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

HumanRNA_postnatal.data <- Read10X(data_dir <- "/home/gpfs/wulab15/HumanPFC/snRNA/postnatal/outs_aggr/outs/filtered_feature_bc_matrix/")
setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/postnatal/raw/")
samples<-c('PFC_Y1_1','PFC_Y3_1','PFC_Y4_1',
                             'PFC_Y4M11_1','PFC_Y6_1',
                             'PFC_Y8_1','PFC_Y9_1',
                             'PFC_Y10_1','PFC_Y12_1',
                             'PFC_Y13_1','PFC_Y14_1','PFC_Y15_1','PFC_Y16_1')
cell_doublets<-c()
for(i in 1:length(samples)){
  findi<-paste0('-',i,'$')
  posi<-grep(findi,colnames(HumanRNA_postnatal.data))
  print(samples[i])
  print(length(posi))
  PFCi <- HumanRNA_postnatal.data[,posi]
  file_name =paste0("HumanRNA_postnatal_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_doublets <- colnames(PFCi)[posi]
  cell_doublets<-c(cell_doublets,PFCi_doublets)
}
length(cell_doublets)
TOTALpos <- which(colnames(HumanRNA_postnatal.data)%in%cell_doublets)
length(TOTALpos)
write.table(as.data.frame(cell_doublets),sep='\t',quote=F,file='HumanRNA_postnatal_cell_doublet.txt')

HumanRNA_postnatal_NEW <- HumanRNA_postnatal.data[,-TOTALpos]
dim(HumanRNA_postnatal_NEW) 

Human_postnatal <- CreateSeuratObject(counts =HumanRNA_postnatal_NEW, project = "Human_postnatal", min.cells = 10)
sample_name <- rownames(Human_postnatal@meta.data)
Human_postnatal@meta.data['sample_name']<-sample_name
Human_postnatal@meta.data['sample_Human']<-NA
Human_postnatal@meta.data['sample_age']<-NA
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,Human_postnatal@meta.data$sample_name)
  Human_postnatal@meta.data$sample_Human[posi]<-as.character(samples[i])
  Human_postnatal@meta.data$sample_age[posi]<-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))}
table(Human_postnatal@meta.data$sample_Human)
table(Human_postnatal@meta.data$sample_age)

head(Human_postnatal@meta.data)
tail(Human_postnatal@meta.data)

Human_postnatal[["percent.mt"]] <- PercentageFeatureSet(Human_postnatal, pattern ="^MT-" )
VlnPlot(Human_postnatal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

Human_postnatal$sample_Human=factor(x=Human_postnatal$sample_Human,levels= c(samples))
Idents(Human_postnatal) <-'sample_Human'
VlnPlot(Human_postnatal, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(Human_postnatal, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(Human_postnatal, features = c("percent.mt"),pt.size=0.1)

Human_postnatal <- subset(Human_postnatal, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 5)
table(Human_postnatal@meta.data$sample_Human)

Human_postnatal <- NormalizeData(Human_postnatal, normalization.method = "LogNormalize", scale.factor = 10000)
Human_postnatal <- FindVariableFeatures(Human_postnatal, selection.method = "vst", nfeatures = 3000)

Human_postnatal <- ScaleData(Human_postnatal, vars.to.regress = c("percent.mt","nCount_RNA"))
Human_postnatal <- RunPCA(Human_postnatal, features = VariableFeatures(object = Human_postnatal),npcs = 100)
ElbowPlot(Human_postnatal, ndims = 100, reduction = "pca")
batch_info <- as.factor(Human_postnatal@meta.data$sample_Human)
table(batch_info)
pca_coord <- Human_postnatal@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
Human_postnatal@reductions$pca@cell.embeddings <- pca_corrected$corrected

Human_postnatal <- FindNeighbors(Human_postnatal, dims = 1:50)
Human_postnatal <- FindClusters(Human_postnatal, resolution = 5)
Human_postnatal <- RunUMAP(Human_postnatal, dims = 1:50)
DimPlot(Human_postnatal, reduction = "umap",label=TRUE)
saveRDS(Human_postnatal, file = "/home/gpfs/wulab15/HumanPFC/snRNA/Human_postnatal_afterQC.rds")

###############Herring postnatal

library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

samples<-c('PFC_M1_1','PFC_M2_1','PFC_M3_1','PFC_M4_1','PFC_M6_1','PFC_M10_1',
           'PFC_Y1_2','PFC_1M6_1','PFC_Y2_1','PFC_Y4_2','PFC_Y6M6_1','PFC_Y10_2',
           'PFC_Y12_2','PFC_Y16_2','PFC_Y17_1','PFC_Y20_1','PFC_Y20_2','PFC_Y25_1')
for(i in 1:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/h5/")
  Herringh5name <-paste0("HumanRNA_postnatal_",samples[i],"_snRNAseq_filtered_feature_bc_matrix.h5")
  PFCi <- Read10X_h5(file=Herringh5name, use.names = TRUE, unique.features = TRUE)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/raw/")
  file_name =paste0("HumanRNA_postnatal_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}


library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/raw/")
samples<-c('PFC_M1_1','PFC_M2_1','PFC_M3_1','PFC_M4_1','PFC_M6_1','PFC_M10_1',
           'PFC_Y1_2','PFC_1M6_1','PFC_Y2_1','PFC_Y4_2','PFC_Y6M6_1','PFC_Y10_2',
           'PFC_Y12_2','PFC_Y16_2','PFC_Y17_1','PFC_Y20_1','PFC_Y20_2','PFC_Y25_1')
for(i in 1:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/h5/")
  Herringh5name <-paste0("HumanRNA_postnatal_",samples[i],"_snRNAseq_filtered_feature_bc_matrix.h5")
  PFCi <- Read10X_h5(file=Herringh5name, use.names = TRUE, unique.features = TRUE)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/raw/")
  file_name =paste0("HumanRNA_postnatal_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_noDoublets.data <- PFCi[,-posi]
  PFCsamplei <- CreateSeuratObject(counts = PFCi_noDoublets.data, project = "HerringData", min.cells = 1)
  PFCsamplei@meta.data['sample_Human'] <-as.character(samples[i])
  PFCsamplei@meta.data['sample_age'] <-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/rds/")
  rdsfilename<-paste0("HumanRNA_postnatal_",samples[i],".rds")
  saveRDS(PFCsamplei,file= rdsfilename)
}

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/rds/")
Hdata1<-readRDS(file = 'HumanRNA_postnatal_PFC_M1_1.rds')
for(i in 2:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/postnatal/rds/")
  rdsfilename<-paste0("HumanRNA_postnatal_",samples[i],".rds")
  PFCsamplei <-readRDS(file=rdsfilename)
  Hdata1<-merge(Hdata1,PFCsamplei)}

HerringData <- CreateSeuratObject(counts = Hdata1@assays$RNA@counts, project = "HerringData_postnatal", min.cells = 10,meta.data = Hdata1@meta.data)
dim(HerringData) 
head(HerringData@meta.data)
HerringData@meta.data['sample_name']<-rownames(HerringData@meta.data)
HerringData[["percent.mt"]] <- PercentageFeatureSet(HerringData, pattern ="^MT-" )
VlnPlot(HerringData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

HerringData$sample_Human=factor(x=HerringData$sample_Human,levels= c(samples))
Idents(HerringData) <-'sample_Human'
VlnPlot(HerringData, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(HerringData, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(HerringData, features = c("percent.mt"),pt.size=0.1)

HerringData <- subset(HerringData, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 5)
table(HerringData@meta.data$sample_Human)

HerringData <- NormalizeData(HerringData, normalization.method = "LogNormalize", scale.factor = 10000)
HerringData <- FindVariableFeatures(HerringData, selection.method = "vst", nfeatures = 3000)
HerringData <- ScaleData(HerringData, vars.to.regress = c("percent.mt","nCount_RNA"))
HerringData <- RunPCA(HerringData, features = VariableFeatures(object = HerringData),npcs = 100)
ElbowPlot(HerringData, ndims = 100, reduction = "pca")
batch_info <- as.factor(HerringData@meta.data$sample_Human)
table(batch_info)
pca_coord <- HerringData@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
HerringData@reductions$pca@cell.embeddings <- pca_corrected$corrected

HerringData <- FindNeighbors(HerringData, dims = 1:50)
HerringData <- FindClusters(HerringData, resolution = 5)
HerringData <- RunUMAP(HerringData, dims = 1:50)
DimPlot(HerringData, reduction = "umap",label=TRUE)
saveRDS(HerringData, file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/Human_postnatal_afterQC.rds")


#############Herring embryo
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

samples<-c('PFC_GW22_2','PFC_GW24_1','PFC_GW34_1')
for(i in 1:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/h5/")
  Herringh5name <-paste0("HumanRNA_embryo_",samples[i],"_snRNAseq_filtered_feature_bc_matrix.h5")
  PFCi <- Read10X_h5(file=Herringh5name, use.names = TRUE, unique.features = TRUE)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/raw/")
  file_name =paste0("HumanRNA_embryo_",samples[i],"_raw.txt")
  write.table(as.data.frame(PFCi),sep='\t',quote=F,file=file_name)
}


library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/raw/")
samples<-c('PFC_GW22_2','PFC_GW24_1','PFC_GW34_1')
for(i in 1:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/h5/")
  Herringh5name <-paste0("HumanRNA_embryo_",samples[i],"_snRNAseq_filtered_feature_bc_matrix.h5")
  PFCi <- Read10X_h5(file=Herringh5name, use.names = TRUE, unique.features = TRUE)
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/raw/")
  file_name =paste0("HumanRNA_embryo_",samples[i],"_doublet_position.csv")
  PFCi_pos<-read.csv(file_name,header=F)
  posi <- grep('^1$',PFCi_pos$V1)
  print(length(posi))
  PFCi_noDoublets.data <- PFCi[,-posi]
  PFCsamplei <- CreateSeuratObject(counts = PFCi_noDoublets.data, project = "HerringData", min.cells = 1)
  PFCsamplei@meta.data['sample_Human'] <-as.character(samples[i])
  PFCsamplei@meta.data['sample_age'] <-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/rds/")
  rdsfilename<-paste0("HumanRNA_embryo_",samples[i],".rds")
  saveRDS(PFCsamplei,file= rdsfilename)
}

setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/rds/")
Hdata1<-readRDS(file = 'HumanRNA_embryo_PFC_GW22_2.rds')
for(i in 2:length(samples)){
  setwd("/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/embryo/rds/")
  rdsfilename<-paste0("HumanRNA_embryo_",samples[i],".rds")
  PFCsamplei <-readRDS(file=rdsfilename)
  Hdata1<-merge(Hdata1,PFCsamplei)}

HerringData <- CreateSeuratObject(counts = Hdata1@assays$RNA@counts, project = "HerringData_embryo", min.cells = 10,meta.data = Hdata1@meta.data)
dim(HerringData) 
head(HerringData@meta.data)
HerringData@meta.data['sample_name']<-rownames(HerringData@meta.data)
HerringData[["percent.mt"]] <- PercentageFeatureSet(HerringData, pattern ="^MT-" )
VlnPlot(HerringData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

HerringData$sample_Human=factor(x=HerringData$sample_Human,levels= c(samples))
Idents(HerringData) <-'sample_Human'
VlnPlot(HerringData, features = c("nFeature_RNA"),pt.size=0.1)
VlnPlot(HerringData, features = c("nCount_RNA"),pt.size=0.1)
VlnPlot(HerringData, features = c("percent.mt"),pt.size=0.1)

HerringData <- subset(HerringData, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 5)
table(HerringData@meta.data$sample_Human)

HerringData <- NormalizeData(HerringData, normalization.method = "LogNormalize", scale.factor = 10000)
HerringData <- FindVariableFeatures(HerringData, selection.method = "vst", nfeatures = 3000)
HerringData <- ScaleData(HerringData, vars.to.regress = c("percent.mt","nCount_RNA"))
HerringData <- RunPCA(HerringData, features = VariableFeatures(object = HerringData),npcs = 100)
ElbowPlot(HerringData, ndims = 100, reduction = "pca")
batch_info <- as.factor(HerringData@meta.data$sample_Human)
table(batch_info)
pca_coord <- HerringData@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
HerringData@reductions$pca@cell.embeddings <- pca_corrected$corrected

HerringData <- FindNeighbors(HerringData, dims = 1:50)
HerringData <- FindClusters(HerringData, resolution = 5)
HerringData <- RunUMAP(HerringData, dims = 1:50)
DimPlot(HerringData, reduction = "umap",label=TRUE)
saveRDS(HerringData, file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/Human_embryo_afterQC.rds")



###########mergeHumanRNA EP
HumanE<-readRDS(file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/Human_embryo_afterQC.rds")
HumanP<-readRDS(file="/home/gpfs/wulab15/HumanPFC/snRNA/Human_postnatal_afterQC.rds")
HumanE_Herring<-readRDS(file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/Human_embryo_afterQC.rds")
HumanP_Herring<-readRDS(file = "/home/gpfs/home/wulab15/HumanPFC/snRNA/HerringData/Human_postnatal_afterQC.rds")

head(HumanE@meta.data)
head(HumanP@meta.data)

HumanE@meta.data['time']<-'Human_embryo'
HumanP@meta.data['time']<-'Human_postnatal'
HumanE_Herring@meta.data['time']<-'Herring_embryo'
HumanP_Herring@meta.data['time']<-'Herring_postnatal'

HumanEmbryo<- merge(HumanE,HumanE_Herring)
HumanPostnatal<- merge(HumanP,HumanP_Herring)

combineddata<-merge(HumanEmbryo,HumanPostnatal)
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

batch_info <- as.factor(dataMERGE@meta.data$sample_Human)
table(batch_info)
pca_coord <- dataMERGE@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
dataMERGE@reductions$pca@cell.embeddings <- pca_corrected$corrected

dataMERGE <- FindNeighbors(dataMERGE, dims = 1:55)
dataMERGE <- FindClusters(dataMERGE, resolution = 5)
dataMERGE <- RunUMAP(dataMERGE,dims = 1:55, reduction='pca')
DimPlot(dataMERGE,label=TRUE)

Idents(dataMERGE)<-'time'
DimPlot(dataMERGE)

saveRDS(dataMERGE, file = "/home/gpfs/wulab15/HumanPFC/HumanRNA_EPmerge.rds")

###############subtypeUMAP
HumanEP<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/202305/combinedHumanDatadim55.rds")

HumanEP@meta.data$subtype<-factor(HumanEP@meta.data$subtype,levels=c('AST_FB', 'AST_PP','OPC','NOL','MOL','Microglia','Neural Progenitor Cell',
                                                                     'InN_PV', 'InN_SST','InN_VIP','InN_ID2',
                                                                     'ExN_L2_IT','ExN_L2/3_IT','ExN_L3/4_IT','ExN_L4_IT',
                                                                     'ExN_L4/5_IT', 'ExN_L5_IT', 'ExN_L5_PT','ExN_L5/6_NP', 'ExN_L5/6_IT', 
                                                                     'ExN_L6_IT','ExN_L6_CT', 'ExN_L6b',
                                                                     'Endothelial Cell'))
Idents(HumanEP)<-'subtype'
DimPlot(HumanEP,label = T,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312","grey",
                                 "#044D2B","palegreen4","olivedrab","forestgreen",
                                 "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                 "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                 "#9C964A"))

Idents(HumanEP)<-'subtype'
p<-DimPlot(HumanEP,label = F,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312","grey",
                                    "#044D2B","palegreen4","olivedrab","forestgreen",
                                    "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                    "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                    "#9C964A"),raster=F)
pdf("/home/gpfs/home/wulab15/Figures/UMAP/HumanAllEP_Subtype_UMAP.pdf",height = 4.5,width = 8)
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/UMAP/HumanAllEP_Subtype_UMAP.tiff",height = 450,width = 700)
p
dev.off()

























