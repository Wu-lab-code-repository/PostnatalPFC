library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)


RNAY8 <- Read10X(data_dir <- "/~/PFCdata/RNAY8/outs/filtered_feature_bc_matrix/")
write.table(as.data.frame(RNAY8),sep='\t',quote=F,file="/~/PFCdata/RNAY8/RNAY8_raw.txt")

RNAY8doublets<-read.csv("/~/PFCdata/RNAY8/RNAY8_doublet_position.csv",header=F)
pos <- grep('^1$',RNAY8doublets$V1)
RNAY8_noDoublets.data <- RNAY8[,-pos]
RNAY8 <- CreateSeuratObject(counts = RNAY8_noDoublets.data, project = "HumanPFC", min.cells = 5)
RNAY8@meta.data['sample_name']<-rownames(RNAY8@meta.data)
RNAY8[["percent.mt"]] <- PercentageFeatureSet(RNAY8, pattern ="^MT-" )
VlnPlot(RNAY8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

RNAY8 <- subset(RNAY8, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 5)
RNAY8 <- NormalizeData(RNAY8, normalization.method = "LogNormalize", scale.factor = 10000)
RNAY8 <- FindVariableFeatures(RNAY8, selection.method = "vst", nfeatures = 3000)
RNAY8 <- ScaleData(RNAY8)
RNAY8 <- RunPCA(RNAY8, features = VariableFeatures(object = RNAY8),npcs = 100)
ElbowPlot(RNAY8, ndims = 100, reduction = "pca")
RNAY8 <- RunUMAP(RNAY8, dims = 1:55)
RNAY8 <- FindNeighbors(RNAY8, dims = 1:55)
RNAY8 <- FindClusters(RNAY8, resolution = c(0.5,1,2,5,10))

DimPlot(RNAY8, reduction = "umap",label=TRUE)

FeaturePlot(RNAY8,features = c("AQP4","GFAP","ALDH1L1","SLC1A2",
                                "C3","CX3CR1","PTPRC","PDGFRA","OPALIN","TCF7L2",
                                "MOG","MAL","MBP","SLC17A7","GAD1","GAD2","SLC32A1","FN1","CLDN5"),order=TRUE)
FeaturePlot(RNAY8, features = c("CUX1","CUX2","COL5A2","POU6F2","RORB","ETV1","SULF2","SLC38A11",
                                 "FEZF2","FOXP2","HS3ST4","OPRK1","THEMIS","DCSTAMP","HTR2C",
                                 "OPRK1","NR4A2","TLE4","SLC30A3","POU3F1","RIT2","SLIT3","CPLX1","CPLX3",
                                 "LAMP5","RELN","PVALB","SST","VIP","ID2","SV2C"),order=TRUE)
            
            
            
            