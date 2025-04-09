library(Seurat)
library(ggplot2)

XY1M6.obj<-LoadXenium(data.dir='/~/PFCdata/output_Y1M6/',fov='fov')
XY1M6<-subset(XY1M6.obj,subset = nCount_Xenium > 0)
XY1M6<-SCTransform(XY1M6,assay='Xenium')
XY1M6<-RunPCA(XY1M6,npcs=100,features=rownames(XY1M6))
XY1M6<-RunUMAP(XY1M6,dims=1:50)
XY1M6<-FindNeighbors(XY1M6,reduction='pca',dims=1:50)
XY1M6<-FindClusters(XY1M6,resolution=c(1,1.2,1.5,2,5))
DimPlot(XY1M6,reduction = 'umap',label=TRUE)

ImageFeaturePlot(XY1M6,features = c("nFeature_Xenium","nCount_Xenium"),dark.background = F)
FeaturePlot(XY1M6,features = c("nFeature_Xenium","nCount_Xenium"),order=TRUE)

DefaultAssay(XY1M6)<-'SCT'

FeaturePlot(XY1M6,features = c("AQP4","GFAP","ALDH1L1","SLC1A2",
                               "C3","CX3CR1","PTPRC","PDGFRA","OPALIN","TCF7L2",
                               "MOG","MAL","MBP","SLC17A7","GAD1","GAD2","SLC32A1",'CLDN5'),order=TRUE)

Idents(XY1M6) <- 'seurat_clusters'
cluster.ids <- c(rep("Astrocyte",14),
                 rep("OPC",8),
                 rep("Oligodendrocyte",20),
                 rep("Microglia",7),
                 rep("Inhibitory Neuron",12),
                 rep("Excitatory Neuron",17),                 
                 rep("Endothelial Cell",14),
                 rep("Low Quality",2))
names(cluster.ids) <- c(71,14,15,84,72,45,21,66,22,3,30,26,0,52,
                        37,68,40,85,82,48,44,90,
                        60,27,88,1,4,47,59,6,5,35,24,18,8,17,89,36,64,33,13,81,
                        16,28,70,54,63,11,25,
                        58,80,29,61,79,93,69,65,92,57,42,83,
                        73,2,62,43,20,32,78,34,23,74,75,19,67,49,39,76,87,
                        55,50,77,7,12,56,10,38,51,86,91,53,46,41,
                        9,31)

XY1M6 <- RenameIdents(XY1M6, cluster.ids)
XY1M6@meta.data$celltype <- Idents(XY1M6)
XY1M6@meta.data$celltype<-factor(XY1M6@meta.data$celltype,levels = c("Astrocyte","OPC","Oligodendrocyte","Microglia",
                                                                     "Inhibitory Neuron","Excitatory Neuron","Endothelial Cell","Low Quality"))
Idents(XY1M6)<-'celltype'
DimPlot(XY1M6,label=TRUE,raster = F)+NoLegend()


################imputation
RNAY1M6<-readRDS(file='/~/PFCdata/Human_snRNA_Y1M6around.rds')
DefaultAssay(RNAY1M6)<-'RNA'
RNAY1M6<-SCTransform(RNAY1M6)
XY1M6<-readRDS(file='/~/PFCdata/XeniumHuman_Y1M6.rds')
ImageDimPlot(XY1M6)

DefaultAssay(XY1M6)<-'SCT'
DefaultAssay(RNAY1M6)<-'SCT'
genes<-c(rownames(RNAY1M6@assays$SCT@data))

anchors1 <- FindTransferAnchors(reference = RNAY1M6, query = XY1M6,
                                reference.assay = 'SCT',query.assay = 'SCT',normalization.method = 'SCT',
                                dims = 1:20 ,reference.reduction = "pca",k.anchor = 50,features = genes,
                                recompute.residuals = F)

predictions1 <- TransferData(anchorset = anchors1, refdata =RNAY1M6@assays$SCT@data,query=XY1M6,
                             dims = 1:10,prediction.assay=TRUE)

XY1M6[['predSCT']]<-predictions1@assays$id

DefaultAssay(XY1M6)<-'predSCT'
FeaturePlot(XY1M6,features = c("AQP4","GFAP","ALDH1L1","SLC1A2","PDGFRA","OPALIN","TCF7L2",
                               "MOG","MAL","MBP","CUX1","CUX2","COL5A2","POU6F2","RORB","ETV1","SULF2","SLC38A11",
                               "FEZF2","FOXP2","HS3ST4","OPRK1","THEMIS","DCSTAMP","HTR2C",
                               "OPRK1","NR4A2","TLE4","SLC30A3","POU3F1","RIT2","SLIT3","CPLX1","CPLX3",
                               "LAMP5","RELN","PVALB","SST","VIP","ID2","SV2C"),order=TRUE)
Idents(XY1M6) <- 'seurat_clusters'
cluster.ids <- c(rep("AST_FB",6),
                 rep("AST_PP",8),
                 rep("OPC",8),
                 rep("NOL",3),
                 rep("MOL",17),
                 rep("Microglia",7),                 
                 rep("InN_PV",3),rep("InN_SST",4),rep("InN_VIP",2),rep("InN_ID2",3),                 
                 rep("ExN_L2_IT",2), rep("ExN_L2/3_IT",2), rep("ExN_L3/4_IT",1), 
                 rep("ExN_L4_IT",1), rep("ExN_L4/5_IT",1), rep("ExN_L5_IT",2), 
                 rep("ExN_L5_PT",1), rep("ExN_L5/6_NP",1), rep("ExN_L5/6_IT",1), 
                 rep("ExN_L6_IT",2), rep("ExN_L6_CT",1), rep("ExN_L6b",2), 
                 rep("Endothelial Cell",14),
                 rep("Low Quality",2))
names(cluster.ids) <- c(22,3,30,26,0,52,
                        71,14,15,84,72,45,21,66,
                        37,68,40,85,82,48,44,90,
                        60,27,88,
                        1,4,47,59,6,5,35,24,18,8,17,89,36,64,33,13,81,
                        16,28,70,54,63,11,25,                        
                        57,42,83,
                        61,69,65,92,
                        80,29,
                        79,93,58,                        
                        73,2,62,43,20,32,78,34,23,87,76,74,75,19,39,67,49,                        
                        55,50,77,7,12,56,10,38,51,86,91,53,46,41,
                        9,31)

XY1M6 <- RenameIdents(XY1M6, cluster.ids)
XY1M6@meta.data$subtype <- Idents(XY1M6)
Idents(XY1M6)<-'subtype'
DimPlot(XY1M6)
