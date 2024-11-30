library(Seurat)
library(ggplot2)

samples<-c('P0','M6','Y1','Y2','Adult')
for(i in 1:length(samples)){
  outputpath<-paste0('/home/gpfs/wulab15/MacacaPFC/Xenium/outs/output-',samples[i],'/')
  Xeniumi <- LoadXenium(data.dir=outputpath,fov='fov')
  Xeniumi <- subset(Xeniumi,subset = nCount_Xenium > 0)
  Xeniumi <- SCTransform(Xeniumi,assay='Xenium')
  Xeniumi <- RunPCA(Xeniumi,npcs=100,features=rownames(Xeniumi))
  Xeniumi <- RunUMAP(Xeniumi,dims=1:30)
  Xeniumi <- FindNeighbors(Xeniumi,reduction='pca',dims=1:30)
  Xeniumi <- FindClusters(Xeniumi,resolution=0.3)
  Xeniumi@meta.data['sample_Macaca']<-as.character(samples[i])
  Xeniumi@meta.data['sample_names']<-paste0(samples[i],'_',rownames(samples[i]@meta.data))
  setwd("/home/gpfs/wulab15/MacacaPFC/Xenium/")
  filenamei<-paste0('MacacaXenium_',samples[i],'_obj.rds')
  saveRDS(Xeniumi,file=filenamei)
}


samples<-c('P0','M6','Y1','Y2','Adult')
for(i in 1:length(samples)){
  setwd("/home/gpfs/wulab15/MacacaPFC/Xenium/")
  filenamei<- paste0('MacacaXenium_',samples[i],'_obj.rds')
  Xeniumi <- readRDS(file=filenamei)
  setwd("/home/gpfs/wulab15/MacacaPFC/snRNA/SCT/")
  filenamei<- paste0('MacacaRNA_',samples[i],'_SCT.rds')
  RNAi <- readRDS(file=filenamei)
  DefaultAssay(Xeniumi)<-'SCT'
  DefaultAssay(RNAi)<-'SCT'
  genes<-c(rownames(RNAi@assays$SCT@data))
  anchors <- FindTransferAnchors(reference = RNAi, query = Xeniumi,
                                reference.assay = 'SCT',query.assay = 'SCT',normalization.method = 'SCT',#reduction = 'cca',
                                dims = 1:30 ,reference.reduction = "pca",k.anchor = 50,features = genes,
                                recompute.residuals = F)
  predictions <- TransferData(anchorset = anchors, refdata =RNAi@assays$SCT@data,query=Xeniumi,
                             dims = 1:30,prediction.assay=TRUE)
  Xeniumi[['predSCT']]<-predictions@assays$id
  setwd("/home/gpfs/wulab15/MacacaPFC/Xenium/")
  filenamei<-paste0('MacacaXenium_',samples[i],'_predSCT.rds')
  saveRDS(Xeniumi,file=filenamei)
}

MXP0<-readRDS('/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_P0_predSCT.rds')
MXM6<-readRDS('/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_M6_predSCT.rds')
MXY1<-readRDS('/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_Y1_predSCT.rds')
MXY2<-readRDS('/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_Y2_predSCT.rds')
MXYA<-readRDS('/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_Adult_predSCT.rds')

Allmerge<-merge(MXP0,c(MXM6,MXY1,MXY2,MXYA))
DefaultAssay(Allmerge)<-'SCT'
Allmerge<-RunPCA(Allmerge,npcs = 100,features = rownames(Allmerge))
batch_info <- as.factor(Allmerge@meta.data$sample_Macaca)
table(batch_info)
pca_coord <- Allmerge@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
Allmerge@reductions$pca@cell.embeddings <- pca_corrected$corrected
Allmerge<-RunUMAP(Allmerge,dims = c(1:40))
Allmerge<-FindNeighbors(Allmerge,reduction = 'pca',dims=c(1:40))
Allmerge<-FindClusters(Allmerge,resolution = 5)
DimPlot(Allmerge,reduction = 'umap')

saveRDS(Allmerge,file='/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_allmerge_predSCT.rds')

###############subtypeUMAP
MacacaXenium<-readRDS(file="/home/gpfs/wulab15/MacacaPFC/Xenium/MacacaXenium_allmerge_0612.rds'")

MacacaXenium@meta.data$subtype<-factor(MacacaXenium@meta.data$subtype,levels=c('AST_FB', 'AST_PP','OPC','NOL','MOL','Microglia',
                                                                     'InN_PV', 'InN_SST','InN_VIP','InN_ID2',
                                                                     'ExN_L2_IT','ExN_L2/3_IT','ExN_L3/4_IT','ExN_L4_IT',
                                                                     'ExN_L4/5_IT', 'ExN_L5_IT', 'ExN_L5_PT','ExN_L5/6_NP', 'ExN_L5/6_IT', 
                                                                     'ExN_L6_IT','ExN_L6_CT', 'ExN_L6b',
                                                                     'Endothelial Cell'))
Idents(MacacaXenium)<-'subtype'
DimPlot(MacacaXenium,label = T,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312",
                                 "#044D2B","palegreen4","olivedrab","forestgreen",
                                 "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                 "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                 "#9C964A"))

Idents(MacacaXenium)<-'subtype'
p<-DimPlot(MacacaXenium,label = F,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312",
                                    "#044D2B","palegreen4","olivedrab","forestgreen",
                                    "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                    "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                    "#9C964A"),raster=F)
pdf("/home/gpfs/home/wulab15/Figures/UMAP/MacacaXenium_Subtype_UMAP.pdf",height = 4.5,width = 8)
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/UMAP/MacacaXenium_Subtype_UMAP.tiff",height = 450,width = 700)
p
dev.off()






