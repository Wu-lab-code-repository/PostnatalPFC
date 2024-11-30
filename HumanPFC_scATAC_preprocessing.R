library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)
library(GenomicRanges)
library(batchelor)

samples<-c('PFC_M1','PFC_M3','PFC_M6','PFC_M10','PFC_Y1',
           'PFC_Y2','PFC_Y4','PFC_Y6','PFC_Y8','PFC_Y10',
           'PFC_Y14','PFC_Y16','PFC_Y20','PFC_Y25')
for(i in 1:length(samples)){
  filepath<-paste0('/home/gpfs/home/wulab15/HumanPFC/scATAC/',samples[i],'/')
  setwd(filepath)
  load(file='/scATAC-cluster-call-peak-mat-binary.Robj')
  metadata <- read.csv(
    file = "snATAC_singlecell.csv",
    header = TRUE,
    row.names = 1)
  cell <- colnames(mat)
  samplesub<-paste0(samples[i],'#')
  cell <- gsub(samplesub,'',cell)
  pos <- which(rownames(metadata)%in%cell)
  metadata <- metadata[pos,]
  pos <- match(colnames(mat),rownames(metadata))
  grep("TRUE",is.na(pos))
  metadata <- metadata[pos,]
  chrom_assay <- CreateChromatinAssay(
    counts = mat,
    sep = c("-", "-"),
    genome = 'hg19',
    fragments = 'snATAC_fragments.tsv.gz',
    min.cells = 10,
    min.features = 200)
  obj =CreateSeuratObject(counts = chrom_assay,assay = "peaks",meta.data = metadata)
  obj[['peaks']]
  obj@meta.data['sample_name'] <- rownames(obj@meta.data)
  obj@meta.data['sample_Human']<-as.character(samples[i])
  obj@meta.data['sample_age']<-as.character(samples[i])
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg19"
  Annotation(obj) <- annotations
  obj <- NucleosomeSignal(object = obj)
  obj <- TSSEnrichment(object = obj, fast = FALSE)
  obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
  obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments
  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')
  filenamei<-paste0('HumanATAC_',samples[i],'_1207.rds')
  saveRDS(obj,file=filenamei)
}


obj_M1<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_M1_1207.rds")
obj_M3<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_M3_1207.rds")
obj_M6<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_M6_1207.rds")
obj_M10<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_M10_1207.rds")
obj_Y1<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y1_1207.rds")
obj_Y2<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y2_1207.rds")
obj_Y4<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y4_1207.rds")
obj_Y6<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y6_1207.rds")
obj_Y8<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y8_1207.rds")
obj_Y10<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y10_1207.rds")
obj_Y14<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y14_1207.rds")
obj_Y16<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y16_1207.rds")
obj_Y20<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y20_1207.rds")
obj_Y25<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_PFC_Y25_1207.rds")
combined.peaks <- UnifyPeaks(object.list = list(obj_M1, obj_M3,obj_M6,obj_M10,
                                                obj_Y1,obj_Y2,obj_Y4,obj_Y6,obj_Y8,obj_Y10,
                                                obj_Y14,obj_Y16,obj_Y20,obj_Y25), mode = "reduce")
combined.peaks
saveRDS(combined.peaks,file="/home/gpfs/home/wulab15/HumanPFC/scATAC/allcombinedpeaks.rds")
obj_M1.counts <- FeatureMatrix(
  fragments = Fragments(obj_M1),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_M1)
)
obj_M3.counts <- FeatureMatrix(
  fragments = Fragments(obj_M3),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_M3)
)
obj_M6.counts <- FeatureMatrix(
  fragments = Fragments(obj_M6),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_M6)
)
obj_M10.counts <- FeatureMatrix(
  fragments = Fragments(obj_M10),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_M10)
)
obj_Y1.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y1),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y1)
)
obj_Y2.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y2),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y2)
)
obj_Y4.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y4),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y4)
)
obj_Y6.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y6),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y6)
)
obj_Y8.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y8),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y8)
)
obj_Y10.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y10),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y10)
)
obj_Y14.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y14),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y14)
)
obj_Y16.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y16),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y16)
)
obj_Y20.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y20),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y20)
)
obj_Y25.counts <- FeatureMatrix(
  fragments = Fragments(obj_Y25),
  features = combined.peaks,
  sep = c("-", "-"),
  cells = colnames(obj_Y25)
)

obj_M1[['peaks']] <- CreateAssayObject(counts = obj_M1.counts)
obj_M3[['peaks']] <- CreateAssayObject(counts = obj_M3.counts)
obj_M6[['peaks']] <- CreateAssayObject(counts = obj_M6.counts)
obj_M10[['peaks']] <- CreateAssayObject(counts = obj_M10.counts)
obj_Y1[['peaks']] <- CreateAssayObject(counts = obj_Y1.counts)
obj_Y2[['peaks']] <- CreateAssayObject(counts = obj_Y2.counts)
obj_Y4[['peaks']] <- CreateAssayObject(counts = obj_Y4.counts)
obj_Y6[['peaks']] <- CreateAssayObject(counts = obj_Y6.counts)
obj_Y8[['peaks']] <- CreateAssayObject(counts = obj_Y8.counts)
obj_Y10[['peaks']] <- CreateAssayObject(counts = obj_Y10.counts)
obj_Y14[['peaks']] <- CreateAssayObject(counts = obj_Y14.counts)
obj_Y16[['peaks']] <- CreateAssayObject(counts = obj_Y16.counts)
obj_Y20[['peaks']] <- CreateAssayObject(counts = obj_Y20.counts)
obj_Y25[['peaks']] <- CreateAssayObject(counts = obj_Y25.counts)

combined_obj <- merge(x = obj_M1, y = list(obj_M3,obj_M6,obj_M10,
                                       obj_Y1,obj_Y2,obj_Y4,obj_Y6,obj_Y8,obj_Y10,
                                       obj_Y14,obj_Y16,obj_Y20,obj_Y25),
                  add.cell.ids = c("HumanM1","HumanM3","HumanM6", "HumanM10",
                                   "HumanY1","HumanY2","HumanY4","HumanY6","HumanY8","HumanY10",
                                   "HumanY14","HumanY16","HumanY20","HumanY25"))

VlnPlot(
  object = combined_obj,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal','nFeature_peaks','unmapped'),
  pt.size = 0.1,
  ncol = 6
)
combined_obj@meta.data$sample_Human=factor(combined_obj@meta.data$sample_Human,levels=samples)
Idents(combined_obj)<-'sample_Human'
VlnPlot(object = combined_obj, features = c('pct_reads_in_peaks'),pt.size = 0.1)
VlnPlot(object = combined_obj, features = c('peak_region_fragments'),pt.size = 0.1)
VlnPlot(object = combined_obj, features = c('TSS.enrichment'),pt.size = 0.1)
VlnPlot(object = combined_obj, features = c('nucleosome_signal'),pt.size = 0.1)
VlnPlot(object = combined_obj, features = c('nFeature_peaks'),pt.size = 0.1)
VlnPlot(object = combined_obj, features = c('unmapped'),pt.size = 0.1)

obj <- subset(
  x = combined_obj,
  subset = peak_region_fragments > 350 &
    peak_region_fragments < 35000 &
    pct_reads_in_peaks > 20 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_peaks > 4000 &
    unmapped < 5000
)

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj)
obj <- RunSVD(obj)

batch_info <- as.factor(obj@meta.data$sample_Human)
table(batch_info)
pca_coord <- obj@reductions$lsi@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
obj@reductions$lsi@cell.embeddings <- pca_corrected$corrected

obj <- RunUMAP(combined_obj, dims = 2:20, reduction = 'lsi')
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:20)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3,resolution = 1)

saveRDS(obj,file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_combinedall_1211.rds")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)
library(GenomicRanges)
library(batchelor)

HumanRNA<-readRDS(file='/home/gpfs/wulab15/HumanPFC/HumanRNA_EPmerge_dsforatac.rds')
HumanATAC<-readRDS(file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_combinedall_1211.rds")
HumanATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(HumanATAC) <- "ACTIVITY"
HumanATAC <- NormalizeData(
  object = HumanATAC,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(HumanATAC$nCount_peaks))
HumanATAC <- ScaleData(HumanATAC, features = rownames(HumanATAC))
transfer.anchors <- FindTransferAnchors(reference = HumanRNA, query = HumanATAC, features = rownames(HumanRNA),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca",k.filter = 5000)
genes.use <- rownames(HumanRNA)
refdata <- GetAssayData(HumanRNA, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = HumanATAC[["lsi"]],dims = 2:30)
HumanATAC[["predictedRNA"]] <- imputation
DefaultAssay(HumanATAC)<-'predictedRNA'
HumanATAC <- FindVariableFeatures(HumanATAC, selection.method = "vst", nfeatures = 3000)
HumanATAC <- ScaleData(HumanATAC,vars.to.regress = c("nFeature_peaks"))
HumanATAC <- RunPCA(HumanATAC, features =VariableFeatures(HumanATAC) ,npcs = 100) 
ElbowPlot(HumanATAC, ndims = 100, reduction = "pca")

batch_info <- as.factor(HumanATAC@meta.data$sample_Human)
table(batch_info)
pca_coord <- HumanATAC@reductions$pca@cell.embeddings
dim(pca_coord)
pca_corrected <- reducedMNN(pca_coord,batch = batch_info)
HumanATAC@reductions$pca@cell.embeddings <- pca_corrected$corrected

HumanATAC <- FindNeighbors(HumanATAC, dims = 2:22)
HumanATAC <- FindClusters(HumanATAC, resolution = 5)
HumanATAC <- RunUMAP(HumanATAC, dims = 2:22)
DimPlot(HumanATAC, reduction = "umap",label=TRUE)
Idents(HumanATAC)<-'sample_Human'
DimPlot(HumanATAC, reduction = "umap",label=TRUE)
DimPlot(HumanATAC, reduction = "umap",label=TRUE)+NoLegend()

saveRDS(HumanATAC,file="/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC.rds")


#########subtypeUMAP
HumanATAC<-readRDS(file='/home/gpfs/home/wulab15/HumanPFC/scATAC/HumanATAC_0109.rds')
HumanATAC@meta.data$subtype<-factor(HumanATAC@meta.data$subtype,levels=c('AST_FB', 'AST_PP','OPC','NOL','MOL','Microglia',
                                                                       'InN_PV', 'InN_SST','InN_VIP','InN_ID2',
                                                                       'ExN_L2_IT','ExN_L2/3_IT','ExN_L3/4_IT','ExN_L4_IT',
                                                                       'ExN_L4/5_IT', 'ExN_L5_IT', 'ExN_L5_PT','ExN_L5/6_NP', 'ExN_L5/6_IT', 
                                                                       'ExN_L6_IT','ExN_L6_CT', 'ExN_L6b',
                                                                       'Endothelial Cell'))
Idents(HumanATAC)<-'subtype'
p<-DimPlot(HumanATAC,label = F,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312", "#044D2B","palegreen4","olivedrab","forestgreen",
                                      "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376",
                                      "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                      "#9C964A"),raster=F)
p
pdf("/home/gpfs/home/wulab15/Figures/UMAP/HumanATAC_Subtype_UMAP_ATAC.pdf",height = 4.5,width = 7.5)
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/UMAP/HumanATAC_UMAP_ATAC.tiff",height = 450,width = 700)
p
dev.off()