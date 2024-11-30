library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(BSgenome.Mfascicularis.NCBI.5.0)
library(AnnotationDbi)
library(GenomicRanges)

df <- read.table("/home/gpfs/home/wulab15/cellranger_atac_software/macaca_reference/Macaca_ref/fasta/genome.fa.fai")
df$isCircular <- FALSE
df[df$V1=="MT",]$isCircular <- TRUE
Macaca <- Seqinfo(seqnames= (as.character(df$V1))[1:21],
                  seqlengths=df$V2[1:21],
                  isCircular=df$isCircular[1:21],
                  genome="Macaca")

counts <- Read10X_h5(filename = "/home/gpfs/wulab15/MacacaPFC/scATAC/outs_aggr/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "/home/gpfs/wulab15/MacacaPFC/scATAC/outs_aggr/outs/singlecell.csv", header = TRUE, row.names = 1)
pos<- which(rownames((metadata))%in%colnames(counts))
length(pos)
metadata_filtered <- metadata[pos,]
chrom_assay <- CreateChromatinAssay(
  counts=counts,
  sep = c(":", "-"),
  genome=Macaca,
  fragments = '/home/gpfs/wulab15/MacacaPFC/scATAC/outs_aggr/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200)

obj =CreateSeuratObject(counts = chrom_assay,assay = "peaks",meta.data = metadata_filtered)
obj[['peaks']]

samples<-c('PFC_P0_1','PFC_P0_2','PFC_P0_3',
                             'PFC_M6_1','PFC_M6_2',
                             'PFC_Y1_1','PFC_Y1_2',
                             'PFC_Y2_1',
                             'PFC_Y4_1','PFC_Y4_2')

sample_name <- rownames(obj@meta.data)
obj@meta.data['sample_name']<-sample_name
obj@meta.data['sample_Macaca']<-NA
obj@meta.data['sample_age']<-NA
for(i in 1:length(samples)){
  samplei<-paste0('-',i,'$')
  posi<-grep(samplei,obj@meta.data$sample_name)
  obj@meta.data$sample_Macaca[posi]<-as.character(samples[i])
  obj@meta.data$sample_age[posi]<-as.character(paste0(strsplit(samples[i],'_')[[1]][1],'_',strsplit(samples[i],'_')[[1]][2]))}
table(obj@meta.data$sample_Macaca)
table(obj@meta.data$sample_age)

gtf <- read.table('/home/gpfs/home/wulab15/macaque_gtf_0925/macaque_gene_gtf_edited.txt',sep='\t',stringsAsFactors =F)
colnames(gtf)<-c('seqnames','start','end','strand','tx_id','gene_name','gene_biotype','type','exon_id')
gtf$exon_id[gtf$exon_id=='not_exon'] <- gtf$tx_id[gtf$exon_id=='not_exon']
gtf$type <- tolower(gtf$type)
pos <- grep('exon$|cds$|utr$|gap$',gtf$type)
gtf <- gtf[pos,]
gtf$type <- as.factor(gtf$type)
gtf$seqnames <- gsub('^Chr','',gtf$seqnames)
GR_gtf <- makeGRangesFromDataFrame(gtf,keep.extra.columns=T)
names(GR_gtf)<-GR_gtf$tx_id
GR_gtf

pos <- grep('^1$|^2$|^3$|^4$|^5$|^6$|^7$|^8$|^9$|^10$|^11$|^12$|^13$|^14$|^15$|^16$|^17$|^18$|^19$|^20$|^X$',seqnames(GR_gtf))
length(pos)
GR_gtf <- GR_gtf[pos,]
Annotation(obj) <- GR_gtf

obj <- NucleosomeSignal(object = obj)
obj <- TSSEnrichment(object = obj, fast = FALSE)
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments
obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(obj, group.by = 'high.tss') + NoLegend()
obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = obj, group.by = 'nucleosome_group')
VlnPlot(object = obj, features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal','nFeature_peaks'), pt.size = 0.1, ncol = 5)

obj@meta.data$sample_Macaca=factor(obj@meta.data$sample_Macaca,levels=c(samples))
Idents(obj)<-'sample_Macaca'
VlnPlot(object = obj, features = c('pct_reads_in_peaks'),pt.size = 0.1)
VlnPlot(object = obj, features = c('peak_region_fragments'),pt.size = 0.1)
VlnPlot(object = obj, features = c('TSS.enrichment'),pt.size = 0.1)
VlnPlot(object = obj, features = c('nucleosome_signal'),pt.size = 0.1)
VlnPlot(object = obj, features = c('nFeature_peaks'),pt.size = 0.1)

MacacaATAC <- subset(x = obj,
  subset = peak_region_fragments > 300 &
    peak_region_fragments < 35000 &
    pct_reads_in_peaks > 10 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nFeature_peaks  >1000)

MacacaATAC <- RunTFIDF(MacacaATAC)
MacacaATAC <- FindTopFeatures(MacacaATAC, min.cutoff = 'q0')
MacacaATAC <- RunSVD(MacacaATAC)
DepthCor(MacacaATAC)
MacacaATAC <- RunUMAP(object = MacacaATAC, reduction = 'lsi', dims = 2:25)
MacacaATAC <- FindNeighbors(object = MacacaATAC, reduction = 'lsi', dims = 2:25)
MacacaATAC <- FindClusters(object = MacacaATAC, verbose = FALSE, algorithm = 3,resolution=1)
DimPlot(object = MacacaATAC, label = TRUE) 

gene.activities <- GeneActivity(MacacaATAC)
MacacaATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
MacacaATAC <- NormalizeData(object = MacacaATAC,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(PFCaggrloose$nCount_RNA))

saveRDS(MacacaATAC, file = "/home/gpfs/wulab15/MacacaPFC/scATAC/MacacaATAC.rds")


###############subtypeUMAP
MacacaATAC<-readRDS(file="/home/gpfs/home/wulab15/MacacaPFC/202301/MacacaATAC_0103.rds")
MacacaATAC@meta.data$subtype<-factor(MacacaATAC@meta.data$subtype,levels=c('AST_FB', 'AST_PP','OPC','NOL','MOL','Microglia',
                                                                       'InN_PV', 'InN_SST','InN_VIP','InN_ID2',
                                                                       'ExN_L2_IT','ExN_L2/3_IT','ExN_L3/4_IT','ExN_L4_IT',
                                                                       'ExN_L4/5_IT', 'ExN_L5_IT', 'ExN_L5_PT','ExN_L5/6_NP', 'ExN_L5/6_IT', 
                                                                       'ExN_L6_IT','ExN_L6_CT', 'ExN_L6b',
                                                                       'Endothelial Cell'))
Idents(MacacaATAC)<-'subtype'
DimPlot(MacacaATAC,label = T,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312",
                                  "#044D2B","palegreen4","olivedrab","forestgreen",
                                  "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                  "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                  "#9C964A"))

Idents(MacacaATAC)<-'subtype'
p<-DimPlot(MacacaATAC,label = F,cols=c("#8E2043","#BC7A8F","#EBCC2A","#CAAB6C","#F2AD00","#C93312",
                                     "#044D2B","palegreen4","olivedrab","forestgreen",
                                     "#C6CDF7", "#8EA7DF", "#6986C4",  "#586BA4",  "#1F4B9D",  "#324376", 
                                     "#4B5695","#7979CE","#8079CE","#9F8DD3","#BCA5D9","#D2C6E5",
                                     "#9C964A"),raster=F)
pdf("/home/gpfs/home/wulab15/Figures/UMAP/MacacaATAC_Subtype_UMAP.pdf",height = 4.5,width = 8)
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/UMAP/MacacaATAC_Subtype_UMAP.tiff",height = 450,width = 700)
p
dev.off()

















