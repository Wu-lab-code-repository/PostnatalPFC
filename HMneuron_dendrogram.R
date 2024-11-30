library(Seurat)
library(ggplot2)
library(dendextend)
library(dplyr)
library(reshape2)

exnMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_ExN_subcluster.rds')
innMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_InN_subcluster.rds')

DefaultAssay(exnMERGE)<-'RNA'
DefaultAssay(innMERGE)<-'RNA'
HMNmerge<-merge(exnMERGE,innMERGE)

genes<-c(rownames(HMNmerge@assays$integrated))
Idents(HMNmerge)<-'subcluster'
data.avg <- AverageExpression(object = HMNmerge,assays = 'integrated',features =c(genes) )
data.dist <- dist(x = t(x = data.avg$integrated[genes, ])) 
dend <- data.dist %>%hclust() %>%  as.dendrogram()  
plot(dend)
par(mar=c(17,0.1,0.1,0.1))
reorder.dend <- function(dend, l.rank,verbose=FALSE)
{
  tmp.dend = dend
  sc=sapply(1:length(dend), function(i){
    l = dend[[i]] %>% labels
    mean(l.rank[dend[[i]] %>% labels])
  })
  ord = order(sc)
  if(verbose){
    print(sc)
    print(ord)
  }
  if(length(dend)>1){
    for(i in 1:length(dend)){
      if(ord[i]!=i){
        dend[[i]]= tmp.dend[[ord[i]]]
      }
      if(length(dend[[i]])>1){
        dend[[i]]=reorder.dend(dend[[i]],l.rank)
      }
    }
  }
  return(dend)
}


dendro_data = as.ggdend(dend)
alllabel=as.character(dendro_data$labels$label)

dendlabel=c("ExN_L2_IT_CUX2_CHN2","ExN_L2_IT_CUX2_PRSS12","ExN_L2_IT_PDGFD_NELL1","ExN_L2/3_IT_COL5A2_CBLN2",
            "ExN_L2/3_IT_CUX2_CHST9",
            "ExN_L2/3_IT_PDGFD_GLP2R","ExN_L2/3_IT_COL5A2_FTH1",
            "ExN_L2/3_IT_COL5A2_ADAMTS12",
            "ExN_L2/3_IT_COL5A2_TMEM163",
            "ExN_L2/3_IT_COL5A2_PDE1C","ExN_L2/3_IT_PDGFD_YBX1","ExN_L2/3_IT_COL5A2_VSNL1",
            "ExN_L3/4_IT_RORB_NTNG1","ExN_L3/4_IT_RORB_ARHGAP15","ExN_L3/4_IT_RORB_COL22A1",
            "ExN_L4_IT_RORB_LRRC16A","ExN_L4_IT_RORB_EPHA3",
            "ExN_L4/5_IT_SULF2_SLC38A11_SHISA6","ExN_L5_IT_ETV1_MKX","ExN_L5_IT_ETV1_BMP6","ExN_L5_IT_ETV1_PLD5",
            "ExN_L5_IT_ETV1_EPHB1","ExN_L5_IT_ETV1_PVRL3","ExN_L5_IT_ETV1_TMSB10",
            "ExN_L5_IT_ETV1_PCDH15","ExN_L5_IT_ETV1_THSD7B","ExN_L5_PT_SULF2_PCSK6_SCN4B",
            "ExN_L5/6_IT_OPRK1_DCSTAMP","ExN_L4/5_IT_SULF2_SLC38A11_INPP4B",
            "ExN_L6_IT_OPRK1_TESPA1_SNCG","ExN_L6_IT_OPRK1_TESPA1_PCSK1",
            "ExN_L6_IT_OPRK1_PDE3A","ExN_L6_IT_OPRK1_TESPA1_KCTD8","ExN_L6_IT_OPRK1_ABI3BP",
            "ExN_L5/6_NP_ETV1_HTR2C_IL26","ExN_L5/6_NP_ETV1_HTR2C_TNS3",
            "ExN_L6_CT_TLE4_ABO_NXPH1","ExN_L6_CT_TLE4_ABO_GPR155","ExN_L6_CT_TLE4_PDE8B",
            "ExN_L6b_TLE4_SCUBE1","ExN_L6b_TLE4_DPYSL3","ExN_L6b_TLE4_CELSR1","ExN_L6b_TLE4_NPFFR2",            
            
            "InN_PVALB_SCUBE3_MDGA2","InN_PVALB_SCUBE3_SOX4",
            "InN_PVALB_MEPE","InN_PVALB_SULF1_OSTN","InN_PVALB_WFDC2_IQGAP2","InN_PVALB_GABRG1",
            "InN_PVALB_SULF1_AGBL1","InN_PVALB_WFDC2_SORCS1",
            "InN_SST_CALB1_DPP10","InN_SST_STK32A","InN_SST_ADGRG6","InN_SST_CALB1_OPRM1",
            "InN_SST_CALB1_PTPRU","InN_SST_THSD7B","InN_SST_ITGA8","InN_SST_KCNH8","InN_SST_TH",
            "InN_SST_TXK","InN_SST_BTBD11","InN_SST_HPGD",
            "InN_SST_NPY","InN_SST_CALB1_CNTN6","InN_SST_EYA2",
            "InN_SST_FRZB","InN_SST_InN_SST_PRKCQ","InN_SST_MME",
            "InN_VIP_CHRM2","InN_VIP_SERPINF1","InN_VIP_CBLN1_GRIK3","InN_VIP_CBLN1_RGS6",
            "InN_VIP_HS3ST3A1","InN_VIP_LBH","InN_VIP_PDE1A","InN_VIP_PENK",
            "InN_VIP_OPRM1","InN_VIP_KCNIP1","InN_VIP_SPAG17","InN_VIP_ADAM12","InN_VIP_HTR1F",
            "InN_VIP_TYR",
            "InN_ID2_HTR7","InN_ID2_SYT10","InN_ID2_MC4R","InN_ID2_PAX6","InN_ID2_SYT6",
            "InN_ID2_NMBR",
            "InN_ID2_DPY19L1","InN_ID2_SCGN","InN_ID2_LHX6",
            "InN_ID2_LPHN2","InN_ID2_LCP2")

HMNmerge@meta.data['subcluster_dend']<-HMNmerge@meta.data$subcluster
HMNmerge@meta.data$subcluster_dend=factor(HMNmerge@meta.data$subcluster_dend,levels = c(dendlabel))

Idents(HMNmerge)<-'subcluster_dend'
cluster.ids <- c(1:94)
names(cluster.ids) <- c(names(table(HMNmerge@meta.data$subcluster_dend)))

HMNmerge <- RenameIdents(HMNmerge, cluster.ids)
HMNmerge@meta.data$lrank <- Idents(HMNmerge)
HMNmerge@meta.data$lrank=factor(HMNmerge@meta.data$lrank,levels=c(as.character(c(1:94))))
head(HMNmerge@meta.data)

l.rank =setNames(as.integer(HMNmerge@meta.data$lrank), HMNmerge@meta.data$subcluster_dend)
dend1 <- reorder.dend(dend,l.rank)
plot(dend1)

pdf("/home/gpfs/home/wulab15/Figures/HMNeuron_dendrogram_0705.pdf",height = 6,width = 30)
par(mar=c(17,0.1,0.1,0.1))
plot(dend1)
dev.off()

saveRDS(dend1,'/home/gpfs/home/wulab15/species_merge/HMNeuron_dendrogram_0705.rds')