library(Seurat)
library(WGCNA)
library(ggplot2)
library(dplyr)
library(igraph)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# loading dataset

## loading human data

human_PFC <- readRDS(file='/gpfs2/wulab15/share/PFC/Human/HumanPFC_withEmbryoData_20230428.rds')
human_PFC$type_year = paste0(as.character(human_PFC$subtype),'_',as.character(human_PFC$sample_year))
type_list = as.data.frame(table(human_PFC$type_year))
sub_cell = ''

for(i in 1:nrow(type_list)){
    pos <- which(human_PFC$type_year==type_list$Var1[i])
    if(length(pos)>3){
        pos1 = sample(length(pos),ceiling((length(pos)/3)))
        cell = colnames(human_PFC)[pos[pos1]]
        
    }else{
        cell = colnames(human_PFC)[pos]
        
    }
    sub_cell = c(sub_cell,cell)
    
}
human_PFC$window = as.character(human_PFC$sample_year)

pos <- grep('^PFC_GW',human_PFC$window)
human_PFC$window[pos]<-'Embryo'
pos <- grep('^PFC_M',human_PFC$window)
human_PFC$window[pos]<-'I'
pos <- grep('^PFC_Y1$|^PFC_Y1M6$|^PFC_Y2$|^PFC_Y3$',human_PFC$window)
human_PFC$window[pos]<-'II'
pos <- grep('^PFC_Y4$|^PFC_Y4M11$',human_PFC$window)
human_PFC$window[pos]<-'III'
pos <- grep('^PFC_Y6$|^PFC_Y6M6$|^PFC_Y8$|^PFC_Y9$|^PFC_Y10$',human_PFC$window)
human_PFC$window[pos]<-'IV'
pos <- grep('^PFC_Y12$|^PFC_Y13$|^PFC_Y14$',human_PFC$window)
human_PFC$window[pos]<-'V'
pos <- grep('^PFC_Y15$|^PFC_Y16$|^PFC_Y17$',human_PFC$window)
human_PFC$window[pos]<-'VI'
pos <- grep('^PFC_Y20$|^PFC_Y25$',human_PFC$window)
human_PFC$window[pos]<-'VII'

sub_cell = sub_cell[-1]
human_PFC = subset(human_PFC,cells = sub_cell)
Idents(human_PFC)<-'window'

human_PFC_win1 <- subset(human_PFC,idents = 'I')
DefaultAssay(human_PFC_win1)<-'RNA'
human_PFC_win1 <- FindVariableFeatures(human_PFC_win1,nfeatures = 10000)
datExpr <- t(as.matrix(human_PFC_win1@assays$RNA@data))
km = kmeans(t(datExpr),200)

## Add the celltype
human_PFC_win1$NEWlayer = as.character(human_PFC_win1$subtype)
pos <- grep('^ExN_L2_IT$',human_PFC_win1$subtype)
human_PFC_win1$NEWlayer[pos] <- 'Layer2'

pos <- grep('^ExN_L2/3_IT$',human_PFC_win1$subtype)
human_PFC_win1$NEWlayer[pos] <- 'Layer3'

pos <- grep('^ExN_L3/4_IT_COL5A2_PHLDB2$',human_PFC_win1$NEWlabel)
human_PFC_win1$NEWlayer[pos] <- 'Layer3'

pos <- grep('^ExN_L4_IT$|ExN_L3/4_IT',human_PFC_win1$subtype)
human_PFC_win1$NEWlayer[pos] <- 'Layer4'

pos <- grep('ExN_L3/4_IT_RORB_NTNG1',human_PFC_win1$newlabel)
human_PFC_win1$NEWlayer[pos] <- 'Layer4'

pos <- grep('^ExN_L4/5_IT$|ExN_L5_IT|ExN_L5_PT|ExN_L5/6_NP',human_PFC_win1$subtype)
human_PFC_win1$NEWlayer[pos] <- 'Layer5'

pos <- grep('^ExN_L5/6_IT$|ExN_L6_IT|ExN_L6_CT|ExN_L6b',human_PFC_win1$subtype)
human_PFC_win1$NEWlayer[pos] <- 'Layer6'

txt = data.frame(gene = colnames(datExpr),colors = as.character(km$cluster))
DefaultAssay(human_PFC_win1)<-'RNA'

for(i in 1:length(unique(txt$colors))){
    
    human_PFC_win1 = AddModuleScore(human_PFC_win1,features = list(txt[txt$colors==unique(txt$colors)[i],]$gene),name = paste0('c_',unique(txt$colors)[i]))
  
}

colors = paste0('c_',unique(txt$colors),'1')
DotPlot(human_PFC_win1,features =colors[1:length(colors)],group.by = 'NEWlayer',col.min = 0)+theme(axis.text.x = element_text(angle=45,h=1))

enrich = human_PFC_win1@meta.data[c(colors[1:length(colors)],'NEWlayer')]
txt$enrich  =''

for(i in 1:length(unique(txt$colors))){
    col_use = paste0('c_',unique(txt$colors)[i],'1')
  t= enrich %>%
  group_by(NEWlayer) %>%
  summarise_at(vars(col_use), list(name = mean))
  txt[txt$colors ==unique(txt$colors)[i],]$enrich = t$NEWlayer[which.max(t$name)]
 }
emb = txt

t_test_d = data.frame(disease = '',module = '',pvalue = 0,len1 = 0,len2 = 0,type1= '',type2= '')

for(i in 1:length(unique(emb$colors))){
    for(j in 1:1){
        overlap = length(intersect(emb[emb$colors==unique(emb$colors)[i],]$gene,d_set[[j]]))
        list1 = length(d_set[[j]])
        list2 = length(emb[emb$colors==unique(emb$colors)[i],]$gene)
          popsize = 10001
        pvalue = phyper(overlap-1,list2,(popsize-list2),list1,lower.tail = FALSE, log.p = FALSE)
        tmp = data.frame(disease = names(d_set)[[j]],
                                     module = paste0('emb_',unique(emb$colors)[i]),pvalue = pvalue,len1 = list1,len2 = list2,
                        type1 =names(d_set)[[j]],type2= unique(emb[emb$colors==unique(emb$colors)[i],]$enrich))
        t_test_d = rbind(t_test_d,tmp)
    }
    
}

t_test_d = t_test_d[-1,]
t_test_d = t_test_d[t_test_d$pvalue<0.05,]

for(i in 1:length(unique(emb$colors))){
    for(j in 2:2){
        overlap = length(intersect(emb[emb$colors==unique(emb$colors)[i],]$gene,d_set[[j]]))
        list1 = length(d_set[[j]])
        list2 = length(emb[emb$colors==unique(emb$colors)[i],]$gene)
          popsize = 10001
        pvalue = phyper(overlap-1,list2,(popsize-list2),list1,lower.tail = FALSE, log.p = FALSE)
        tmp = data.frame(disease = names(d_set)[[j]],
                                     module = paste0('emb_',unique(emb$colors)[i]),pvalue = pvalue,len1 = list1,len2 = list2,
                        type1 =names(d_set)[[j]],type2= unique(emb[emb$colors==unique(emb$colors)[i],]$enrich))
        t_test_d = rbind(t_test_d,tmp)
    }
    
}

t_test_d = t_test_d[-1,]
t_test_d = t_test_d[t_test_d$pvalue<0.05,]
emb = txt
t_test_d = data.frame(disease = '',module = '',pvalue = 0,len1 = 0,len2 = 0,type1= '',type2= '')

for(i in 1:length(unique(emb$colors))){
    for(j in 3:3){
        overlap = length(intersect(emb[emb$colors==unique(emb$colors)[i],]$gene,d_set[[j]]))
        list1 = length(d_set[[j]])
        list2 = length(emb[emb$colors==unique(emb$colors)[i],]$gene)
          popsize = 10001
        pvalue = phyper(overlap-1,list2,(popsize-list2),list1,lower.tail = FALSE, log.p = FALSE)
        tmp = data.frame(disease = names(d_set)[[j]],
                                     module = paste0('emb_',unique(emb$colors)[i]),pvalue = pvalue,len1 = list1,len2 = list2,
                        type1 =names(d_set)[[j]],type2= unique(emb[emb$colors==unique(emb$colors)[i],]$enrich))
        t_test_d = rbind(t_test_d,tmp)
    }
    
}

t_test_d = t_test_d[-1,]
t_test_d = t_test_d[t_test_d$pvalue<0.05,]
emb = txt
d_set[['anxiety']] = anxiety_gene$x
t_test_d = data.frame(disease = '',module = '',pvalue = 0,len1 = 0,len2 = 0,type1= '',type2= '')
t_test_d = data.frame(disease = '',module = '',pvalue = 0,len1 = 0,len2 = 0,type1= '',type2= '')

for(i in 1:length(unique(emb$colors))){
    for(j in 4:4){
        overlap = length(intersect(emb[emb$colors==unique(emb$colors)[i],]$gene,d_set[[j]]))
        list1 = length(d_set[[j]])
        list2 = length(emb[emb$colors==unique(emb$colors)[i],]$gene)
          popsize = 10001
        pvalue = phyper(overlap-1,list2,(popsize-list2),list1,lower.tail = FALSE, log.p = FALSE)
        tmp = data.frame(disease = names(d_set)[[j]],
                                     module = paste0('emb_',unique(emb$colors)[i]),pvalue = pvalue,len1 = list1,len2 = list2,
                        type1 =names(d_set)[[j]],type2= unique(emb[emb$colors==unique(emb$colors)[i],]$enrich))
        t_test_d = rbind(t_test_d,tmp)
    }
    
}

t_test_d = t_test_d[-1,]
t_test_d = t_test_d[t_test_d$pvalue<0.05,]
d_set[['depression']] = depression_gene$x
t_test_d = data.frame(disease = '',module = '',pvalue = 0,len1 = 0,len2 = 0,type1= '',type2= '')

for(i in 1:length(unique(emb$colors))){
    for(j in 5:5){
        overlap = length(intersect(emb[emb$colors==unique(emb$colors)[i],]$gene,d_set[[j]]))
        list1 = length(d_set[[j]])
        list2 = length(emb[emb$colors==unique(emb$colors)[i],]$gene)
          popsize = 10001
        pvalue = phyper(overlap-1,list2,(popsize-list2),list1,lower.tail = FALSE, log.p = FALSE)
        tmp = data.frame(disease = names(d_set)[[j]],
                                     module = paste0('emb_',unique(emb$colors)[i]),pvalue = pvalue,len1 = list1,len2 = list2,
                        type1 =names(d_set)[[j]],type2= unique(emb[emb$colors==unique(emb$colors)[i],]$enrich))
        t_test_d = rbind(t_test_d,tmp)
    }
    
}

t_test_d = t_test_d[-1,]
t_test_d = t_test_d[t_test_d$pvalue<0.05,]
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

hg19_genes = readRDS(file='~/hs_ge')
promoters<- promoters(hg19_genes, 5000, 5000)
peaks = rownames(atac@assays$celltypepeaks@data)
peaks <- strsplit(peaks,'-')
peaks = do.call(rbind,peaks)
peaks = as.data.frame(peaks)
colnames(peaks) <- c('chr','start','end')
peaks$strand = '*'
peaks = makeGRangesFromDataFrame(peaks)
hits <- findOverlaps(query = peaks, subject = promoters)

promoter.peaks = peaks[queryHits(hits)]
promoter.peaks_df = data.frame(chr = seqnames(promoter.peaks),start =granges((promoter.peaks))@ranges@start,end = granges((promoter.peaks))@ranges@start+granges((promoter.peaks))@ranges@width-1)
promoter.peaks_chr = paste0(promoter.peaks_df$chr,'-',promoter.peaks_df$start,'-',promoter.peaks_df$end)
DefaultAssay(atac) <-'celltypepeaks'

 match_gene = ClosestFeature(
    object = atac,
    regions = promoter.peaks_chr)

pos <- which(atac_g_deg$gene%in%match_gene$gene_name)
atac_g_deg = atac_g_deg[pos,]

pos <- which(atac_g_deg$gene%in%tf_family$Symbol)
atac_g_deg = atac_g_deg[pos,]

pos <- match(atac_g_deg$gene,tf_family$Symbol)
grep('TRUE',is.na(pos))

atac_g_deg$tf_family = ''
atac_g_deg$tf_family = tf_family$Family[pos]
avg_atac_h = AverageExpression(atac,features =  atac_g_deg$gene,group.by = 'NEWlayer',assays = 'RNA')
avg_atac_h = as.data.frame(avg_atac_h$RNA)
avg_atac_h$gene = rownames(avg_atac_h)
avg_atac_h_l = melt(avg_atac_h,id='gene')

pos <- match(avg_atac_h_l$gene,tf_family$Symbol)
grep('TRUE',is.na(pos))

avg_atac_h_l$family = ''
avg_atac_h_l$family =tf_family$Family[pos]
atac_g_deg$avg = ''

for(i in 1:nrow(atac_g_deg)){arr
    pos <- which(avg_atac_h_l$variable==atac_g_deg$cluster[i])
    pos1 <- which(avg_atac_h_l$gene[pos]==atac_g_deg$gene[i])
    atac_g_deg[i,]$avg = avg_atac_h_l[pos[pos1],]$value
    
}
    
atac_g_deg_w =data.frame(tf_family = atac_g_deg$tf_family,cluster = atac_g_deg$cluster,value=atac_g_deg$avg)
atac_g_deg_w$value = as.numeric(atac_g_deg_w$value)
tmp_h <- atac_g_deg_w %>% group_by(tf_family, cluster)%>% 
    summarize(mean=mean(value))
tmp_h$zscore = scale(tmp_h$mean)
tmp_h$log = log2(tmp_h$mean)

ggplot(tmp_h,aes(x=cluster,y=tf_family,fill=log))+geom_tile(color='gray80')+
theme(axis.text.x = element_text(angle=45,h=1))+scale_fill_gradientn(colors = c('gray90','lightblue','purple','black'))


tmp_h_atac_g_prom = tmp_h
DefaultAssay(atac) <-'celltypepeaks'

atac <- AddMotifs(
  object = atac,
  genome =BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)   

atac$NEWlayer = as.character(atac$subtype)
Idents(atac)<-'NEWlayer'
DefaultAssay(atac)<-'celltypepeaks'
atac.list <- SplitObject(atac,split.by = 'NEWlayer')
atac.list <- lapply(X = atac.list, FUN = function(x) {
x <- AddMotifs(
  object = x,
  genome =BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)    
x = RegionStats(
  object = x,
  assay = 'celltypepeaks',
  genome = BSgenome.Hsapiens.UCSC.hg19
)
    
})

tf_famile = data.frame(motif = rep('',length(pfm)),tf = rep('',length(pfm)),family = rep('',length(pfm)))
for (i in 1:length(pfm)){
    print(i)
    tf_famile$motif[i] = names(pfm@listData)[[i]]
    tf_famile$tf[i] = pfm@listData[[i]]@name
    if (length(pfm@listData[[i]]@tags$family)==1){
    tf_famile$family[i] = pfm@listData[[i]]@tags$family
    }
    else{
        
       tf_famile$family[i] = paste(combine(pfm@listData[[i]]@tags$family),collapse = '_')
    }
}

library(survcomp)
kmean_peak_cl_order  = kmean_peak_cl
da_peaks <- read.csv(file='/gpfs2/wulab10/adult_PFC/human_layer_da_peak.csv')
pos <- which(names(kmean_peak_cl_order)%in%da_peaks$gene)
kmean_peak_cl_order = kmean_peak_cl_order[pos]
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Layer2',]$gene))] = 'Layer2'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Layer3',]$gene))] = 'Layer3'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Layer4',]$gene))] = 'Layer4'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Layer5',]$gene))] = 'Layer5'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Layer6',]$gene))] = 'Layer6'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='InN_SST',]$gene))] = 'InN_SST'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='InN_VIP',]$gene))] = 'InN_VIP'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='InN_ID2',]$gene))] = 'InN_ID2'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='InN_PV',]$gene))] = 'InN_PV'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='AST_FB',]$gene))] = 'AST'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='AST_PP',]$gene))] = 'AST'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='AST',]$gene))] = 'AST'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Endothelial Cell',]$gene))] = 'Endothelial Cell'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='Microglia',]$gene))] = 'Microglia'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='NOL',]$gene))] = 'OLIG'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='MOL',]$gene))] = 'OLIG'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='OLIG',]$gene))] = 'OLIG'
kmean_peak_cl_order[(grep('TRUE',names(kmean_peak_cl_order)%in%da_peaks[da_peaks$cluster=='OPC',]$gene))] = 'OPC'

cre_module = data.frame(group = c(rep('Layer2',1),rep('Layer3',1),rep('Layer4',1),
                                 rep('Layer5',1),rep('Layer6',1),'InN_SST','InN_PV','InN_VIP','InN_ID2',
                                 'AST','Endothelial Cell','Microglia','OLIG','OPC'),
                        subtype = c('Layer2','Layer3','Layer4','Layer5','Layer6',
                                    'InN_SST','InN_PV','InN_VIP','InN_ID2',
                                    'AST','Endothelial Cell',
                                   'Microglia','OLIG','OPC'))

DefaultAssay(atac) <-'celltypepeaks'
link_peak <- data.frame(X=names(kmean_peak_cl_order),x=kmean_peak_cl_order)
meta.feature <- GetAssayData(atac.list[[i]], assay = "celltypepeaks", slot = "meta.features")

DefaultAssay(atac.list[[i]]) <-'celltypepeaks'
cre_g = cre_module[cre_module$subtype==names(atac.list)[[i]],]$group
tmp_peak = link_peak[link_peak$x==cre_module[cre_module$subtype==names(atac.list)[[i]],]$group,]$X
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[tmp_peak, ],
  n = 100
)

enriched.motifs <- FindMotifs(
  object = atac.list[[i]],
  features = tmp_peak,
    background=peaks.matched,
    p.adjust = T,
    p.adjust.methods='BH'
)

pos <- match(tf_famile$motif,enriched.motifs$motif)
enriched.motifs$tf=''
enriched.motifs$tf[pos] <- tf_famile$motif
enriched.motifs$family = ''
enriched.motifs$family[pos] <- tf_famile$family
names(atac.list)[[i]] =  gsub('\\/','_',names(atac.list)[[i]])
enriched.motifs$p_adj = p.adjust(enriched.motifs$pvalue,'BH')#
enriched.motifs$subtype = names(atac.list)[[i]] 
for (i in 1:1){
  
Idents(atac.list[[i]]) <-'NEWlayer'
open.peaks <- AccessiblePeaks(atac.list[[i]], idents = names(atac.list)[[i]],min.cells = 5)
meta.feature <- GetAssayData(atac.list[[i]], assay = "celltypepeaks", slot = "meta.features")

DefaultAssay(atac.list[[i]]) <-'celltypepeaks'
cre_g = cre_module[cre_module$subtype==names(atac.list)[[1]],]$group
tmp_peak = link_peak[link_peak$x==cre_module[cre_module$subtype==names(atac.list)[[1]],]$group,]$X
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[tmp_peak, ],
  n = 150
)
    
enriched.motifs <- FindMotifs(
  object = atac.list[[i]],
  features = tmp_peak,
    background=peaks.matched,
    p.adjust = T,
    p.adjust.methods='BH'
)

#enriched.motifs$p_adj <- p.adjust(enriched.motifs$pvalue,'BH') 
pos <- match(tf_famile$motif,enriched.motifs$motif)
enriched.motifs$tf=''
enriched.motifs$tf[pos] <- tf_famile$motif
enriched.motifs$family = ''
enriched.motifs$family[pos] <- tf_famile$family
names(atac.list)[[i]] =  gsub('\\/','_',names(atac.list)[[i]])
enriched.motifs$p_adj = p.adjust(enriched.motifs$pvalue,'BH')#
enriched.motifs$subtype = names(atac.list)[[1]]   
        }

enriched.motifs_all = enriched.motifs
for (i in 2:length(atac.list)){
  print(i)
    print(names(atac.list)[[i]])
Idents(atac.list[[i]]) <-'NEWlayer'
open.peaks <- AccessiblePeaks(atac.list[[i]], idents = names(atac.list)[[i]],min.cells = 5)
meta.feature <- GetAssayData(atac.list[[i]], assay = "celltypepeaks", slot = "meta.features")

DefaultAssay(atac.list[[i]]) <-'celltypepeaks'
cre_g = cre_module[cre_module$subtype==names(atac.list)[[i]],]$group
tmp_peak = link_peak[link_peak$x==cre_module[cre_module$subtype==names(atac.list)[[i]],]$group,]$X
bg_n = cre_module[cre_module$subtype==names(atac.list)[[i]],]$bg_n
    peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[tmp_peak, ],
  n = bg_n
)

enriched.motifs <- FindMotifs(
  object = atac.list[[i]],
  features = tmp_peak,
    background=peaks.matched,
    p.adjust = T,
    p.adjust.methods='BH',
)
pos <- match(tf_famile$motif,enriched.motifs$motif)
enriched.motifs$tf=''
enriched.motifs$tf[pos] <- tf_famile$motif
enriched.motifs$family = ''
enriched.motifs$family[pos] <- tf_famile$family
names(atac.list)[[i]] =  gsub('\\/','_',names(atac.list)[[i]])

enriched.motifs$p_adj = p.adjust(enriched.motifs$pvalue,'BH')#
enriched.motifs$subtype = names(atac.list)[[i]]   
enriched.motifs_all <- rbind(enriched.motifs_all,enriched.motifs)   
    print(names(atac.list)[i])
    
    }
pos <- grep('FALSE',is.na(enriched.motifs_all$p_adj))
enriched.motifs_all$p_adj_log = -log10(enriched.motifs_all$p_adj+1e-320)
enriched.motifs_all$p_adj_log_z =scale(enriched.motifs_all$p_adj_log)
enriched.motifs_all[enriched.motifs_all$p_adj>=0.05,]$p_adj=NA
enriched.motifs_all_w = data.frame(motif = enriched.motifs_all$motif.name,pvalue = enriched.motifs_all$p_adj_log,subtype = enriched.motifs_all$subtype)

library(reshape2)
enriched.motifs_all_w  = reshape(enriched.motifs_all_w, idvar = "motif", timevar = "subtype", direction = "wide")
rownames(enriched.motifs_all_w) = enriched.motifs_all_w$motif
enriched.motifs_all_w = enriched.motifs_all_w[,-1]
enriched.motifs_all_w <- as.matrix(enriched.motifs_all_w)
kmean_peak = kmeans(enriched.motifs_all_w,20)
kmean_peak_cl = kmean_peak$cluster
kmean_peak_cl_order = kmean_peak_cl[order(kmean_peak_cl,decreasing = F)]
row_order = unique(row_order)
row_order = c(7,8,14,16,10,13,17,3,20,12,6,1,2,18,15,4,5,9,11,19)
kmean_peak_cl_order = kmean_peak_cl_order[order(match(kmean_peak_cl_order,row_order))]
enriched.motifs_all = enriched.motifs_all[order(match(enriched.motifs_all$motif.name,names(kmean_peak_cl_order))),]
enriched.motifs_all$motif.name = factor(enriched.motifs_all$motif.name,levels=c(names(kmean_peak_cl_order)))
pos <-match(enriched.motifs_all$motif.name,names(kmean_peak_cl_order))
grep('TRUE',is.na(pos))
enriched.motifs_all$kmean = kmean_peak_cl_order
enriched.motifs_all_w = enriched.motifs_all_w[order(match(rownames(enriched.motifs_all_w),names(kmean_peak_cl_order))),]

library(ComplexHeatmap)
library(circlize)

ha_left = HeatmapAnnotation(kmean = kmean_peak_cl_order,
    col = list(kmean = c("1" = "red", "2" = "green", "3" = "blue",
                        "4"='#9A29D6','5'='#613C75',
                        '6'='#3245B7','7'='#66B732',
                        '8'='#5C8442','9'='#847F42',
                        '10'='goldenrod2','11'='#D844D1',
                        '12'='#6BE8E8','13'='#4C8E8E','14'='#45EDBD',
                        '15'='#19470F','16'='beige','17'='tan3','18'='#19470F',
                         '19'='#603611','20'='#AF5546')),which='row')

f1 <- circlize::colorRamp2(seq(0,10,length=5),c('gray90','lightblue','purple','#230444','black'))#used
enriched.motifs_all$subtype = factor(enriched.motifs_all$subtype,levels=c('Layer2','Layer3','Layer4','Layer5','Layer6','InN_PV', 'InN_SST','InN_VIP','InN_ID2',
    'AST','OPC','OLIG','Microglia','Endothelial Cell'))
colnames(enriched.motifs_all_w) <- gsub('pvalue.','',colnames(enriched.motifs_all_w))
pos <- match(c('Layer2','Layer3','Layer4','Layer5','Layer6','InN_PV', 'InN_SST','InN_VIP','InN_ID2',
    'AST','OPC','OLIG','Microglia','Endothelial Cell'),colnames(enriched.motifs_all_w))
pos <- match(c('Layer2','Layer3','Layer4','Layer5','Layer6','InN_PV', 'InN_SST','InN_VIP','InN_ID2',
    'AST','OPC','OLIG','Microglia','Endothelial Cell'),enriched.motifs_all$subtype)
ggplot(enriched.motifs_all[p,],aes(x=subtype,y=motif.name,fill=p_adj_log_l))+geom_tile(color='gray80')+
theme(axis.text.x = element_text(angle=45,h=1))+scale_fill_gradientn(colors = c('gray90','lightblue','purple','black'))

enriched.motifs_all$subtype = as.character(enriched.motifs_all$subtype)

