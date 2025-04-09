library(CellChat)
library(Seurat)
library(ggplot2)

cellchat = readRDS(file='~/HumanPastneuron_cellchat_20231219.rds')

df.net <- subsetCommunication(cellchat)
df.net$interaction_name <- df.net$interaction_name_2
df.net[df.net$pval==0,]$pval=0.0000000001
df.net$pval_log <- -log10(df.net$pval)
df.net$interaction_pair <- paste0(df.net$source,'-',df.net$target)

pos<-which(df.net$interaction_pair%in%c('AST-ExN_L2_IT','AST-ExN_L2/3_IT',"AST-ExN_L3/4_IT","AST-ExN_L4_IT",
                                        'AST-ExN_L4/5_IT','AST-ExN_L5_IT',"AST-ExN_L5_PT","AST-ExN_L5/6_IT",
                                        "AST-ExN_L5/6_NP","AST-ExN_L6_CT","AST-ExN_L6_IT","AST-ExN_L6b",
                                        "AST-InN_ID2","AST-InN_PV","AST-InN_SST","AST-InN_VIP"))


df.net_sub <- df.net[pos,]
#df.net_sub$interaction_pair <- factor(df.net_sub$interaction_pair,levels = c('nociceptor_satellite glia cell','C-LTMR_satellite glia cell'))
df.net_sub$interaction_pair <- factor(df.net_sub$interaction_pair,levels = c("AST-InN_PV","AST-InN_SST","AST-InN_VIP","AST-InN_ID2",
                                                                             'AST-ExN_L2_IT','AST-ExN_L2/3_IT',"AST-ExN_L3/4_IT","AST-ExN_L4_IT",
                                                                             'AST-ExN_L4/5_IT','AST-ExN_L5_IT',"AST-ExN_L5_PT",
                                                                             "AST-ExN_L5/6_NP","AST-ExN_L5/6_IT","AST-ExN_L6_IT","AST-ExN_L6_CT","AST-ExN_L6b"))

#subset LR pairs between AST and  IT_NP   ,  InN    , ExN_L6b ,  PT_CT 
pos <- grep('IT|NP',df.net_sub$target)
mech_pos = ''
for (i in 1:length(strsplit(df.net_sub[pos,]$receptor,'_'))){
    
    if(strsplit(df.net_sub[pos,]$receptor,'_')[[i]][1]%in%type_deg[type_deg$cluster=='IT_NP',]$gene
||strsplit(df.net_sub[pos,]$receptor,'_')[[i]][2]%in%type_deg[type_deg$cluster=='IT_NP',]$gene){
        print(i)
         mech_pos <- c(mech_pos,i)
    }
   
}

mech_pos <- mech_pos[-1]
mech_pos = as.numeric(mech_pos)
df.net_sub_mech = df.net_sub[pos[mech_pos],]
pos <- grep('PT|CT',df.net_sub$target)

mech_pos = ''
for (i in 1:length(strsplit(df.net_sub[pos,]$receptor,'_'))){
    
    if(strsplit(df.net_sub[pos,]$receptor,'_')[[i]][1]%in%type_deg[type_deg$cluster=='PT_CT',]$gene
||strsplit(df.net_sub[pos,]$receptor,'_')[[i]][2]%in%type_deg[type_deg$cluster=='PT_CT',]$gene){
        print(i)
         mech_pos <- c(mech_pos,i)
    }
   
}

mech_pos <- mech_pos[-1]
mech_pos = as.numeric(mech_pos)
df.net_sub_pt_ct = df.net_sub[pos[mech_pos],]
pos <- grep('L6b',df.net_sub$target)
mech_pos = ''
for (i in 1:length(strsplit(df.net_sub[pos,]$receptor,'_'))){
    
    if(strsplit(df.net_sub[pos,]$receptor,'_')[[i]][1]%in%type_deg[type_deg$cluster=='ExN_L6b',]$gene
||strsplit(df.net_sub[pos,]$receptor,'_')[[i]][2]%in%type_deg[type_deg$cluster=='ExN_L6b',]$gene){
        print(i)
         mech_pos <- c(mech_pos,i)
    }
   
}

mech_pos <- mech_pos[-1]
mech_pos = as.numeric(mech_pos)
df.net_sub_l6b = df.net_sub[pos[mech_pos],]
pos <- grep('InN',df.net_sub$target)

mech_pos = ''
for (i in 1:length(strsplit(df.net_sub[pos,]$receptor,'_'))){
    
    if(strsplit(df.net_sub[pos,]$receptor,'_')[[i]][1]%in%type_deg[type_deg$cluster=='InN',]$gene
||strsplit(df.net_sub[pos,]$receptor,'_')[[i]][2]%in%type_deg[type_deg$cluster=='InN',]$gene){
        print(i)
         mech_pos <- c(mech_pos,i)
    }
   
}

mech_pos <- mech_pos[-1]
mech_pos = as.numeric(mech_pos)
df.net_sub_InN = df.net_sub[pos[mech_pos],]
df.net_deg <- rbind(df.net_sub_mech,df.net_sub_pt_ct,df.net_sub_l6b,df.net_sub_InN)

#df.net_deg <-df.net_sub
interaction_pair = unique(df.net_deg$interaction_pair)
interaction_pair <- as.character(interaction_pair)
strsplit(interaction_pair[1],'-')[[1]][2]
receptor_list = df.net_deg$receptor
receptor_list <- strsplit(receptor_list,'_')
receptor_list <- as.data.frame(do.call(rbind,receptor_list))
receptor_list$V1 = paste0(receptor_list$V1,'_',seq(1:nrow(receptor_list)))
receptor_list$V2 = paste0(receptor_list$V2,'_',seq(1:nrow(receptor_list)))
receptor_list <- unique(c(receptor_list$V1,receptor_list$V2))

#load scRNA-seq data to compute degs between neuron subtypes 
HumanEP<-readRDS(file="/gpfs2/wulab15/share/PFC/Human/cellchat/HumanASTNeuron_20231219.rds")
HumanP = HumanEP
rm(HumanEP)
gc()

HumanP$deg_type = as.character(HumanP$Astsubtype)
Idents(HumanP) <-'deg_type'
DefaultAssay(HumanP) = 'RNA'
type_deg = FindAllMarkers(HumanP,only.pos = T)
Idents(HumanP)<-'temptype'

for(i in 1:1){
  tmp = df.net_deg[df.net_deg$interaction_pair==interaction_pair[i],]
  receptor_list = tmp$receptor
  receptor_list <- strsplit(receptor_list,'_')
  receptor_list <- as.data.frame(do.call(rbind,receptor_list))
  receptor_list$V1 = paste0(receptor_list$V1,'_',seq(1:nrow(receptor_list)))
  receptor_list$V2 = paste0(receptor_list$V2,'_',seq(1:nrow(receptor_list)))
  receptor_list <- unique(c(receptor_list$V1,receptor_list$V2))
  receptor_list <- strsplit(receptor_list,'_')
  receptor_list <- as.data.frame(do.call(rbind,receptor_list))
  tmp = tmp[unique(as.numeric(receptor_list$V2)),]
  tmp$avg_expr = 0
  for (j in 1:nrow(tmp)){
    avg_ligand = AverageExpression(HumanP,features = tmp$ligand[j],group.by='temptype')$RNA[,strsplit(interaction_pair[i],'-')[[1]][1]]
    avg_receptor = sum(AverageExpression(HumanP,features = c(toupper(receptor_list[receptor_list$V2==j,]$V1)),group.by='temptype')$RNA[,strsplit(interaction_pair[i],'-')[[1]][2]])
    avg_expr = log2(avg_ligand+avg_receptor+1)#
    tmp$avg_expr[j] = avg_expr
    
  }
}

for(i in 2:length(interaction_pair)){
  tmp1 = df.net_deg[df.net_deg$interaction_pair==interaction_pair[i],]
  receptor_list = tmp1$receptor
  receptor_list <- strsplit(receptor_list,'_')
  receptor_list <- as.data.frame(do.call(rbind,receptor_list))
  receptor_list$V1 = paste0(receptor_list$V1,'_',seq(1:nrow(receptor_list)))
  receptor_list$V2 = paste0(receptor_list$V2,'_',seq(1:nrow(receptor_list)))
  receptor_list <- unique(c(receptor_list$V1,receptor_list$V2))
  receptor_list <- strsplit(receptor_list,'_')
  receptor_list <- as.data.frame(do.call(rbind,receptor_list))
  tmp1 = tmp1[unique(as.numeric(receptor_list$V2)),]
  print(nrow(tmp1))
  if (nrow(tmp1)>0){
    tmp1$avg_expr = 0
    for (j in 1:nrow(tmp1)){
      avg_ligand = AverageExpression(HumanP,features = tmp1$ligand[j],group.by='temptype')$RNA[,strsplit(interaction_pair[i],'-')[[1]][1]]
      avg_receptor = sum(AverageExpression(HumanP,features = toupper(receptor_list[receptor_list$V2==j,]$V1),group.by='temptype')$RNA[,strsplit(interaction_pair[i],'-')[[1]][2]])
      avg_expr = log2(avg_ligand+avg_receptor+1)
      tmp1$avg_expr[j] = avg_expr
    }
    
    tmp = rbind(tmp,tmp1)
  }
  else{
    tmp = rbind(tmp,tmp1)
  }
}

tmp_deg=tmp
tmp_deg$map_name = paste0(tmp_deg$interaction_pair,'_',tmp_deg$interaction_name_2)
df.net_deg$map_name = paste0(df.net_deg$interaction_pair,'_',df.net_deg$interaction_name_2)
pos <- match(df.net_deg$map_name,tmp_deg$map_name)
df.net_deg$avg_expr = tmp_deg$avg_expr[pos]
write.csv(df.net_deg,file='~/interaction_deg_table.csv')

tmp <- data.frame(interaction = df.net_deg$interaction_pair,lr = df.net_deg$interaction_name_2,value=df.net_deg$avg_expr)
tmp_w = reshape(tmp, idvar = "interaction", timevar = "lr", direction = "wide")
tmp_w[is.na(tmp_w)]=0
rownames(tmp_w)<- tmp_w$interaction
tmp_w <- tmp_w[,-1]
tmp_wt <- t(tmp_w)

pca_mat_re <- irlba::prcomp_irlba(tmp_wt, n=7)$x
rownames(pca_mat_re) <- rownames(tmp_wt)
x <- as.matrix(pca_mat_re)

g_umap1 <- umap_tbl <- uwot::umap(pca_mat_re[,c(1:7)],n_neighbors = 50) %>%
{colnames(.) <- c('UMAP_1', 'UMAP_2'); .} %>%
  as_tibble(rownames='gene')

cl = kmeans(pca_mat_re, 30, nstart = 10)
table(cl$cluster)
g_umap1$kmean = as.character(cl$cluster)
g_umap1$cluster = as.character(cl$cluster)
g_umap1$gene = rownames(pca_mat_re)
colnames(g_umap1) = c('gene','UMAP_1','UMAP_2','cluster')

kmean_list = unique(g_umap1$cluster)
lr_module = g_umap1[g_umap1$cluster
                    ==kmean_list[1],]$gene
lr_module <- gsub('value.','',lr_module)
pos <- which(df.net_deg$interaction_name_2%in%lr_module)

df.net_sub_k3 = df.net_deg[pos,]
lr_module = as.data.frame(df.net_sub_k3 %>%
                            group_by(interaction_pair) %>%
                            summarise_at(vars(avg_expr), list(name = mean)))

lr_module$y=kmean_list[1]

for(i in 2:length(kmean_list)){
    
    k2 = g_umap1[g_umap1$cluster==kmean_list[i],]$gene
k2 <- gsub('value.','',k2)
pos <- which(df.net_deg$interaction_name_2%in%k2)
length(pos)
df.net_sub_k2 = df.net_deg[pos,]
k2_prop = as.data.frame(df.net_sub_k2 %>%
  group_by(interaction_pair) %>%
  summarise_at(vars(avg_expr), list(name = mean)))
k2_prop$y=kmean_list[i]
    #k2_prop$name = scale(k2_prop$name)
    
    lr_module = rbind(lr_module,k2_prop)
}

lr_module_p = g_umap1[g_umap1$cluster
                    ==kmean_list[1],]$gene
lr_module_p <- gsub('value.','',lr_module_p)
pos <- which(df.net_deg$interaction_name_2%in%lr_module_p)
length(pos)

df.net_sub_k3 = df.net_deg[pos,]
lr_module_p = as.data.frame(df.net_sub_k3 %>%
                            group_by(interaction_pair) %>%
                            summarise_at(vars(prob), list(name = mean)))

lr_module_p$y=kmean_list[1]

for(i in 2:length(kmean_list)){
  
  k2 = g_umap1[g_umap1$cluster==kmean_list[i],]$gene
  k2 <- gsub('value.','',k2)
  
  pos <- which(df.net_deg$interaction_name_2%in%k2)
  length(pos)
  
  df.net_sub_k2 = df.net_deg[pos,]
  
  k2_prop = as.data.frame(df.net_sub_k2 %>%
                            group_by(interaction_pair) %>%
                            summarise_at(vars(prob), list(name = mean)))
  
  k2_prop$y=kmean_list[i]
  #k2_prop$name = scale(k2_prop$name)
  
  lr_module_p = rbind(lr_module_p,k2_prop)
}

colnames(lr_module)<- c('interaction_pair','avg_expr','y')
lr_module$p = lr_module_p$name
lr_module$y <- gsub('c','LR_',lr_module$y)
normalize <- function(x){(x-min(x))/(max(x)-min(x))}
lr_module$p = normalize(lr_module_p$name)
lr_module$interaction_pair = as.character(lr_module$interaction_pair)
lr_module_s <- 
  lr_module %>%
  group_by(y) %>%
  mutate_at(vars(avg_expr,p),list(scaled=scale))
lr_module_s$p_scaled_norm = normalize(lr_module_s$p_scaled)
lr_module_s$avg_expr_scaled_norm = normalize(lr_module_s$avg_expr_scaled)
lr_module_s$y = as.character(lr_module_s$y)
lr_module_s$interaction_pair = factor(lr_module_s$interaction_pair,levels = c( 'AST-ExN_L2_IT','AST-ExN_L2/3_IT','AST-ExN_L3/4_IT','AST-ExN_L4_IT','AST-ExN_L4/5_IT','AST-ExN_L5_IT','AST-ExN_L5/6_IT','AST-ExN_L6_IT','AST-ExN_L5/6_NP','AST-ExN_L5_PT','AST-ExN_L6_CT','AST-ExN_L6b','AST-InN_ID2','AST-InN_VIP','AST-InN_PV','AST-InN_SST'))
lr_module_s$y = factor(lr_module_s$y,levels = rev(c('c1','c2','c3','c4','c5','c6','c7','c8')))

#pdf('~/LR_clustering.pdf',height = 6,width = 7)
ggplot(lr_module_s,aes(x=interaction_pair,y=y,color=p_scaled_norm))+geom_point(aes(size=avg_expr_scaled_norm))+theme(axis.text.x = element_text(angle=45,h=1))+
scale_color_gradientn(colours=c('#E3E5E7','#C9E0EC','#B1A6EB','#2D0C44'))+theme_bw()+theme(axis.text.x = element_text(angle=45,h=1))+scale_size(range=c(1,10))
#dev.off()

