library(Seurat)

protein_coding <- read.csv(file='~/protein_coding.csv')
human_PFC <- readRDS(file='~/HumanPFC_withnewlabel_20230103.rds')
macaque_PFC <- readRDS(file='~/MacacaPFC_withnewlabel_20230103.rds')

#embryonic model

protein_coding <- read.csv(file='~/protein_coding.csv')
gene_use <- intersect(rownames(human_PFC@assays$RNA@data),protein_coding$HGNC)
gene_use <- intersect(rownames(macaque_PFC@assays$RNA@data),gene_use)
length(gene_use)
human_age <- data.frame(cell=colnames(human_PFC),age=as.character(human_PFC$sample_year))
table(human_age$age)
human_age$age <- gsub('PFC_GW','',human_age$age)
human_age$age <- as.numeric(human_age$age)*7
table(human_age$age)

library(glmnet)
celltype_list <- unique(human_PFC$sample_year)
cell_use = ''
celltype_list <- unique(human_PFC$sample_year)
set.seed(111)
for(i in 1:length(celltype_list)){
    cell <- colnames(human_PFC[,human_PFC$sample_year==celltype_list[i]])
    
        cell = cell[sample(length(cell),(length(cell)/3))]
        cell_use = c(cell_use,cell)
        
   
    
}
pos <- which(colnames(human_PFC)%in%cell_use)
table(human_PFC$sample_year[pos])

human_PFC = subset(human_PFC,cells = cell_use)
human_age <- data.frame(cell=colnames(human_PFC),age_year = as.character(human_PFC$sample_year),age=as.character(human_PFC$sample_year))
human_age$age <- gsub('PFC_GW','',human_age$age)
human_age$age <- as.numeric(human_age$age)*7
huamn_count_early <- as.matrix(human_PFC@assays$RNA@counts[gene_use,])

library(glmnet)

human_age$age_trans <- log2(human_age$age)
lambda <- cv.glmnet(x=t(huamn_count_early),y=human_age$age_trans)
fit <- glmnet(t(huamn_count_early),y=human_age$age_trans,lambda=lambda$lambda.1se,alpha=0.01,nlambda = 100,standardize = T)
cell_use_h <- colnames(huamn_count_early)
results <- predict.glmnet(fit, t(huamn_count_early), type="response")
results <- as.data.frame(results)
human_age$pred <- results$s0
human_age$pred_day <- 2^human_age$pred
human_age_mean_early = data.frame(age = unique(human_age_sub$sample),pred_mean = 0, real_mean = 0)
human_age_list <- c(unique(as.character(human_age_sub$sample)))

for(i in 1:length(human_age_list)){
    human_age_mean_early[human_age_mean_early$age==human_age_list[i],]$real_mean=mean(human_age_sub[human_age_sub$sample==human_age_list[i],]$age)
    human_age_mean_early[human_age_mean_early$age==human_age_list[i],]$pred_mean=mean(human_age_sub[human_age_sub$sample==human_age_list[i],]$pred_day)
    
}

library(ggplot2)

ggplot(human_age_mean_early,aes(x=real_mean,y=pred_mean))+geom_point()+geom_abline(slope = 1)+ylim(c(0,300))+xlim(c(0,300))
human_age_em = human_age_mean_early
Idents(macaque_PFC)<-'sample_year'
macaque_PFC <- subset(macaque_PFC,idents = c('PFC_E90','PFC_E110'))
celltype_list <- unique(macaque_PFC$sample_year)
cell_use = ''
celltype_list <- unique(macaque_PFC$sample_year)
for(i in 1:length(celltype_list)){
   
    cell <- colnames(macaque_PFC[,macaque_PFC$sample_year==celltype_list[i]])
     set.seed(12345)
        cell = cell[sample(length(cell),(length(cell)/3))]
        cell_use = c(cell_use,cell)
}

pos <- which(colnames(macaque_PFC)%in%cell_use)
macaque_count1 <- macaque_PFC@assays$RNA@data[gene_use,pos]
macaque_PFC_sub1 = subset(macaque_PFC,cells = colnames(macaque_count1))
macaque_age <- data.frame(cell=colnames(macaque_PFC_sub1),age=as.character(macaque_PFC_sub1$sample_year))
macaque_age$age <- gsub('PFC_E','',macaque_age$age)
macaque_age$age <- as.numeric(macaque_age$age)
macaque_age$age_trans <- log2(macaque_age$age)

results_M<-predict.glmnet(fit, t(macaque_count1), type="response")
results_M <- as.data.frame(results_M)
macaque_age$pred <- results_M$s0
macaque_age$pred_day <- 2^macaque_age$pred
macaque_age_list <- c(unique(as.character(macaque_PFC$sample_Macaca)))
macaque_age_mean <- data.frame(age=macaque_age_list,pred_mean=0,real_mean=0)
pos <- match(macaque_age$cell,colnames(macaque_PFC_sub1))
grep('TRUE',is.na(pos))
macaque_age$sample =as.character(macaque_PFC_sub1$sample_Macaca)[pos]
macaque_age_list <- c(unique(as.character(macaque_PFC_sub1$sample_Macaca)))
for(i in 1:length(macaque_age_list)){
    macaque_age_mean[macaque_age_mean$age==macaque_age_list[i],]$real_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$age)
    macaque_age_mean[macaque_age_mean$age==macaque_age_list[i],]$pred_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$pred_day)
    
}

#early model
DefaultAssay(human_PFC)<-'RNA'
gene_use <- intersect(rownames(human_PFC),protein_coding$HGNC)
DefaultAssay(macaque_PFC)<-'RNA'
gene_use <- intersect(rownames(macaque_PFC),gene_use)
human_age <- data.frame(cell=colnames(human_PFC),age_year = as.character(human_PFC$sample_year),age=as.character(human_PFC$sample_year))
human_age$age <- gsub('^PFC_M1$',2,human_age$age)
human_age$age <- gsub('^PFC_M2$',34,human_age$age)
human_age$age <- gsub('^PFC_M3$',86,human_age$age)
human_age$age <- gsub('^PFC_M4$',118,human_age$age)
human_age$age <- gsub('^PFC_M6$',179,human_age$age)
human_age$age <- gsub('^PFC_M10$',301,human_age$age)
human_age$age <- gsub('^PFC_Y1$',422,human_age$age)
human_age$age <- gsub('^PFC_Y1M6$',627,human_age$age)
human_age$age <- gsub('^PFC_Y2$',785,human_age$age)
human_age$age <- gsub('^PFC_Y3$',3*365,human_age$age)
human_age$age <- gsub('^PFC_Y4$',4*365,human_age$age)
human_age$age <- gsub('^PFC_Y4M11$',ceiling(((4+(11/12))*365)),human_age$age)
human_age$age <- gsub('^PFC_Y6$',6*365,human_age$age)
human_age$age <- gsub('^PFC_Y6M6$',ceiling(6.5*365),human_age$age)
human_age$age <- gsub('^PFC_Y8$',8*365,human_age$age)
human_age$age <- gsub('^PFC_Y9$',9*365,human_age$age)
human_age$age <- gsub('^PFC_Y10$',10*365,human_age$age)
human_age$age <- gsub('^PFC_Y12$',12*365,human_age$age)
human_age$age <- gsub('^PFC_Y13$',13*365,human_age$age)
human_age$age <- gsub('^PFC_Y14$',14*365,human_age$age)
human_age$age <- gsub('^PFC_Y15$',15*365,human_age$age)
human_age$age <- gsub('^PFC_Y16$',16*365,human_age$age)
human_age$age <- gsub('^PFC_Y17$',17*365,human_age$age)
human_age$age <- gsub('^PFC_Y20$',20*365,human_age$age)
human_age$age <- gsub('^PFC_Y25$',25*365,human_age$age)
human_age$age <- as.numeric(human_age$age)
human_age$age_day = human_age$age

# human M1-Y8 as referece predict macaque
pos <- grep('^PFC_M1$|^PFC_M2$|^PFC_M3$|^PFC_M4$|^PFC_M6$|^PFC_M10$|^PFC_Y1$|^PFC_Y1M6$|^PFC_Y2$|^PFC_Y3$|^PFC_Y4$|^PFC_Y4M11$|^PFC_Y6$|^PFC_Y6M6$|^PFC_Y8$',human_age$age_year)
cell_use = ''
sample(length(cell),(length(cell)/3))
sample(10,3)
cell_use = ''
celltype_list <- unique(human_PFC$sample_year[pos])
for(i in 1:length(celltype_list)){
    set.seed(111)
    cell <- colnames(human_PFC)[human_PFC$sample_year==celltype_list[i]]
    print(length(cell))
    cell1 = cell[sample(length(cell),(length(cell)/3))]
    print(length(cell1))
    cell_use = c(cell_use,cell1)
        
   
}
pos <- which(colnames(human_PFC)%in%cell_use)
huamn_count_early <- as.matrix(human_PFC@assays$RNA@counts[gene_use,pos])
pos <- which(human_age$cell%in%cell_use)
human_age_early <- human_age[pos,]
human_age_early$age_trans <- log2(human_age_early$age_day)
library(glmnet)
lambda_early <- cv.glmnet(x=t(huamn_count_early),y=human_age_early$age_trans)
fit_early <- glmnet(t(huamn_count_early),y=human_age_early$age_trans,lambda=lambda_early$lambda.1se,alpha=0.01,nlambda = 100,standardize = T)
cell_use_early_h <- colnames(huamn_count_early)

results_early <-predict.glmnet(fit_early, t(huamn_count_early), type="response")
results_early <- as.data.frame(results_early)
human_age_early$pred <- results_early$s0
human_age_early$pred_day <- 2^human_age_early$pred
human_age_early$pred_year =human_age_early$pred_day/365
human_age_list <- c(unique(as.character(human_PFC$sample_Human)))
human_age_mean_early <- data.frame(age=human_age_list,pred_mean=0,real_mean=0)
human_age_mean_early[human_age_mean_early$age==human_age_list[1],]
pos <- match(human_age_sub$cell,colnames(human_PFC))
human_age_sub$sample <- as.character(human_PFC$sample_Human[pos])

for(i in 1:length(human_age_list)){
    human_age_mean_early[human_age_mean_early$age==human_age_list[i],]$real_mean=mean(human_age_sub[human_age_sub$sample==human_age_list[i],]$age)
    human_age_mean_early[human_age_mean_early$age==human_age_list[i],]$pred_mean=mean(human_age_sub[human_age_sub$sample==human_age_list[i],]$pred_day)
    
}

library(ggplot2)
ggplot(human_age_mean_early,aes(x=real_mean,y=pred_mean))+geom_point()+geom_abline(slope = 1)+ylim(c(0,3000))+xlim(c(0,3000))



# predict macaque age
Idents(macaque_PFC)<-'sample_year'
macaque_PFC <- subset(macaque_PFC,idents = c('PFC_P0','PFC_M6','PFC_Y1'))
macaque_count <- macaque_PFC@assays$RNA@data[gene_use,]
macaque_age <- data.frame(cell=colnames(macaque_PFC),age=macaque_PFC$sample_year)
macaque_age$age <- gsub('PFC_P0',0,macaque_age$age)
macaque_age$age <- gsub('PFC_M6',0.5*365,macaque_age$age)
macaque_age$age <- gsub('PFC_Y1',1*365,macaque_age$age)

macaque_age$age <- as.numeric(macaque_age$age)
macaque_age$age_day <- macaque_age$age*365

results_early_m <-predict.glmnet(fit_early, t(macaque_count), type="response")
results_early_m <- as.data.frame(results_early_m)
pos <- match(rownames(results_early_m),macaque_age$cell)
macaque_age <- macaque_age[pos,]
macaque_age$pred <- results_early_m$s0
macaque_age$pred_day <- 2^macaque_age$pred
macaque_age$pred_year = macaque_age$pred_day/365
macaque_age_list <- c(unique(as.character(macaque_PFC$sample_Macaca)))
pos <- match(macaque_age$cell,colnames(macaque_PFC))
macaque_age$sample =as.character(macaque_PFC$sample_Macaca)[pos]
macaque_age_list <- c(unique(macaque_age$sample))

#macaque_age_list
macaque_age_mean_early <- data.frame(age=macaque_age_list,pred_mean=0,real_mean=0)
macaque_age_mean_early[macaque_age_mean_early$age==macaque_age_list[1],]

for(i in 1:length(macaque_age_list)){
    macaque_age_mean_early[macaque_age_mean_early$age==macaque_age_list[i],]$real_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$age)
    macaque_age_mean_early[macaque_age_mean_early$age==macaque_age_list[i],]$pred_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$pred_year)
    
}
human_age_mean_early$group = 'Human'
macaque_age_mean_early$group = 'Macaque'
age_mean = rbind(human_age_mean_early,macaque_age_mean_early)

ggplot(age_mean,aes(x=real_mean,y=pred_mean,group=group,color=group))+geom_point()+geom_abline(slope=1,lty='dashed')+#+xlim(c(0,240))+ylim(c(0,240))+
theme_bw()

#late model

pos <- grep('^PFC_M1$|^PFC_M2$|^PFC_M3$|^PFC_M4$|^PFC_M6$|^PFC_M10$|^PFC_Y1$|^PFC_Y1M6$|^PFC_Y2$|^PFC_Y3$|^PFC_Y4$|^PFC_Y4M11$|^PFC_Y6$|^PFC_Y6M6$|^PFC_Y8$',human_age$age_year)
cell_use = ''
celltype_list <- unique(human_PFC$sample_year[-pos])
for(i in 1:length(celltype_list)){ 
    cell <- colnames(human_PFC)[human_PFC$sample_year==celltype_list[i]]
    print(length(cell))
    set.seed(111)
    cell1 = cell[sample(length(cell),(length(cell)/3))]
    print(length(cell1))
    cell_use = c(cell_use,cell1)
        
   
}
pos <- which(colnames(human_PFC)%in%cell_use)

huamn_count_late <- as.matrix(human_PFC@assays$RNA@counts[gene_use,pos])
pos <- which(human_age$cell%in%cell_use)
human_age_late <- human_age[pos,]
human_age_late$age_trans <- log2(human_age_late$age_day)
library(glmnet)
lambda_late <- cv.glmnet(x=t(huamn_count_late),y=human_age_late$age_trans)
fit_late <- glmnet(t(huamn_count_late),y=human_age_late$age_trans,lambda=lambda_late$lambda.1se,alpha=0.01,nlambda = 100,standardize = T)
cell_use_late_h <- colnames(huamn_count_late)
results_late <-predict.glmnet(fit_late, t(huamn_count_late), type="response")
results_late <- as.data.frame(results_late)
human_age_late$pred <- results_late$s0
human_age_late$pred_day <- 2^human_age_late$pred
human_age_late$pred_year =human_age_late$pred_day/365
pos = match(human_age_late$cell,colnames(human_PFC))
human_age_late$sample = human_PFC$sample_Human[pos]


# predict macaque age

Idents(macaque_PFC)<-'sample_year'
macaque_PFC <- subset(macaque_PFC,idents = c('PFC_P0','PFC_M6','PFC_Y1'),invert=T)
macaque_count <- macaque_PFC@assays$RNA@data[gene_use,]
macaque_age <- data.frame(cell=colnames(macaque_PFC),age=macaque_PFC$sample_year)
macaque_age$age <- gsub('PFC_Y','',macaque_age$age)
macaque_age$age <- as.numeric(macaque_age$age)
macaque_age$age_day <- macaque_age$age*365
results_late <-predict.glmnet(fit_late, t(macaque_count), type="response")
results_late <- as.data.frame(results_late)
pos <- match(rownames(results_late),macaque_age$cell)
macaque_age <- macaque_age[pos,]
macaque_age$pred <- results_late$s0
macaque_age$pred_day <- 2^macaque_age$pred
macaque_age$pred_year = macaque_age$pred_day/365
macaque_age_list <- c(unique(as.character(macaque_PFC$sample_Macaca)))
pos <- match(macaque_age$cell,colnames(macaque_PFC))
macaque_age$sample =as.character(macaque_PFC$sample_Macaca)[pos]
macaque_age_list <- c(unique(macaque_age$sample))

#macaque_age_list
macaque_age_mean_late <- data.frame(age=macaque_age_list,pred_mean=0,real_mean=0)
macaque_age_mean_late[macaque_age_mean_late$age==macaque_age_list[1],]

for(i in 1:length(macaque_age_list)){
    macaque_age_mean_late[macaque_age_mean_late$age==macaque_age_list[i],]$real_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$age)
    macaque_age_mean_late[macaque_age_mean_late$age==macaque_age_list[i],]$pred_mean=mean(macaque_age[macaque_age$sample==macaque_age_list[i],]$pred_year)
    
}
macaque_age_mean_late = macaque_age_mean_late


