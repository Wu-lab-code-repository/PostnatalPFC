library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(reshape2)
library(ggplot2)
library(SeuratDisk)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)

exnMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_ExN_subcluster.rds')
innMERGE<-readRDS(file='/home/gpfs/home/wulab15/species_merge/HMintegrate_InN_subcluster.rds')

DefaultAssay(exnMERGE)<-'RNA'
DefaultAssay(innMERGE)<-'RNA'
HMNmerge<-merge(exnMERGE,innMERGE)

DefaultAssay(HMNmerge)<-'RNA'
HMNmerge <- NormalizeData(HMNmerge, normalization.method = "LogNormalize", scale.factor = 10000)
dend1=readRDS(file="/home/gpfs/home/wulab15/species_merge/HMNeuron_dendrogram_0705.rds")

Idents(HMNmerge)<-'species'
table(HMNmerge@meta.data$species)
Idents(HMNmerge)<-'species'
Hmerge<-subset(HMNmerge,ident=c("Human"))
Mmerge<-subset(HMNmerge,ident=c("Macaca"))

Hmerge@meta.data$sample=factor(Hmerge@meta.data$sample,levels= c(paste0("H_",c("PFC_M1_1","PFC_M2_1","PFC_M3_1","PFC_M4_1","PFC_M6_1",PFC_M10_1",
                                                                                     "PFC_Y1_1","PFC_Y1_1","PFC_Y1M6_1","PFC_Y2_1","PFC_Y3_1",
                                                                                     "PFC_Y4_1","PFC_Y4_1","PFC_Y4M11_1",
                                                                                     "PFC_Y6_1","PFC_Y6M6_1","PFC_Y8_1","PFC_Y9_1","PFC_Y10#2","PFC_Y10_1",
                                                                                     "PFC_Y12_1","PFC_Y12_1","PFC_Y13_1","PFC_Y14_1","PFC_Y15_1",
                                                                                     "PFC_Y16_1","PFC_Y16_1","PFC_Y17_1","PFC_Y20_1","PFC_Y20_2","PFC_Y25_1"))))

Idents(Hmerge)<-'sample'
cluster.ids <- c("I","I","I","I","I","I",
                 "II","II","II","II","II",
                 "III","III","III",
                 "IV","IV","IV","IV","IV","IV",
                 "V","V","V","V",
                 "VI","VI","VI","VI",
                 "VII","VII","VII")
names(cluster.ids) <- c(paste0("H_",c("PFC_M1_1","PFC_M2_1","PFC_M3_1","PFC_M4_1","PFC_M6_1",PFC_M10_1",
                                                                                     "PFC_Y1_1","PFC_Y1_1","PFC_Y1M6_1","PFC_Y2_1","PFC_Y3_1",
                                                                                     "PFC_Y4_1","PFC_Y4_1","PFC_Y4M11_1",
                                                                                     "PFC_Y6_1","PFC_Y6M6_1","PFC_Y8_1","PFC_Y9_1","PFC_Y10#2","PFC_Y10_1",
                                                                                     "PFC_Y12_1","PFC_Y12_1","PFC_Y13_1","PFC_Y14_1","PFC_Y15_1",
                                                                                     "PFC_Y16_1","PFC_Y16_1","PFC_Y17_1","PFC_Y20_1","PFC_Y20_2","PFC_Y25_1")))

Hmerge <- RenameIdents(Hmerge, cluster.ids)
Hmerge@meta.data$window <- Idents(Hmerge)
Hmerge@meta.data$window <- factor(Hmerge@meta.data$window,levels = c("I","II","III","IV","V","VI","VII"))
Idents(Hmerge)<-'window'
DimPlot(Hmerge,cols = c("#4575B4","#74ADD1","#ABD9E9","#D4DCBA","#FEE08B","#FDC775","#FDAE61"),order = TRUE)

Idents(Hmerge)<-'sample'
table(Hmerge@meta.data$subcluster)
table(Hmerge@meta.data$sample)
HNeuron=as.data.frame.matrix(Hmerge@meta.data$subcluster%>%table(Hmerge@meta.data$sample))
HNeuron1=t(t(HNeuron)/colSums(HNeuron))
HNeuron1prop=HNeuron1/rowSums(HNeuron1)
rowSums(HNeuron1prop)
cluster_data=HNeuron1prop
colnames(cluster_data)<-colnames(HNeuron1prop)
rownames(cluster_data)<-rownames(HNeuron1prop)

I=(cluster_data[,1]+cluster_data[,2]+cluster_data[,3]+cluster_data[,4]+cluster_data[,5]+cluster_data[,6])/6
II=(cluster_data[,7]+cluster_data[,8]+cluster_data[,9]+cluster_data[,10]+cluster_data[,11])/5
III=(cluster_data[,12]+cluster_data[,13]+cluster_data[,14])/3
IV=(cluster_data[,15]+cluster_data[,16]+cluster_data[,17]+cluster_data[,18]+cluster_data[,19]+cluster_data[,20])/6
V=(cluster_data[,21]+cluster_data[,22]+cluster_data[,23]+cluster_data[,24])/4
VI=(cluster_data[,25]+cluster_data[,26]+cluster_data[,27]+cluster_data[,28])/4
VII=(cluster_data[,29]+cluster_data[,30]+cluster_data[,31])/3

cluster_data_win=cbind(I,II,III,IV,V,VI,VII)
samplesum=apply(cluster_data_win,1,sum)
df=cbind(cluster_data_win,samplesum)
df=data.frame(df)
df$Iprop=df$I/df$samplesum
df$IIprop=df$II/df$samplesum
df$IIIprop=df$III/df$samplesum
df$IVprop=df$IV/df$samplesum
df$Vprop=df$V/df$samplesum
df$VIprop=df$VI/df$samplesum
df$VIIprop=df$VII/df$samplesum

ctable=cbind(df$Iprop,df$IIprop,df$IIIprop,df$IVprop,df$Vprop,df$VIprop,df$VIIprop)
colnames(ctable)=c("I","II","III","IV","V","VI","VII")
rownames(ctable)<-c(rownames(cluster_data_win))
cluster_data=as.data.frame(ctable)
cluster_data
window=c()
clusternumber=length(rownames(cluster_data))
for (i in colnames(cluster_data)){
  window=c(window,rep(i,clusternumber))}
window

cluster_data<-as.matrix(cluster_data)
prop=c(cluster_data)
cluster=c(rep(rownames(cluster_data),length(colnames(cluster_data))))
proportion=as.numeric(prop)
clusterdf=data.frame(cbind(window,cluster,prop))
clusterdf
corder=c(labels(dend1))
clusterdforder=clusterdf
clusterdforder$cluster=factor(clusterdforder$cluster,levels = c(corder))
clusterdforder$window=factor(clusterdforder$window,levels = c("I","II","III","IV","V","VI","VII"))
clusterdforder$prop<-as.numeric(clusterdforder$prop)

clusterdforder
p<-ggplot(clusterdforder,mapping=aes(x=cluster,y=prop,fill=window)) +
  geom_bar(stat = 'identity',width = 0.9) +
  scale_fill_manual(values=c("#4575B4","#74ADD1","#ABD9E9","#D4DCBA","#FEE08B","#FDC775","#FDAE61"),limits=c("I","II","III","IV","V","VI","VII"))+   
  theme_bw()+ggtitle("Proportion(window)")+
  xlab("HumanPFC label")+ylab("Proportion")+
  theme(panel.grid =element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  theme(panel.border = element_blank())+theme(plot.title = element_text(size=10,face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10,face = "bold"),
                                              axis.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 10,face = "bold"),
                                              axis.title.x.bottom = element_text(angle = 0, hjust = 0, vjust = 0,size = 0),
                                              axis.title.y.left =  element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 0))+theme(plot.margin =margin(t = 1, r = 1, b = 1,l = 1,unit = "cm"))
p


pdf("/home/gpfs/home/wulab15/Figures/HumanNeuron_winpropbarplot.pdf",height = 6,width = 22)  
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/HumanNeuron_winpropbarplot.tiff",height = 400,width = 1500)
p
dev.off()

####macaca
Mmerge@meta.data$sample<-factor(Mmerge@meta.data$sample,levels = c('M_PFC_P0_1', 'M_PFC_P0_2', 'M_PFC_P0_3',
                                                                         'M_PFC_M6_1', 'M_PFC_M6_2',
                                                                         'M_PFC_Y1_1', 'M_PFC_Y1_2', 
                                                                         'M_PFC_Y2_1', 'M_PFC_Y2_2', 
                                                                         'M_PFC_Y4_1', 'M_PFC_Y4_2'))
table(Mmerge@meta.data$sample)
MNeuron=as.data.frame.matrix(Mmerge@meta.data$subcluster%>%table(Mmerge@meta.data$sample))
MNeuron1=t(t(MNeuron)/colSums(MNeuron))
MNeuron1prop=MNeuron1/rowSums(MNeuron1)
rowSums(MNeuron1prop)
cluster_data=MNeuron1prop
colnames(cluster_data)<-colnames(MNeuron1prop)
rownames(cluster_data)<-rownames(MNeuron1prop)
cluster_data
P0=(cluster_data[,1]+cluster_data[,2]+cluster_data[,3])/3
M6=(cluster_data[,4]+cluster_data[,5])/2
Y1=(cluster_data[,6]+cluster_data[,7])/2
Y2=(cluster_data[,8]+cluster_data[,9])/2
Y4=(cluster_data[,10]+cluster_data[,11])/2
cluster_data_win=cbind(P0,M6,Y1,Y2,Y4)

samplesum=apply(cluster_data_win,1,sum)
df=cbind(cluster_data_win,samplesum)
df=data.frame(df)
df$P0prop=df$P0/df$samplesum
df$M6prop=df$M6/df$samplesum
df$Y1prop=df$Y1/df$samplesum
df$Y2prop=df$Y2/df$samplesum
df$Y4prop=df$Y4/df$samplesum
ctable=cbind(df$P0prop,df$M6prop,df$Y1prop,df$Y2prop,df$Y4prop)
colnames(ctable)=c("P0","M6","Y1","Y2","Y4")
rownames(ctable)<-rownames(MNeuron1prop)
cluster_data=as.data.frame(ctable)
age=c()
clusternumber=length(rownames(cluster_data))
for (i in colnames(cluster_data)){
  age=c(age,rep(i,clusternumber))}
age

cluster_data<-as.matrix(cluster_data)
prop=c(cluster_data)
cluster=c(rep(rownames(cluster_data),length(colnames(cluster_data))))
proportion=as.numeric(prop)
clusterdf=data.frame(cbind(age,cluster,prop))
corder=c(labels(dend1))
clusterdforder=clusterdf
clusterdforder$cluster=factor(clusterdforder$cluster,levels = c(corder))
clusterdforder$age=factor(clusterdforder$age,levels = c("P0","M6","Y1","Y2","Y4"))
clusterdforder$prop<-as.numeric(clusterdforder$prop)
clusterdforder

p<-ggplot(clusterdforder,mapping=aes(x=cluster,y=prop,fill=age)) +
  geom_bar(stat = 'identity',width = 0.9) +
  scale_fill_manual(values=c("#4575B4","#74ADD1","#ABD9E9","#FEE08B","#FDAE61"),limits=c("P0","M6","Y1","Y2","Y4"))+  
  theme_bw()+ggtitle("Proportion(age)")+
  xlab("MacacaPFC label")+ylab("Proportion")+
  theme(panel.grid =element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  theme(panel.border = element_blank())+theme(plot.title = element_text(size=10,face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10,face = "bold"),
                                              axis.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 10,face = "bold"),
                                              axis.title.x.bottom = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 0),
                                              axis.title.y.left =  element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 0))+theme(plot.margin =margin(t = 1, r = 1, b = 1,l = 1,unit = "cm"))
p
pdf("/home/gpfs/home/wulab15/Figures/MacacaNeuron_agepropbarplot.pdf",height = 6,width = 22)   
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/MacacaNeuron_agepropbarplot.tiff",height = 400,width = 1500)
p
dev.off()

###
HMNmerge@meta.data$species<-factor(HMNmerge@meta.data$species,levels = c('Human',"Macaca"))
table(HMNmerge@meta.data$species)
MNeuron=as.data.frame.matrix(HMNmerge@meta.data$subcluster%>%table(HMNmerge@meta.data$species))
MNeuron1=t(t(MNeuron)/colSums(MNeuron))

MNeuron1prop=MNeuron1/rowSums(MNeuron1)
rowSums(MNeuron1prop)
cluster_data=MNeuron1prop
colnames(cluster_data)<-colnames(MNeuron1prop)
rownames(cluster_data)<-rownames(MNeuron1prop)
ctable<-cluster_data
colnames(ctable)=c('Human',"Macaca")
rownames(ctable)<-rownames(MNeuron1prop)
cluster_data=as.data.frame(ctable)
cluster_data
species=c()
clusternumber=length(rownames(cluster_data))
for (i in colnames(cluster_data)){
  species=c(age,rep(i,clusternumber))}
species

cluster_data<-as.matrix(cluster_data)
prop=c(cluster_data)
cluster=c(rep(rownames(cluster_data),length(colnames(cluster_data))))
proportion=as.numeric(prop)
clusterdf=data.frame(cbind(age,cluster,prop))
clusterdf

corder=c(labels(dend1))
clusterdforder$cluster=factor(clusterdforder$cluster,levels = c(corder))
clusterdforder$species=factor(clusterdforder$species,levels = c('Human',"Macaca"))
clusterdforder$prop<-as.numeric(clusterdforder$prop)
clusterdforder

p<-ggplot(clusterdforder,mapping=aes(x=cluster,y=prop,fill=species)) +
  geom_bar(stat = 'identity',width = 0.9) +
  scale_fill_manual(values=c("#3D578E","plum3"),limits=c('Human',"Macaca"))+ 
  theme_bw()+ggtitle("Proportion(species)")+
  xlab("label")+ylab("Proportion")+
  theme(panel.grid =element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  theme(panel.border = element_blank())+theme(plot.title = element_text(size=10,face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10,face = "bold"),
                                              axis.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 10,face = "bold"),
                                              axis.title.x.bottom = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 0),
                                              axis.title.y.left =  element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 0))+theme(plot.margin =margin(t = 1, r = 1, b = 1,l = 1,unit = "cm"))

p

pdf("/home/gpfs/home/wulab15/Figures/HMmergeNeuron_speciespropbarplot.pdf",height = 6,width = 22)   
p
dev.off()
tiff("/home/gpfs/home/wulab15/Figures/HMmergeNeuron_speciespropbarplot.tiff",height = 400,width = 1500)
p
dev.off()
