

library(Seurat)
library(data.table)
library("sscVis")
library("sscClust")
library("scPip")
library(ggplot2)
library(ggridges)
library(dplyr)
library(cowplot)
library(PCAtools)
library(ggpubr)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")


sub_sce<-readRDS("3.Cluster/13.Annotation/3.PC_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("PN" ,   "PT", "mBrain", "mLN" ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("PN" ,   "PT", "mBrain", "mLN" ))
#dir.create("3.Cluster/13.plot/3.myeloid/")
sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
sub_sce$SubCluster<-gsub("PVL","PC",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(      "PC_C1_mPC-ACTG2" ,
                                                             "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif"))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c(     "mPC-ACTG2" ,
                       "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tSNE")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-unique(sort(sce$SubCluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"),
               RColorBrewer::brewer.pal(8,"Set1"))
names(g.colSet1)<-Cluster
#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(SubCluster1) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )


umap_p1<-ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c(   "PC_C1_mPC-ACTG2" ,
                                  "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif" ),
                     values = g.colSet1)+
  ggtitle("Clustering of PCs")+
  geom_text(aes(label = SubCluster1), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 9))+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=1))


umap_p1


write.table(umap_p1$data[,c("SubCluster","umap_1","umap_2")],"source_data_NC/Figure_4A.txt",sep = "\t",row.names = F,col.names = T)



#Merged plot
sub_sce<-readRDS("3.Cluster/13.Annotation/3.PC_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("PN" ,   "PT", "mBrain", "mLN"  ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("PN" ,   "PT", "mBrain", "mLN" ))
#dir.create("3.Cluster/13.plot/3.myeloid/")
sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
sub_sce$SubCluster<-gsub("PVL","PC",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(      "PC_C1_mPC-ACTG2" ,
                                                             "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif"))       


Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c(  "mPC-ACTG2" ,
                    "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tSNE")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-unique(sort(sce$SubCluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"),
               RColorBrewer::brewer.pal(8,"Set1"))
names(g.colSet1)<-Cluster
#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(SubCluster1) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )


umap_p2<-ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c(  "PC_C1_mPC-ACTG2" ,
                                 "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif" ),
                     values = g.colSet1)+
  ggtitle("Clustering of stromal cells")+
  #geom_text(aes(label = SubCluster1), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10),
        legend.position = "none")+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=2))+
  facet_wrap(~Type,ncol=2)

umap_p2





#VlnPlot
sub_sce<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
sub_sce$SubCluster<-gsub("PVL","PC",sub_sce$SubCluster)
    

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(      "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                                                             "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,    
                                                             "PC_C1_mPC-ACTG2" ,
                                                             "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif" ))       

VlnP1<-VlnPlot(sub_sce,group.by = "SubCluster",pt.size = 0,stack = T,flip = T,
               features = c("COL1A1","FAP","ACTA2","PI16","APOD","RGS5","MYH11"))



##MIlo
library(miloR)
milo.obj<-readRDS("4.characteristics/7.3.Milo_PC.milo.obj.RDS")
milo.res<-readRDS("4.characteristics/7.3.Milo_PC.milo.res.RDS")
nh_graph_endo<- plotNhoodGraphDA(milo.obj, milo.res, layout="umap",alpha=0.1) +
  theme_classic()+
  labs(title = "Nhoods of stromal cells (tumor vs. normal)",x="umap_1",y="umap_2")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 9))

milo.res$SubCluster<-as.character(milo.res$SubCluster)
library(limma)
milo.res$SubCluster<-strsplit2(milo.res$SubCluster,"-")[,1]
order_tab<-dplyr::summarize(group_by(milo.res,SubCluster),FC=mean(logFC))
order_tab<-order_tab[order(order_tab$FC),]
milo.res$SubCluster<-factor(milo.res$SubCluster,levels=order_tab$SubCluster)

beeswarm_endo<-plotDAbeeswarm(milo.res, group.by = "SubCluster")+
  theme_bw()+
  labs(title = "Stromal cells (tumor vs. normal)")+
  geom_hline(yintercept = c(0),linetype=c("solid"),size=0.5,color="grey")+
  theme(axis.text.y = element_text(size = 9), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))

milo.res$SubCluster<-gsub("PVL","PC",milo.res$SubCluster)
write.table(milo.res,"source_data_NC/Figure_4C-D.txt",sep = "\t",row.names = F,col.names = T)






sub_sce<-readRDS("3.Cluster/13.Annotation/3.PC_annotation.rds")
sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
sub_sce$SubCluster<-gsub("PVL","PC",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(      "PC_C1_mPC-ACTG2" ,
                                                             "PC_C2_mPC-MYH11" , "PC_C3_imPC-CD36"  ,"PC_C4_imPC-MCAM", "PC_C5-prolif"))       


Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c(       "mPC-ACTG2" ,
                     "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif" ) 


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
dot_p<-DotPlot(sub_sce, features = c("ACTG2","DES","MYLK",
                                          "MYH11","ADIRF","ACTA2",
                                          "RGS5","CCL19","CD36","COL4A1","GJA4",
                                          "MCAM","THY1","MKI67","STMN1"),
               dot.scale=5)+RotatedAxis()+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(5,"BuPu"))+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(legend.position = "none")+
  ggtitle("Markers of PCs")+
  #scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 9))

dot_p


write.table(dot_p$data,"source_data_NC/Figure_4B.txt",sep = "\t",row.names = F,col.names = T)



##GO_BP CD8
Enrich_type<-"Hallmark"
Type="Stromal"
total_table<-data.frame()

cluster="PT"
path<-paste("4.characteristics/6.Enrichment/3.PVL/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
if (file.exists(path)){
  enrich_table<-data.table::fread(path)
  enrich_table$Cluster<-cluster
  
  if (nrow(enrich_table)>6){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,6)
  }
  enrich_table$Type<-"Up"
  total_table<-rbind(total_table,enrich_table)
}


cluster="PN"
path<-paste("4.characteristics/6.Enrichment/3.PVL/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
if (file.exists(path)){
  enrich_table<-data.table::fread(path)
  enrich_table$Cluster<-cluster
  enrich_table<-enrich_table[-5,]
  if (nrow(enrich_table)>5){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,5)
  }
  enrich_table$Type<-"Down"
  enrich_table$Count<- -(enrich_table$Count)
  total_table<-rbind(total_table,enrich_table)
}




total_table$Description<-gsub("HALLMARK_","",total_table$Description)
total_table$Description<-factor(total_table$Description,levels=rev(unique(total_table$Description)))
#total_table$Cluster<-factor(total_table$Cluster,levels=Select_cluster_list)

total_table$text_y<- 0
total_table$pos_y<- 1
total_table$pos_y[which(total_table$Type=="Down")]<- 0
total_table$Description<-gsub("HALLMARK_","",total_table$Description)

GSEA_P1<-ggplot(total_table, aes(x=reorder(Description, Count), y=Count,fill=Type))+
  geom_bar(stat="identity",width=0.8,alpha=0.6)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  geom_text(aes(x=reorder(Description, Count),y=text_y,hjust =pos_y,label = Description),size=2.5)+
  theme_minimal()+
  ggtitle(paste("Hallmarks enrichment of PCs (PT vs. PN)",sep="-")) +
  theme()+
  xlab("Pathways")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+
  coord_flip()
print(GSEA_P1)

total_table$Count<-abs(total_table$Count)
write.table(total_table,"source_data_NC/Figure_4F.txt",sep = "\t",row.names = F,col.names = T)


##boxplot
combined1<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("PN","PT"))

combined$SubCluster<-factor(combined$SubCluster,levels=c(   "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                                                            "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,       "PVL_C1_mPC-ACTG2" ,
                                                               "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36"  ,"PVL_C4_imPC-MCAM", "PVL_C5-prolif" ))       

Idents(combined)<-combined$SubCluster
new.cluster.id<-c(   "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                     "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,    "mPC-ACTG2" ,
                       "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) 



names(new.cluster.id)<-levels(combined)
combined<-RenameIdents(combined,new.cluster.id)
combined$SubCluster1<-Idents(combined)


##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster1,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c(     "mPC-ACTG2" ,
                                                         "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) ),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(     "mPC-ACTG2" ,
                                                          "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) )

color_list<-c("#E64B35", "#4DBBD5")
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c("PN","PT"   )
sum_table$Type<-factor(sum_table$Type,levels=c("PN","PT"  ))
select_table<-sum_table
box_p<-ggplot(select_table, aes(Type, per, fill = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "t.test",label.y.npc=0.9)+
  theme_classic()+
  scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 9),
        strip.background = element_blank())+ 
  facet_wrap(~Cluster,ncol=3)

print(box_p)



write.table(select_table,"source_data_NC/Figure_4E.txt",sep = "\t",row.names = F,col.names = T)



##boxplot
combined1<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("PN","PT"))

combined$SubCluster<-factor(combined$SubCluster,levels=c(   "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                                                            "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,       "PVL_C1_mPC-ACTG2" ,
                                                            "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36"  ,"PVL_C4_imPC-MCAM", "PVL_C5-prolif" ))       

Idents(combined)<-combined$SubCluster
new.cluster.id<-c(   "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                     "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,    "mPC-ACTG2" ,
                     "mPC-MYH11" , "imPC-CD36"  ,"imPC-MCAM", "PC-prolif"   ) 



names(new.cluster.id)<-levels(combined)
combined<-RenameIdents(combined,new.cluster.id)
combined$SubCluster1<-Idents(combined)


##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster1,Type=combined$Type,Cancer.type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,Cancer.type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)
sum_table<-sum_table[sum_table$Cancer.type!="HNSCC",]
sum_table<-sum_table[which(sum_table$Cluster %in% c(  "imPC-MCAM") ),]

tab1<-dplyr::summarise(group_by(sum_table,Cancer.type),total_num=median(per))
tab1<-tab1[order(tab1$total_num,decreasing = T),]
sum_table$Cancer.type<-factor(sum_table$Cancer.type,levels=c(     tab1$Cancer.type ) )

color_list<-c("#4DBBD5","#E64B35")
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c("PN","PT"   )
sum_table$Type<-factor(sum_table$Type,levels=c("PN","PT"  ))
select_table<-sum_table
box_p11<-ggplot(select_table, aes(Type, per, color = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0.8),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,position=position_dodge(0.8))+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "t.test",label.y.npc=0.9)+
  theme_classic()+
  scale_color_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)',title = "MCAM+ imPC across cancer types") +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 9),
        strip.background = element_blank())+ 
  facet_wrap(~Cancer.type,ncol=6)

print(box_p11)


write.table(select_table,"source_data_NC/Figure_4H.txt",sep = "\t",row.names = F,col.names = T)




library(Seurat)
library(data.table)
library("sscVis")
library("sscClust")
library("scPip")
library(ggplot2)
library(ggridges)
library(dplyr)
library(cowplot)
library(PCAtools)
library(monocle)
library(ggpubr)
library(patchwork)
library(cowplot)


do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  #count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


library(plyr)

combined1<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("PN","PT"))

##combined$Type<-factor(combined$Type,levels=c("GN" ,  "GC", "GCPM", "PE" ,"PBMC"  ))
meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb
cellInfo.tb$Type<-paste(cellInfo.tb$Cancer.type,cellInfo.tb$Type,sep="_")

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Type,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 0
pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               legend_title ="O/e",
                               cluster_rows = T,
                               cluster_cols = T,
                               treeheight_row = 10,treeheight_col  = 10,
                               main="Tissue distribution",
                               angle_col = 45
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p2<-plot_grid(pheatmp_p2$gtable)

select_table<-data.frame(Cluster=gsub("PVL","PC",rownames(OR.dist.mtx)),round(OR.dist.mtx,2))
write.table(select_table,"source_data_NC/Figure_4G.txt",sep = "\t",row.names = F,col.names = T)




##survival

library(ggpubr)
library(magrittr)
library(ggsignif)
library(ggplot2)
library(ggsci)
library(RColorBrewer)




Gene_name="MCAM_imPC"
select_gene=Gene_name


RNA_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.htseq_fpkm-uq.tsv",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix)
gene_probe<-data.table::fread("~/data/TCGA/pancancer_expression/gencode.v22.annotation.gene.probeMap",header = T,stringsAsFactors = F)
phenotype_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.basic_phenotype.tsv",header = T,stringsAsFactors = F)
colnames(gene_probe)[1]<-"Ensembl_ID"
colnames(RNA_matrix)[1]<-"Ensembl_ID"
gene_probe<-gene_probe[,c(1,2)]
RNA_matrix<-merge(gene_probe,RNA_matrix,by="Ensembl_ID")
RNA_matrix<-data.frame(RNA_matrix)
RNA_matrix<-RNA_matrix[which(!(duplicated(RNA_matrix$gene))),]
rownames(RNA_matrix)<-RNA_matrix$gene

RNA_matrix<-RNA_matrix[,3:ncol(RNA_matrix)]
colnames(RNA_matrix)<-gsub("\\.","-",colnames(RNA_matrix))

Gene_name1="MCAM_imPC"
Gene_name2="SPP1_TAM"

mymatrix<-as.matrix(RNA_matrix)
Tip_Markers<-c("RGS5","MYL9","EGFL6","ACTA2","MCAM","PGF","ANGPT2","THY1","PDGFRB","COL4A1","ARHGDIB","NOTCH3")
Endo_Markers<-c("CD31","PECAM1")
mysymbol1<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)
mysymbol2<-data.frame(Gene_set="Endo",Gene_symbol=Endo_Markers)
mysymbol<-rbind(mysymbol1,mysymbol2)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)



type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs

library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=20)

##survplot
mt<-es.dif[1,]/es.dif[2,]
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=es.dif[1,],
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)


select_fpkm_matrix<-merge_matrix
colnames(select_fpkm_matrix)[1]<-"Tumor_ID"


library(limma)
select_fpkm_matrix$Tumor_ID<-gsub("\\.","-",select_fpkm_matrix$Tumor_ID)
Type_list<-strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,4]
select_fpkm_matrix<-select_fpkm_matrix[which(!(grepl("^1",Type_list))),]

select_fpkm_matrix$Tumor_ID<-paste(strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,1],strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,2],strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,3],sep="-")



##survival_analysis
surv_table<-read.table("~/data/TCGA/survival/TCGA_survival.txt",header = T,stringsAsFactors = F,sep="\t")
surv_table<-na.omit(surv_table)

library(survival)
library(survminer)
library(patchwork)
##Cut off by OS
total_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")

project_list<-c(  "TCGA-BLCA" ,
                  "TCGA-KIRP", "TCGA-LGG", "TCGA-MESO", 
                  "TCGA-STAD",
                  "TCGA-UVM" )


#project_list<-unique(sort(total_matrix$project_id))


sur_plot_list<-list()
risk_list<-list()
total_out<-data.frame()
for (i in 1:length(project_list)){
  Study_name<-project_list[[i]]
  merged_matrix<-total_matrix[total_matrix$project_id==Study_name,]
  res.cut <- surv_cutpoint(merged_matrix, #数据集
                           time = "OS.Time", #生存状态
                           event = "OS", #生存时间
                           variables = c("gene1") #需要计算的数据列名
  )
  merged_matrix$gene_level<-"Low"
  merged_matrix$gene_level[merged_matrix$gene1>res.cut$cutpoint$cutpoint]<-"High"
  merged_matrix$gene_level<-"Low"
  merged_matrix$gene_level[merged_matrix$gene1>quantile(merged_matrix$gene1,1/3)]<-"Mid"
  merged_matrix$gene_level[merged_matrix$gene1>quantile(merged_matrix$gene1,2/3)]<-"High"
  merged_matrix<-merged_matrix[which(merged_matrix$gene_level!="Mid"),]
  
  out_tab<-merged_matrix[,c("Tumor_ID","project_id","gene1","OS","OS.Time","gene_level")]
  total_out<-rbind(total_out,out_tab)
  
  surv_fit<-survfit(Surv(OS.Time , OS) ~ gene_level,data= merged_matrix)
  
  gg_surv1<-ggsurvplot(surv_fit,
                       conf.int = F,
                       #fun = "cumhaz",
                       linetype =  1, # Change line type by groups
                       size=0.5,
                       censor = F,
                       #surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_classic(),# Change ggplot2 theme
                       
                       palette = c("#E72C19", "#224FA2"),
                       title = paste(Study_name),
                       #font.family = "Arial",
                       #axis
                       xscale = "d_y",
                       pval = T,
                       surv.scale = "percent",
                       xlim = c(0, 2920), 
                       break.time.by=730.5,
                       xlab = "Time from diagnosis (years)",
                       #ylim = c(0, 0.05),
                       break.y.by=NULL,
                       ylab = "OS rate (%)",
                       #legend
                       legend = c(0.2,0.35),
                       legend.title = Gene_name,
                       legend.labs = c("High","Low"),
                       #risk table
                       risk.table = T,# Add risk table
                       font.risk.table=10,
                       fontsize=3,
                       
                       #字体
                       font.tickslab = c(10, "black"),
                       font.x = c(10, "black"),
                       font.y = c(10, "black"),
                       font.main = c(10, "black"),
                       font.legend = c(10, "black"),
  )
  gg_surv1$plot<-gg_surv1$plot+theme(plot.title = element_text(hjust = 0.5),
                                     legend.position = "none") 
  
  gg_surv1$table<-gg_surv1$table+
    theme(plot.title   = element_blank())+
    theme_classic2(
    base_size = 10)
  
  if(i>1&i!=7){
    gg_surv1$plot<-gg_surv1$plot+theme(axis.title.y   = element_blank()) 
     
  }
  if(i<7){
    gg_surv1$plot<-gg_surv1$plot+theme(axis.title.x   = element_blank()) 
    gg_surv1$table<-gg_surv1$table+theme(axis.title.x   = element_blank())
  }
  
  
  sur_plot_list[[length(sur_plot_list)+1]]<-gg_surv1$plot
  risk_list[[length(risk_list)+1]]<-gg_surv1$table
  
  print(gg_surv1)
}

total_out

write.table(total_out,"source_data_NC/Figure_4I.txt",sep = "\t",row.names = F,col.names = T)




##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/Angiogenesis_MCAM.csv",sep=",",header=T)

select_table$Type<-"Ctrl(1201Medium)"
select_table$Type[which(select_table$group=="pericyte_normaixa")]<-"PDGFRB+MCAM+ normoxia CM"
select_table$Type[which(select_table$group=="pericyte_hypoxia")]<-"PDGFRB+MCAM+ hypoxia CM"
select_table$Type[which(select_table$group=="pdgfrb+")]<-"PDGFRB+MCAM- normoxia CM"
select_table$Type[which(select_table$group=="Double_neg")]<-"PDGFRB-MCAM- normoxia CM"
select_table$Type[which(select_table$group=="HUVEC")]<-"HUVEC normoxia CM"

##Nodes
select_table$Type<-factor(select_table$Type,levels=c("Ctrl(1201Medium)","HUVEC normoxia CM",
                                                     "PDGFRB-MCAM- normoxia CM",
                                                     "PDGFRB+MCAM- normoxia CM",
                                                     "PDGFRB+MCAM+ normoxia CM",
                                                     "PDGFRB+MCAM+ hypoxia CM"
                                                     ))
select_table$mg.ml<-select_table$Nb.nodes

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.nodes<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM+ hypoxia CM"),
                                        c("Ctrl(1201Medium)","HUVEC normoxia CM"),
                                        c("PDGFRB+MCAM+ normoxia CM","PDGFRB-MCAM- normoxia CM")
                                        ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("HUVEC normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(
                       c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM- normoxia CM")
                     ),
                     size=3)+
  
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  #scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Number of nodes',title=paste ("Number of nodes")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,1500))

print(Nb.nodes)


write.table(select_table[,1:6],"source_data_NC/Figure_4K.txt",sep = "\t",row.names = F,col.names = T)




##Nb.Junctions
#select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
select_table$mg.ml<-select_table$Nb.Junctions

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.Junctions<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM+ hypoxia CM"),
                                        c("Ctrl(1201Medium)","HUVEC normoxia CM"),
                                        c("PDGFRB+MCAM+ normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("HUVEC normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(
                       c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM- normoxia CM")
                     ),
                     size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  #scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Number of junctions',title=paste ("Number of junctions")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,400))


print(Nb.Junctions)







##Tot.meshes.area
#select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot.meshes.area[which(select_table$Type=="Ctrl(1201Medium)")])
select_table$mg.ml<-select_table$Tot.meshes.area/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot.meshes.area<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM+ hypoxia CM"),
                                        c("Ctrl(1201Medium)","HUVEC normoxia CM"),
                                        c("PDGFRB+MCAM+ normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(
                                        c("HUVEC normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(
                       c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM- normoxia CM")
                     ),
                     size=3)+
  
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  #scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Area (fold change)',title=paste ("Total meshes area")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,10))


print(Tot.meshes.area)




##Tot.meshes.area
#select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot..lenght[which(select_table$Type=="Ctrl(1201Medium)")])
select_table$mg.ml<-select_table$Tot..lenght/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot..lenght<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM+ hypoxia CM"),
                                        c("Ctrl(1201Medium)","HUVEC normoxia CM"),
                                        c("PDGFRB+MCAM+ normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(c("HUVEC normoxia CM","PDGFRB-MCAM- normoxia CM")
                     ),
                     size=3)+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,
                     comparisons = list(
                       c("PDGFRB+MCAM+ normoxia CM","PDGFRB+MCAM- normoxia CM")
                     ),
                     size=3)+
  geom_jitter(width = 0.2 )+
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  labs(x = 'Type', y = 'Length (fold change)',title=paste ("Total length")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5))


print(Tot..lenght)



mid_p<-Nb.nodes+Nb.Junctions+Tot..lenght+Tot.meshes.area+plot_layout(ncol=2)
mid_p





total_p4<-ggarrange(plotlist =sur_plot_list,ncol = 7,nrow = 5,
                    common.legend = T,legend = "bottom")
total_p41<-ggarrange(plotlist =c(sur_plot_list[1:6],risk_list[1:6]),ncol = 6,nrow = 2,heights = c(2,1),
                    common.legend = T,legend = "bottom")
total_p4<-sur_plot_list[[1]]+sur_plot_list[[2]]+sur_plot_list[[3]]+sur_plot_list[[4]]+sur_plot_list[[5]]+sur_plot_list[[6]]+
  risk_list[[1]]+risk_list[[2]]+risk_list[[3]]+risk_list[[4]]+risk_list[[5]]+risk_list[[6]]+
  plot_layout(ncol = 6,nrow =2,heights = c(2,1))


library(patchwork)

total_pn<-umap_p1+nh_graph_endo+beeswarm_endo+plot_layout(ncol=3,widths = c(2,2,1))
total_p1<-ggarrange(dot_p,box_p,GSEA_P1,ncol=3,nrow = 1,widths =c(2,1.2,1.3) )
total_p3<-ggarrange(pheatmp_p2,box_p11,widths = c(1.5,1))
total_p5<-ggarrange(mid_p,ncol = 2,widths = c(2,1))
total_p<-ggarrange(total_pn,total_p1,total_p3,total_p4,total_p5,nrow = 5,ncol = 1,heights = c(3.5,2.8,3.5,3 ,3.5 ))
pdf("14.Figure/4.Figure_4.pdf",width = 12,height = 16.2)
print(total_p)
dev.off()

