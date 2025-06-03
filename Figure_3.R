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
library(patchwork)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")



##correlation
##Subcluster
sub_sce1<-readRDS("3.Cluster/13.Annotation/non_immune.rds")
sub_sce1$Cluster<-sub_sce1$TopCluster

DefaultAssay(sub_sce1)<-"RNA"
markers1<-c("VEGFA","VEGFB","VEGFC","PGF","ANGPT2","HGF","FGF1" , "FGF2" ,"FGF3" , "FGF4",  "FGF5",  "FGF6" , "FGF7" , "FGF8" , "FGF9" , "FGF11", "FGF12" ,"FGF13" ,"FGF14" ,"FGF16", "FGF17", "FGF18" ,"FGF19"  , "FGF20", "FGF21",
            "FGF22", "FGF23" ,"MET")
markers2<-c("FLT1","KDR","FLT4","FGFR1","FGFR2","FGFR3","FGFR4")

marker_list<-c(markers1,markers2)
sub_sce1<-subset(sub_sce1,features=marker_list)


sub_sce2<-readRDS("3.Cluster/13.Annotation/immune.rds")
sub_sce2$Cluster<-sub_sce2$MidCluster

DefaultAssay(sub_sce2)<-"RNA"
markers1<-c("VEGFA","VEGFB","VEGFC","PGF","ANGPT2","HGF","FGF1" , "FGF2" ,"FGF3" , "FGF4",  "FGF5",  "FGF6" , "FGF7" , "FGF8" , "FGF9" , "FGF11", "FGF12" ,"FGF13" ,"FGF14" ,"FGF16", "FGF17", "FGF18" ,"FGF19"  , "FGF20", "FGF21",
            "FGF22", "FGF23" ,"MET")
markers2<-c("FLT1","KDR","FLT4","FGFR1","FGFR2","FGFR3","FGFR4")

marker_list<-c(markers1,markers2)
sub_sce2<-subset(sub_sce2,features=marker_list)



combined<-merge(sub_sce1,sub_sce2)

combined$Cluster<-factor(combined$Cluster,levels=c( "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                    "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"))

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p


library(ggcorrplot2)
cor_p1<-ggcorrplot(cor_value, method = "ellipse",type = "lower",p.mat = cor_test_mat, 
                   insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001),col = rev(c("#CA0020" ,"#F4A582", "#F7F7F7", "#92C5DE" ,"#0571B0")))



cell_pair_list<-list(c("PVL","Endo"),
                     c("Endo","CD4_conv"),
                     c("Endo","CD4_Treg"),
                     c("Endo","CD8"))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=1,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste("Correlation of",Cluster1,"and",Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(plot.title = element_text(size = 10,hjust = 0.5),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.line = element_blank(), 
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_blank(), 
          legend.title = element_text(size = 10),
          legend.position = "none",
          legend.text = element_text(size = 9))+ 
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p2<-p1+p2+p3+p4+plot_layout(ncol=2)
cor_p2<-p1





combined1<-combined
combined1$Cluster<-as.character(combined1$Cluster)
combined1$Cluster[which(grepl("PVL",combined1$Cluster))]<-combined1$SubCluster[which(grepl("PVL",combined1$Cluster))]

tmp_table<-data.frame(Cluster=combined1$Cluster,Type=combined1$Type,sample_ID=combined1$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)



cell_pair_list<-list(c("PVL_C1_mPC-ACTG2","Endo"),
                     c("PVL_C2_mPC-MYH11","Endo"),
                     c("PVL_C3_imPC-CD36","Endo"),
                     c("PVL_C4_imPC-MCAM","Endo"))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=1,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste("Correlation of",Cluster1,"and",Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p3<-p1+p2+p3+p4+plot_layout(ncol=2)






combined1<-combined
combined1$Cluster<-as.character(combined1$Cluster)
combined1$Cluster[which(grepl("Endo",combined1$Cluster))]<-combined1$MidCluster[which(grepl("Endo",combined1$Cluster))]

tmp_table<-data.frame(Cluster=combined1$Cluster,Type=combined1$Type,sample_ID=combined1$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p


library(ggcorrplot2)
cor_p4<-ggcorrplot(cor_value, method = "ellipse",type = "lower",p.mat = cor_test_mat, 
                   insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001))


cell_pair_list<-list(c("PVL","Tip_ECs"),
                     c("PVL","Immature_ECs"),
                     c("PVL","Artery_ECs"),
                     c("PVL","Venous_ECs"),
                     c("PVL","Capillary_ECs"),
                     c("PVL","LECs"))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=1,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste(Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p4<-p1+p2+p3+p4+p5+p6+plot_layout(ncol=3)





library(CellChat)
##Cell_chat
cellchat<-readRDS("6.cell_chat/2.All/cellchat.RDS")
#cellchat@idents<-factor(cellchat@idents,levels=c( "B" ,"Plasma"   ,    "CD4_conv"   ,   "CD4_Treg",  "CD8", "MAIT"  ,   "NK"   ,    "DC"  , "Macro"  ,    "Mono"  ,  
#                                                  "Mast" ,    "Endo", "APC","pAD"  ,  "SMC"            )) 

groupSize <- as.numeric(table(cellchat@idents))
cellchat@idents<-factor(cellchat@idents,levels=c(     "DC" , "Mono"   , "Macro" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi")
)
levels(cellchat@idents)
vertex.receiver = seq(1,4)


pathways.show <- "VEGF"
select_list<-"VEGF|CSF|ANGPT2|ANGPT1"
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])
#pairLR.use$interaction_name<-unique(sort(pairLR.use$interaction_name))
pairLR.use<-pairLR.use[!(grepl("CD99|THBS1|HLA",pairLR.use$interaction_name)),]

pairLR.use<-data.frame(interaction_name=pairLR.use)
netVisual_b1<-netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = "PVL",remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "none",
        legend.text = element_text(size = 8))

netVisual_b11<-netVisual_bubble(cellchat, pairLR.use = pairLR.use, 
                                sources.use = c(     "DC" , "Mono"   , "Macro" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"),
                                targets.use  = "Endo",remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "none",
        legend.text = element_text(size = 8))





library(CellChat)
library(patchwork)

cellchat1<-readRDS("6.cell_chat/2.All/PN.RDS")
cellchat2<-readRDS("6.cell_chat/2.All/PT.RDS")



cellchat <- mergeCellChat(list(cellchat1, cellchat2), add.names = c("PN", "PT"))

cellchat@idents<-factor(cellchat@idents,levels=c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,  
                                                        "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi")
)

groupSize <- as.numeric(table(cellchat1@idents))
levels(cellchat1@idents) 
vertex.receiver = seq(1,4) 




pathways.show <- "VEGF"
select_list<-"VEGF|ANGPT|DLL"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])
pairLR.use<-pairLR.use[!(grepl("CD99|THBS1|HLA",pairLR.use$interaction_name)),]

pairLR.use<-data.frame(interaction_name=pairLR.use)

cell_chat_p2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,  comparison = c(1, 2) , sources.use = "PVL", targets.use = c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,
                                                                                                                                   "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"), remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8))




##correlation
##Subcluster
sub_sce1<-readRDS("3.Cluster/13.Annotation/non_immune.rds")
sub_sce1$Cluster<-sub_sce1$TopCluster

combined<-sub_sce1
if (T){
  sample_ID<-runif(20000,1,ncol(combined))
  combined$merge_var<-"Drop"
  combined$merge_var[sample_ID]<-"Keep"
  Idents(combined)<-combined$merge_var
  combined<-subset(combined,idents="Keep")
}

sub_sce1<-combined


sub_sce2<-readRDS("3.Cluster/13.Annotation/immune.rds")
sub_sce2$Cluster<-sub_sce2$MidCluster

combined<-sub_sce2
if (T){
  sample_ID<-runif(20000,1,ncol(combined))
  combined$merge_var<-"Drop"
  combined$merge_var[sample_ID]<-"Keep"
  Idents(combined)<-combined$merge_var
  combined<-subset(combined,idents="Keep")
}

sub_sce2<-combined


combined<-merge(sub_sce1,sub_sce2)
combined$Cluster<-factor(combined$Cluster,levels=c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                        "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi")
)
Idents(combined)<-combined$Cluster
Vln_exp<-VlnPlot(combined,features = c("ANGPT2","TEK","PGF","FLT1"),pt.size=0,group.by="Cluster",ncol=2) &
  theme(axis.text = element_text(size = 10), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        #axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))




combined<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")


Idents(combined)<-combined$TopCluster
combined<-subset(combined,idents="PVL")

data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("PGF")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Type),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]))


box_p2<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of PGF")+
  xlab("Expression level")+
  scale_fill_manual(values = viridis::viridis(5,option = "D"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlim(c(-2 ,6))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p2



gene_name<- c("ANGPT2")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Type),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]))


box_p3<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of ANGPT2")+
  xlab("Expression level")+
  scale_fill_manual(values = viridis::viridis(5,option = "D"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlim(c(-2 ,6))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p3








##Bevtreatment
data<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/GSE140082_geo.normdata.csv",sep=",",header = T)

colnames(data)[1]<-"ID"

ID_map<-data.table::fread("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/GPL14951-11332.txt")
ID_map<-ID_map[,c("ID","ILMN_Gene")]
data<-merge(ID_map,data,by="ID")

data<-data[which(!(duplicated(data$ILMN_Gene))),]
data<-data.frame(data)
rownames(data)<-data$ILMN_Gene

data<-data[,which(!(grepl("pval",colnames(data))))]
data<-data[,c(-1,-2)]



Gene_name="EGFL6+ imPC"
select_gene=Gene_name


##survival_analysis
surv_table<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/Survival.txt",header = T,stringsAsFactors = F,sep="\t")
surv_table<-surv_table[which(surv_table$Treatment=="bevacizumab"),]
surv_table<-na.omit(surv_table)


mymatrix<-as.matrix(data)

Tip_Markers<-c("RGS5","MYL9","EGFL6","ACTA2","MCAM","PGF","ANGPT2","THY1")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=10)


mt<-es.dif[1,]

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))


##Cut off by OS
merged_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")

library(survminer)
library(survival)
res.cut <- surv_cutpoint(merged_matrix, #数据集
                         time = "OS_Time", #生存状态
                         event = "OS", #生存时间
                         variables = c("gene") #需要计算的数据列名
)
merged_matrix$gene_level<-"Low"
merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"

surv_fit<-survfit(Surv(OS_Time , OS) ~ gene_level,data= merged_matrix)

Bev_gg_surv1<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#E72C19", "#224FA2"),
                         title = paste("OS by ",Gene_name,",\nBev treatments cohort",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "d_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 1470), 
                         break.time.by=365.25,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "OS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("High","Low"),
                         #risk table
                         risk.table = F,# Add risk table
                         
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv1$plot<-Bev_gg_surv1$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_blank(),
                                                           legend.position = "none")
print (Bev_gg_surv1)





surv_fit<-survfit(Surv(PFS_Time , PFS) ~ gene_level,data= merged_matrix)

Bev_gg_surv2<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#E72C19", "#224FA2"),
                         title = paste("PFS by ",Gene_name,",\nBev treatments cohort",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "d_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 1470), 
                         break.time.by=365.25,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "PFS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("High","Low"),
                         #risk table
                         risk.table = F,# Add risk table
                         
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv2$plot<-Bev_gg_surv2$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_blank(),
                                                           legend.position = "right")
print (Bev_gg_surv2)



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
Tip_Markers<-c("MYL9","GUCY1B3","ACTA2","SOD3","CALD1","TPM2","TAGLN","EGFL6","RGS5","PPP1R14A","SEPT4","HIGD1B","MFGE8","GJA4","FOXS1","C1QTNF1","HEYL","GUCY1A3","COX4I2","KCNMB1","ADAP2","RCSD1","BGN","FSCN1","SPARCL1","CRISPLD2","S100A10","SLC38A11","HEY2","NOTCH3","WFDC1","PTP4A3","TPM1","PERP","GEM","COL4A2","GPR4","CSRP2","KCNA5","C11orf96","MYLK","HSPB6","EFHD1","DKK3","MCAM","S100A16","FRZB","ESAM","RASL12","PLXDC1","EDNRA","FOXF1","A2M")
Endo_Markers<-c("PLVAP","AQP1","PECAM1","RNASE1","CD93","EMCN","SRGN","ADGRL4","FBLN1","GIMAP4","PTPRB","IL3RA","ARHGAP29","PITX2","VWF","ESAM","CDH5","S1PR1","CLEC14A","RAMP3","TIE1","ICAM2","FLT1","CXorf36","CD34","CYYR1","PODXL","FLI1","CCL14","CLDN5","RAMP2","ENG","LDB2","PCDH17","HYAL2","NPDC1","ADCY4","SLCO2A1","TMEM255B","MCTP1","RND1","PLTP","GNG11","CD200","ERG","ACKR1","TGFBR2","COL4A1","COL15A1","EGFL7")
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


project_list<-c(  "TCGA-BLCA" ,"TCGA-BRCA", "TCGA-CESC", "TCGA-COAD" , "TCGA-GBM" ,
                  "TCGA-KIRP", "TCGA-LIHC",  "TCGA-MESO",
                  "TCGA-PAAD","TCGA-STAD", 
                  "TCGA-UCEC","TCGA-UVM" )


select_matrix<-merge_matrix
project_id_list<-unique(select_matrix$project_id)
project_id_list<-project_id_list[which((project_id_list %in% c(project_list)))]
select_table<-select_matrix[select_matrix$project_id %in% project_id_list, ]


#绘图
TCGA_cor_p1 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("ssGSEA score (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of PCs and ECs, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=6)
print(TCGA_cor_p1)







library(ggpubr)
library(magrittr)
library(ggsignif)
library(ggplot2)
library(ggsci)
library(RColorBrewer)




Gene_name="VEGFA"
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

Gene_name1="PGF"
Gene_name2="ECs"

mymatrix<-as.matrix(RNA_matrix)
Tip_Markers<-c("MYL9","GUCY1B3","ACTA2","SOD3","CALD1","TPM2","TAGLN","EGFL6","RGS5","PPP1R14A","SEPT4","HIGD1B","MFGE8","GJA4","FOXS1","C1QTNF1","HEYL","GUCY1A3","COX4I2","KCNMB1","ADAP2","RCSD1","BGN","FSCN1","SPARCL1","CRISPLD2","S100A10","SLC38A11","HEY2","NOTCH3","WFDC1","PTP4A3","TPM1","PERP","GEM","COL4A2","GPR4","CSRP2","KCNA5","C11orf96","MYLK","HSPB6","EFHD1","DKK3","MCAM","S100A16","FRZB","ESAM","RASL12","PLXDC1","EDNRA","FOXF1","A2M")
Endo_Markers<-c("PLVAP","AQP1","PECAM1","RNASE1","CD93","EMCN","SRGN","ADGRL4","FBLN1","GIMAP4","PTPRB","IL3RA","ARHGAP29","PITX2","VWF","ESAM","CDH5","S1PR1","CLEC14A","RAMP3","TIE1","ICAM2","FLT1","CXorf36","CD34","CYYR1","PODXL","FLI1","CCL14","CLDN5","RAMP2","ENG","LDB2","PCDH17","HYAL2","NPDC1","ADCY4","SLCO2A1","TMEM255B","MCTP1","RND1","PLTP","GNG11","CD200","ERG","ACKR1","TGFBR2","COL4A1","COL15A1","EGFL7")
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

Gene_name1="PGF"
Gene_name2="PCs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix1<-merge_matrix

Gene_name1="ANGPT2"
Gene_name2="PCs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix2<-merge_matrix


project_id_list<-unique(merge_matrix2$project_id)
total_tab<-data.frame()
for (i in 1:length(project_id_list)){
  project_id<-project_id_list[[i]]
  
  
  merged_table<-merge_matrix1[which(merge_matrix1$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  cor.value1<-cor.value
  p.value1<-p.value
  
  merged_table<-merge_matrix2[which(merge_matrix2$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  
  cor.value2<-cor.value
  p.value2<-p.value
  
  tmp_table<-data.frame(Cancer.type=project_id,
                        PGF.cor.value=cor.value1,
                        PGF.p.value=p.value1,
                        ANGPT2.cor.value=cor.value2,
                        ANGPT2.p.value=p.value2
  )
  total_tab<-rbind(total_tab,tmp_table)
  
}

total_tab<-total_tab[order(total_tab$PGF.cor.value,decreasing = T),]
total_tab$lable<-NA
total_tab$lable[1:5]<-gsub("TCGA-","",total_tab$Cancer.type)[1:3]
total_tab<-total_tab[order(total_tab$ANGPT2.cor.value,decreasing = T),]
total_tab$lable[1:5]<-gsub("TCGA-","",total_tab$Cancer.type)[1:3]

cor_dot_plot<-ggplot(total_tab,aes(x=PGF.cor.value,y=ANGPT2.cor.value,color=Cancer.type))+
  geom_point()+
  geom_text_repel(aes(label = lable))+
  theme_classic()+
  ylim(c(-1,1))+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0,linetype="dashed")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ggtitle("Correlation of PCs and PGF/ANGPT2\nTCGA pan-caner dataset (RNA-seq)")+
  xlab(paste0("Cor of PCs and PGF"))+
  ylab(paste0("Cor of PCs and ANGPT2"))+
  theme(legend.position = 'none',
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
        axis.text.x = element_text(size = 10,vjust =0.5,hjust = 0.5),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.text  = element_text(size = 10),
        axis.line = element_blank()
  )

cor_dot_plot

total_p1<-ggarrange(cor_p1,ggarrange(netVisual_b11+netVisual_b1+plot_layout(ncol=2,widths = c(2,1)),
                                     cell_chat_p2,ncol = 1,nrow = 2,heights = c(1,1.5)),ncol=2)
total_p2<-ggarrange(netVisual_b1,cell_chat_p2,ncol=2,widths = c(1,1.2))
total_p3<-ggarrange(Vln_exp,cor_p4,ncol=2,widths = c(1,1))
total_p4<-box_p2+box_p3+Bev_gg_surv1$plot+Bev_gg_surv2$plot+plot_layout(ncol=4,widths = c(1,1,2,2))

total_p5<-plot_spacer()+cor_dot_plot+plot_layout(ncol = 2,widths =c(3,1) )
total_p<-ggarrange(total_p1,total_p5,total_p3,ncol=1,nrow=4,heights = c(5,3.5,4,3))
pdf("14.Figure/6.Figure_6.pdf",width = 12,height = 15.5)
print(total_p)
dev.off()


pdf("14.Figure/Figure_3.pdf",width = 2,height = 2)
print(cor_p2)
dev.off()
