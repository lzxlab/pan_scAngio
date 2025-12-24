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
                                                    "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PC",     "Endo" ,   "Epi"))

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

cor_data<-cor_p1$data
write.table(cor_data,"source_data_NC/Figure_3A.txt",sep = "\t",row.names = F,col.names = T)



cell_pair_list<-list(c("PC","Endo"))

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
cor_p2<-p1

colnames(merged_table)<-c("Sample_ID","PC_percentage","Endo_percentage")
write.table(merged_table,"source_data_NC/Figure_3B.txt",sep = "\t",row.names = F,col.names = T)




combined1<-combined
combined1$Cluster<-as.character(combined1$Cluster)
combined1$Cluster[which(grepl("PC",combined1$Cluster))]<-combined1$SubCluster[which(grepl("PC",combined1$Cluster))]

tmp_table<-data.frame(Cluster=combined1$Cluster,Type=combined1$Type,sample_ID=combined1$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)



cell_pair_list<-list(c("PC_C1_mPC-ACTG2","Endo"),
                     c("PC_C2_mPC-MYH11","Endo"),
                     c("PC_C3_imPC-CD36","Endo"),
                     c("PC_C4_imPC-MCAM","Endo"))


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


cell_pair_list<-list(c("PC","Tip_ECs"),
                     c("PC","Immature_ECs"),
                     c("PC","Artery_ECs"),
                     c("PC","Venous_ECs"),
                     c("PC","Capillary_ECs"),
                     c("PC","LECs"))
total_outputdata<-data.frame()
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
  total_outputdata<-rbind(total_outputdata,
                          data.frame(Cluster1=Cluster1,Cluster2=Cluster2,
                                     merged_table)
  )
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

colnames(total_outputdata)<-c("Cluster1","Cluster2","sample_ID",
                              "percentage_Cluster1", "percentage_Cluster2")
total_outputdata<-total_outputdata[,c("sample_ID","Cluster1","percentage_Cluster1",
                                      "Cluster2",
                                      "percentage_Cluster2")]

write.table(total_outputdata,"source_data_NC/Figure_3H.txt",sep = "\t",row.names = F,col.names = T)




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

netVisual_b1_data<-netVisual_b1$data
netVisual_b1_data$source<-as.character(netVisual_b1_data$source)
netVisual_b1_data$source.target<-as.character(netVisual_b1_data$source.target)
netVisual_b1_data$source[which(netVisual_b1_data$source=="PVL")]<-"PC"
netVisual_b1_data$source.target[which(netVisual_b1_data$source.target=="PVL")]<-"PC"
netVisual_b1_data$source.target<-gsub("PVL","PC",netVisual_b1_data$source.target)

write.table(netVisual_b1_data,"source_data_NC/Figure_3D_right.txt",sep = "\t",row.names = F,col.names = T)


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


netVisual_b11_data<-netVisual_b11$data
netVisual_b11_data$source<-as.character(netVisual_b11_data$source)
netVisual_b11_data$target<-as.character(netVisual_b11_data$target)
netVisual_b11_data$source.target<-as.character(netVisual_b11_data$source.target)
netVisual_b11_data$source[which(netVisual_b11_data$source=="PVL")]<-"PC"
netVisual_b11_data$target[which(netVisual_b11_data$target=="PVL")]<-"PC"
netVisual_b11_data$source.target[which(netVisual_b11_data$source.target=="PVL")]<-"PC"
netVisual_b11_data$source.target<-gsub("PVL","PC",netVisual_b11_data$source.target)

write.table(netVisual_b11_data,"source_data_NC/Figure_3D_left.txt",sep = "\t",row.names = F,col.names = T)




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


cell_chat_p2_data<-cell_chat_p2$data
cell_chat_p2_data$source<-as.character(cell_chat_p2_data$source)
cell_chat_p2_data$target<-as.character(cell_chat_p2_data$target)
cell_chat_p2_data$source.target<-as.character(cell_chat_p2_data$source.target)
cell_chat_p2_data$source[which(cell_chat_p2_data$source=="PVL")]<-"PC"
cell_chat_p2_data$target[which(cell_chat_p2_data$target=="PVL")]<-"PC"
cell_chat_p2_data$source.target<-gsub("PVL","PC",cell_chat_p2_data$source.target)
cell_chat_p2_data$group.names<-gsub("PVL","PC",cell_chat_p2_data$group.names)

write.table(cell_chat_p2_data,"source_data_NC/Figure_3E.txt",sep = "\t",row.names = F,col.names = T)




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
                                                        "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PC",     "Endo" ,   "Epi")
)
Idents(combined)<-combined$Cluster
Vln_exp<-VlnPlot(combined,features = c("ANGPT2","TEK","PGF","FLT1"),pt.size=0,group.by="Cluster",ncol=2) &
  theme(axis.text = element_text(size = 10), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        #axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))

data1<-data.frame(Cluster=Vln_exp[[1]]$data$ident,
                  ANGPT2=Vln_exp[[1]]$data$ANGPT2,
                  TEK=Vln_exp[[2]]$data$TEK,
                  PGF=Vln_exp[[3]]$data$PGF,
                  FLT1=Vln_exp[[4]]$data$FLT1
                  )
write.table(data1,"source_data_NC/Figure_3G.txt",sep = "\t",row.names = F,col.names = T)







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

write.table(cor_dot_plot$data,"source_data_NC/Figure_3F.txt",sep = "\t",row.names = F,col.names = T)





##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/Angiogenesis_plot.csv",sep=",",header=T)
select_table$Type<-as.character(select_table$group)


##Nodes
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
select_table$mg.ml<-select_table$Nb.nodes

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.nodes<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
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


write.table(select_table[,1:6],"source_data_NC/Figure_3J.txt",sep = "\t",row.names = F,col.names = T)




##Nb.Junctions
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
select_table$mg.ml<-select_table$Nb.Junctions

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.Junctions<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
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
  scale_y_continuous(expand = c(0,0),limits = c(0,500))


print(Nb.Junctions)







##Tot.meshes.area
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot.meshes.area[which(select_table$Type=="normoxia")])
select_table$mg.ml<-select_table$Tot.meshes.area/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot.meshes.area<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
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
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))


print(Tot.meshes.area)




##Tot.meshes.area
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot..lenght[which(select_table$Type=="normoxia")])
select_table$mg.ml<-select_table$Tot..lenght/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot..lenght<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
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
  scale_y_continuous(expand = c(0,0),limits = c(0,1.4))


print(Tot..lenght)



mid_p<-Nb.nodes+Nb.Junctions+Tot..lenght+Tot.meshes.area+plot_layout(ncol=2)
mid_p


library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/pericyte_tumor_curve.csv",sep=",",header = T)

g1<-tab$X19[which(tab$group=="HCT116+Pericyte")]
g2<-tab$X19[which(tab$group=="HCT116")]
color_list<-c("#E41A1C","#377EB8","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))
tab1$Group<-"HCT116"
tab1$Group[tab1$group=="HCT116+Pericyte"]<-"HCT116+PCs"
tab1$Group<-factor(tab1$Group,levels=c("HCT116+PCs" , "HCT116"))
library(dplyr)
result_table<-dplyr::summarize(group_by(tab1,Group,Time),
                               mean_var=mean(Value),
                               SD=sd(Value),
                               N_N=n(),  
                               se=SD/sqrt(N_N),  
                               upper_limit=mean_var+SD,  
                               lower_limit=mean_var-SD)



mice_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
  geom_point(data = result_table,aes(x=Time, y=mean_var),pch=15,size=1.5)+
  geom_line(size=0.5)+
  geom_errorbar(data = result_table,aes(ymin = lower_limit, ymax=upper_limit), #误差条表示均值±标准差
                width=0.2, #误差条末端短横线的宽度
                alpha = 0.7,
                size=0.5) +
  scale_color_manual(values = color_list)+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve\n(HCT116 with or without PCs)")+
  theme_classic2()+
  theme(legend.position=c(0.3,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8),breaks =seq(0,0.8,0.2))
print (mice_line_p1)


write.table(result_table,"source_data_NC/Figure_3K_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(1,3,4)],"source_data_NC/Figure_3K_data.txt",sep = "\t",row.names = F,col.names = T)






##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/pericyte_tumoe_weight.csv",sep=",",header=T)
select_table$Type<-"HCT116"
select_table$Type[select_table$group=="HCT116+Pericyte"]<-"HCT116+PCs"
select_table$Type<-factor(select_table$Type,levels=c("HCT116+PCs" , "HCT116"))


##Nodes
select_table$Type<-factor(select_table$Type,levels=c( "HCT116","HCT116+PCs"))
select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
PC_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c( "#377EB8", "#E41A1C"))+
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8))

print(PC_weight)

write.table(select_table[1:2],"source_data_NC/Figure_3L_top.txt",sep = "\t",row.names = F,col.names = T)



select_table<-read.table("Experiments/pericyte_CD31.csv",sep=",",header=T)
select_table$Type<-"HCT116"
select_table$Type[select_table$group=="HCT116+Pericyte"]<-"HCT116+PCs"
select_table$Type<-factor(select_table$Type,levels=c("HCT116+PCs" , "HCT116"))


##Nodes
select_table$Type<-factor(select_table$Type,levels=c( "HCT116","HCT116+PCs"))
select_table$mg.ml<-select_table$X.area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
PC_CD31<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c( "#377EB8", "#E41A1C"))+
  labs(x = 'Type', y = 'CD31 density',title=paste ("CD31 density")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,2))

print(PC_CD31)

write.table(select_table[1:2],"source_data_NC/Figure_3L_bottom.txt",sep = "\t",row.names = F,col.names = T)



mid_p<-Nb.nodes+Nb.Junctions+Tot..lenght+Tot.meshes.area+plot_layout(ncol=2)
right_p<-PC_weight+PC_CD31+plot_layout(ncol=1)
mid_p


total_p6<-ggarrange(mid_p,mice_line_p1,right_p,ncol=4,widths = c(2,1.5,1.4,0.7))




total_p1<-ggarrange(cor_p1,ggarrange(netVisual_b11+netVisual_b1+plot_layout(ncol=2,widths = c(2,1)),
                                     cell_chat_p2,ncol = 1,nrow = 2,heights = c(1,1.5)),ncol=2)
total_p2<-ggarrange(netVisual_b1,cell_chat_p2,ncol=2,widths = c(1,1.2))
total_p3<-ggarrange(Vln_exp,cor_p4,ncol=2,widths = c(1,1))

total_p5<-plot_spacer()+cor_dot_plot+plot_layout(ncol = 2,widths =c(3,1) )
total_p<-ggarrange(total_p1,total_p5,total_p3,total_p6,ncol=1,nrow=4,heights = c(5,3.5,4,3))
pdf("14.Figure/3.Figure_3.pdf",width = 12,height = 15.5)
print(total_p)
dev.off()


pdf("14.Figure/3.Figure_3B.pdf",width = 2,height = 2)
print(cor_p2)
dev.off()

