##BOXplot
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
library(RColorBrewer)
library(ggpubr)
library(easyGgplot2)
library(ggpubr)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")
combined1<-readRDS("3.Cluster/13.Annotation/non_immune.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("PT"))

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c(  "Epi_C5-NDRG1" ,"Epi_C10-SST"   ) ),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(  "Epi_C5-NDRG1" ,"Epi_C10-SST"    ) )

color_list<-c(RColorBrewer::brewer.pal(9,"PuBuGn")[3:9],RColorBrewer::brewer.pal(9,"YlOrRd")[3:9])

cluster_list<-c(  "Epi_C5-NDRG1" ,"Epi_C10-SST"  )

ps_list<-list()
for (cluster in cluster_list){
  #cluster="Epi_C5-NDRG1"
  select_table<-sum_table[sum_table$Cluster==cluster,]
  order_tab<-dplyr::summarize(group_by(select_table,Type),per=median(per))
  order_tab<-order_tab[order(order_tab$per,decreasing = T),]
  select_table$Type<-factor(select_table$Type,levels=order_tab$Type)
  
  box_p<-ggplot(select_table, aes(Type, per, color = Type))+
    geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
    geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
    #ggrepel::geom_label_repel(aes(label=sample_ID))+
    stat_compare_means(label.y.npc = 0.9,label.x.npc = 0.1)+
    theme_classic()+
    scale_color_manual(values = color_list)+
    labs(x = 'Type of tissue', y = 'Relative Abundance(%)',title = cluster) +
    #theme(panel.grid =a element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust = 1),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          strip.background = element_blank())
  
  ps_list[[length(ps_list)+1]]<-box_p
  
}





##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c("iCAF_C2-MMP3"  ,"PVL_C4_imPC-MCAM" ,"PVL_C5-prolif") ),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c( "iCAF_C2-MMP3"  , "PVL_C4_imPC-MCAM" ,"PVL_C5-prolif"  ) )

color_list<-c(RColorBrewer::brewer.pal(9,"PuBuGn")[3:9],RColorBrewer::brewer.pal(9,"YlOrRd")[3:9])

cluster_list<-c("iCAF_C2-MMP3"  ,"PVL_C4_imPC-MCAM" ,"PVL_C5-prolif" )

#ps_list<-list()
for (cluster in cluster_list){
  #cluster="Tip_ECs"
  select_table<-sum_table[sum_table$Cluster==cluster,]
  order_tab<-dplyr::summarize(group_by(select_table,Type),per=median(per))
  order_tab<-order_tab[order(order_tab$per,decreasing = T),]
  select_table$Type<-factor(select_table$Type,levels=order_tab$Type)
  
  box_p<-ggplot(select_table, aes(Type, per, color = Type))+
    geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
    geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
    #ggrepel::geom_label_repel(aes(label=sample_ID))+
    stat_compare_means(label.y.npc = 0.9,label.x.npc = 0.1)+
    theme_classic()+
    scale_color_manual(values = color_list)+
    labs(x = 'Type of tissue', y = 'Relative Abundance(%)',title = cluster) +
    #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust = 1),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          strip.background = element_blank())
  
  ps_list[[length(ps_list)+1]]<-box_p
  
}

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)
out_p1<-sum_table[which(sum_table$Cluster %in% c("iCAF_C2-MMP3"  ,"PC_C4_imPC-MCAM" ,
                                                 "PC_C5-prolif","Epi_C5-NDRG1" ,"Epi_C10-SST" ) ),]



combined1<-readRDS("3.Cluster/13.Annotation/immune.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("PT"))

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c(  "Mono_C3-IL1B" ,"Mono_C2-RTEN" ,"TAM_C5-SPP1") ),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(  "Mono_C3-IL1B" ,"Mono_C2-RTEN" ,"TAM_C5-SPP1"  ) )

color_list<-c(RColorBrewer::brewer.pal(9,"PuBuGn")[3:9],RColorBrewer::brewer.pal(9,"YlOrRd")[3:9])

cluster_list<-c( "Mono_C3-IL1B" ,"Mono_C2-RTEN" ,"TAM_C5-SPP1" )

#ps_list<-list()
for (cluster in cluster_list){
  #cluster="Tip_ECs"
  select_table<-sum_table[sum_table$Cluster==cluster,]
  order_tab<-dplyr::summarize(group_by(select_table,Type),per=median(per))
  order_tab<-order_tab[order(order_tab$per,decreasing = T),]
  select_table$Type<-factor(select_table$Type,levels=order_tab$Type)
  
  box_p<-ggplot(select_table, aes(Type, per, color = Type))+
    geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
    geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
    #ggrepel::geom_label_repel(aes(label=sample_ID))+
    stat_compare_means(label.y.npc = 0.9,label.x.npc = 0.1)+
    theme_classic()+
    scale_color_manual(values = color_list)+
    labs(x = 'Type of tissue', y = 'Relative Abundance(%)',title = cluster) +
    #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust = 1),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          strip.background = element_blank())
  
  ps_list[[length(ps_list)+1]]<-box_p
  
}

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Cancer.type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
#tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)
out_p2<-sum_table[which(sum_table$Cluster %in% c("Mono_C3-IL1B" ,"Mono_C2-RTEN" ,"TAM_C5-SPP1"  ) ),]

out_p<-rbind(out_p1,out_p2)
write.table(out_p,"source_data_NC/Figure_S7A.txt",sep = "\t",row.names = F,col.names = T)


box_p<-ggarrange(plotlist = ps_list,ncol=3,nrow=3,common.legend=T,legend = "none")


pdf("14.Figure/Figure_S5.pdf",width = 12,height = 8)
print(box_p)
dev.off()
