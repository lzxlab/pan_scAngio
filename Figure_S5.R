
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




##Expresion compare CAF
combined1<-readRDS("3.Cluster/13.Annotation/non_immune.rds")


Idents(combined1)<-combined1$MidCluster
combined<-subset(combined1,idents=c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                   "Other_ECs" ,    "PC" ,          "Tip_ECs"     , 
                                   "Venous_ECs"  ))


combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster[which(combined$MidCluster %in% c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                                     "Other_ECs" ,        "Tip_ECs"     , 
                                                     "Venous_ECs"  ))]<-combined$MidCluster[which(combined$MidCluster %in% c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                                                                                                            "Other_ECs" ,        "Tip_ECs"     , 
                                                                                                                            "Venous_ECs"  ))]

library(limma)

combined$SubCluster<-factor(combined$SubCluster,levels = c(
  "Artery_ECs"     ,"Capillary_ECs"   , "Immature_ECs",  
  "Venous_ECs"  ,   "Tip_ECs" ,    "Other_ECs"   , 
  "PC_C1_mPC-ACTG2" ,"PC_C2_mPC-MYH11",
  "PC_C3_imPC-CD36" ,"PC_C4_imPC-MCAM" ,"PC_C5-prolif"    
  
))



Cluster<-levels(combined$SubCluster)
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster



box_p1<-VlnPlot(combined,features = c("ANGPT2","PGF","MCAM"),
                group.by = "TopCluster",cols = RColorBrewer::brewer.pal(3,"Set1"),
                pt.size = 0) +
  plot_layout(ncol=3) &
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y.npc = 0.9) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p1

endo_data<-round(t(data.frame(combined@assays$integrated@scale.data[c("ANGPT2","PGF","MCAM"),])),2)
endo_data<-data.frame(Cluster=combined$TopCluster,endo_data)

write.table(endo_data,"source_data_NC/Figure_S5C.txt",sep = "\t",row.names = F,col.names = T)



combined$MCAM<-combined@assays$integrated@scale.data["MCAM",]
combined$MCAM_level<-"MCAM-"
combined$MCAM_level[which(combined$MCAM>1)]<-"MCAM+"

combined$Type<-paste0(combined$MCAM_level,combined$TopCluster)

box_p2<-VlnPlot(combined,features = c("ANGPT2","PGF"),
                group.by = "MCAM_level",
                cols = RColorBrewer::brewer.pal(3,"Set1"),
                split.by = "TopCluster",
                pt.size = 0) +
  plot_layout(ncol=2) &
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y.npc = 0.9) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p2

endo_data<-round(t(data.frame(combined@assays$integrated@scale.data[c("ANGPT2","PGF"),])),2)
endo_data<-data.frame(Cluster=combined$TopCluster,MCAM_level=combined$MCAM_level,endo_data)

write.table(endo_data,"source_data_NC/Figure_S5D.txt",sep = "\t",row.names = F,col.names = T)


Idents(combined1)<-combined1$MidCluster
combined<-subset(combined1,idents=c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                    "Other_ECs" ,    "PC" ,          "Tip_ECs"     , 
                                    "Venous_ECs"  ))


combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster[which(combined$MidCluster %in% c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                                     "Other_ECs" ,        "Tip_ECs"     , 
                                                     "Venous_ECs"  ))]<-combined$MidCluster[which(combined$MidCluster %in% c("Artery_ECs"   , "Capillary_ECs", "Immature_ECs"  ,
                                                                                                                             "Other_ECs" ,        "Tip_ECs"     , 
                                                                                                                             "Venous_ECs"  ))]

library(limma)

combined$SubCluster<-factor(combined$SubCluster,levels = c(
  "Artery_ECs"     ,"Capillary_ECs"   , "Immature_ECs",  
  "Venous_ECs"  ,   "Tip_ECs" ,    "Other_ECs"   , 
  "PC_C1_mPC-ACTG2" ,"PC_C2_mPC-MYH11",
  "PC_C3_imPC-CD36" ,"PC_C4_imPC-MCAM" ,"PC_C5-prolif"    
  
))


combined$Type<-as.character(combined$Type)
combined$Type[which(combined$Type!="PN")]<-"Tumor"
combined$Type[which(combined$Type=="PN")]<-"Normal"


Cluster<-levels(combined$SubCluster)
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster


box_p4<-VlnPlot(combined,features = c("ANGPT2","PGF","MCAM"),
                group.by = "TopCluster",
                cols = c( "#4DBBD5","#E64B35"),
                split.by = "Type",
                pt.size = 0) +
  plot_layout(ncol=3) &
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y.npc = 0.9) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p4


endo_data<-round(t(data.frame(combined@assays$integrated@scale.data[c("ANGPT2","PGF","MCAM"),])),2)
endo_data<-data.frame(Cluster=combined$TopCluster,Type=combined$Type,endo_data)

write.table(endo_data,"source_data_NC/Figure_S5E.txt",sep = "\t",row.names = F,col.names = T)



total_p1<-ggarrange(ggarrange(box_p1,box_p2,widths = c(1,1.3)),box_p4,ncol = 1,nrow = 3)
pdf("14.Figure/3.Figure_3_S3_R1.pdf",width = 11,height = 8)
print(total_p1)
dev.off()
