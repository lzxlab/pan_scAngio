
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
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")



combined1<-readRDS("3.Cluster/5.SubAnnotation/5.Epi/sub_sce_annotation.rds")
combined<-combined1

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

#combined$Cluster<-factor(combined$Cluster,levels = rev(c( "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
#                                                          "DC" , "Mono"   , "Macro" ,"Mast" ,"pAD", "CAF"   ,   "PVL",     "Endo" ,   "Epi")))
markers1<-c("VEGFA","VEGFB","VEGFC","PGF","ANGPT2","HGF","FGF1" , "FGF2" ,"FGF3" , "FGF4",  "FGF5",  "FGF6" , "FGF7" , "FGF8" , "FGF9" , "FGF11", "FGF12" ,"FGF13" ,"FGF14" ,"FGF16", "FGF17", "FGF18" ,"FGF19"  , "FGF20", "FGF21",
            "FGF22", "FGF23" ,"MET")
markers2<-c("FLT1","KDR","FLT4","FGFR1","FGFR2","FGFR3","FGFR4")
markers<-list("Angiogenic factors"=markers1,
              "Receptors"=markers2)
Idents(combined)<-combined$Cancer.type
Epi_p<-DotPlot(combined, features = markers,
            cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Expression of angiogenic factors and receptors in epithelial cells")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(axis.title=element_blank(),
        plot.title=element_text(size=10),
        legend.position = "right",
        panel.border = element_rect(size = 0.5,colour = "grey",fill = NA),
        axis.line = element_blank())
Epi_p


endo_data<-Epi_p$data
endo_data$id <-gsub("PVL","PC",endo_data$id)
write.table(endo_data,"source_data_NC/Figure_S6A.txt",sep = "\t",row.names = F,col.names = T)



combined1<-readRDS("3.Cluster/5.SubAnnotation/4.myeloid/sub_sce_annotation.rds")
combined<-subset(combined1,subset=MidCluster %in% c("Mono","Macro"))

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)


Idents(combined)<-combined$Cancer.type
Mye_p<-DotPlot(combined, features = markers,
            cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Expression of angiogenic factors and receptors in monecytes/TAM cells")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(axis.title=element_blank(),
        plot.title=element_text(size=10),
        legend.position = "right",
        panel.border = element_rect(size = 0.5,colour = "grey",fill = NA),
        axis.line = element_blank())


endo_data<-Mye_p$data
endo_data$id <-gsub("PVL","PC",endo_data$id)
write.table(endo_data,"source_data_NC/Figure_S6C.txt",sep = "\t",row.names = F,col.names = T)




combined1<-readRDS("3.Cluster/5.SubAnnotation/2.Stromal/sub_sce_annotation.rds")
combined<-subset(combined1,subset=TopCluster=="CAF")

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)



Idents(combined)<-combined$Cancer.type
Strmal_p<-DotPlot(combined, features = markers,
            cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Expression of angiogenic factors and receptors in CAF")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(axis.title=element_blank(),
        plot.title=element_text(size=10),
        legend.position = "right",
        panel.border = element_rect(size = 0.5,colour = "grey",fill = NA),
        axis.line = element_blank())


endo_data<-Strmal_p$data
endo_data$id <-gsub("PVL","PC",endo_data$id)
write.table(endo_data,"source_data_NC/Figure_S6B.txt",sep = "\t",row.names = F,col.names = T)



combined1<-readRDS("3.Cluster/13.Annotation/3.PC_annotation_new.rds")
combined<-subset(combined1,subset=TopCluster=="PC")

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Idents(combined)<-combined$Cancer.type
PVL_p<-DotPlot(combined, features = markers,
            cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Expression of angiogenic factors and receptors in pericytes")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(axis.title=element_blank(),
        plot.title=element_text(size=10),
        legend.position = "right",
        panel.border = element_rect(size = 0.5,colour = "grey",fill = NA),
        axis.line = element_blank())

endo_data<-PVL_p$data
endo_data$id <-gsub("PVL","PC",endo_data$id)
write.table(endo_data,"source_data_NC/Figure_S6D.txt",sep = "\t",row.names = F,col.names = T)



combined1<-readRDS("3.Cluster/13.Annotation/Endo_annotation.rds")
combined<-combined1

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Idents(combined)<-combined$Cancer.type
Endo_p<-DotPlot(combined, features = markers,
            cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Expression of angiogenic factors and receptors in endothelial cells")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(axis.title=element_blank(),
        plot.title=element_text(size=10),
        legend.position = "right",
        panel.border = element_rect(size = 0.5,colour = "grey",fill = NA),
        axis.line = element_blank())

library(patchwork)
total_p<-Epi_p+Strmal_p+Mye_p+PVL_p+plot_layout(ncol = 1,nrow = 4)
pdf("14.Figure/Figure_S4_angio.pdf",width = 12,height = 13.5)
print(total_p)
dev.off()