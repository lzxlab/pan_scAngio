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







sub_sce<-readRDS("3.Cluster/5.SubAnnotation/4.myeloid/sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("PN" ,   "PT", "mBrain", "mLN" ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("PN" ,   "PT", "mBrain", "mLN" ))
#dir.create("3.Cluster/13.plot/3.myeloid/")
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=rev(c(       "DC_C1_pDC-LILRA4", "DC_C2_cDC1-CCL19" , "DC_C3_cDC1-CPVL" ,  "DC_C4_cDC2-FCER1A",
                                                                  "Mono_C1-S100A9"  ,  "Mono_C2-RTEN"    ,  "Mono_C3-IL1B"   ,   "TAM_C1-MARCO"  ,   
                                                                  "TAM_C2-F13A1"    ,  "TAM_C3-TREM2"   ,   "TAM_C4-APOE"   ,   
                                                                  "TAM_C5-SPP1"    ,   "TAM_C6-CXCL9"    ,  "TAM_C7-FTL"     ,   "TAM_C8-IL32"  ,    
                                                                  "TAM_C9-prolif"  , "TAM_C10-MMP9")))       
markers1<-c("LILRA4","IRF7","PTGDS","CCL19","CCL22","LAMP3","CPVL","CLEC9A","HSPD1",
            "S100A8","S100A9","VCAN","LYZ","RETN","EREG","IL1B","CXCL8","CCL20",
            "FABP4","MARCO","LGALS3",
            "F13A1","FOLR2","STAB1",
            "TREM2","CCL4L2","CCL3L2",
            "CTSS","APOE","PLTP",
            "SPP1","CSTB","ADM",
            "ISG15","ISG20","CXCL10",
            "FTL","HBB","CLEC5A",
            "IL32","CCL5","CD2",
            "STMN1","TUBB","HMGB2",
            "MMP9","CTSK","CKB")

library(viridis)
Idents(sub_sce)<-sub_sce$SubCluster
Myeloid_dot_p1<-DotPlot(sub_sce, features = unique(markers1),assay = "RNA",
                cols = c("grey","blue"))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

Myeloid_dot_p1




endo_data<-Myeloid_dot_p1$data
write.table(endo_data,"source_data_NC/Figure_S4C.txt",sep = "\t",row.names = F,col.names = T)




sub_sce<-readRDS("3.Cluster/5.SubAnnotation/5.Epi//sub_sce_annotation.rds")
markers<-unique(c( "TFF1"  ,  "PIP"    , "CXCL14" , "UBE2C"  , "TUBA1B" , "CDC20"  , "PIGR", "SLPI"  ,  "WFDC2"  , "REG1A",   "MT1G" ,  
                   "REG3A" ,  "S100A2" , "NDRG1"  , "SLC2A1" , "SPRR2E" , "SPRR2A",  "SPRR3"  , "FABP5" ,  "KRT16"  , "LY6D"  ,  "SAA1" ,  
                   "MMP7"  ,  "LTF"    , "TNC"  , "TAGLN" ,  "LAMC2"  , "SST"   ,  "CPB1"   , "GAST"  ,  "HSPA1B" , "ATF3"  ,  "FOS" ,   
                   "SPARC" ,  "VIM"    , "COL1A1" , "IGLC2"  , "SRGN"  ,  "IGKC"  ,  "CXCL10"  ,"ISG15" ,  "IFITM1" , "A2M"   ,  "MIA" ,   
                   "SFRP1" ,  "BCAM"   , "GPC3"   , "KRT15" ))
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=rev(c(    "Epi_C1-TFF1" , "Epi_C2-prolif", "Epi_C3-PIGR"  , "Epi_C4-REG1A" , "Epi_C5-NDRG1" , "Epi_C6-SPRR2A" ,"Epi_C7-LY6D" , 
                                                               "Epi_C8-SAA1" ,  "Epi_C9-TNC"  , "Epi_C10-SST"  , "Epi_C11-ATF3" , "Epi_C12-EMT" ,  "Epi_C13-IGHG1", "Epi_C14-ISG" ,  "Epi_C15-A2M" , 
                                                               "Epi_C16-BCAM" )))     
Idents(sub_sce)<-sub_sce$SubCluster
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

library(viridis)
Idents(sub_sce)<-sub_sce$SubCluster
Epi_dot_p1<-DotPlot(sub_sce, features = unique(markers),assay = "RNA",
                cols = c("grey","blue"))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

Epi_dot_p1


endo_data<-Epi_dot_p1$data
write.table(endo_data,"source_data_NC/Figure_S4A.txt",sep = "\t",row.names = F,col.names = T)









sub_sce<-readRDS("3.Cluster/8.SubAnnotation/6.Neutro/sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("N" ,   "T", "PB"  ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("N" ,   "T", "PB"  ))
#dir.create("3.Cluster/13.plot/3.myeloid/")
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=rev(c(   "Neutro_C1-NIBAN1", "Neutro_C2-CD83"  , "Neutro_C3-TSPO" ,  "Neutro_C4-CFD" ,  
                                                              "Neutro_C5-RESF1" , "Neutro_C6-IFI6"  , "Neutro_C7-THBS1" , "Neutro_C8-PPDPF", 
                                                              "Neutro_C9-EGR1")))       
markers1<-c("NIBAN1"  ,  "GNG10"    , "CLEC12A" ,  "CD83"    ,  "PLIN2"  ,   "CCL3L1"  , 
            "TSPO" ,  "TXN"     ,  "CST7"   ,   "CFD"    ,   "FAM129A",   "KIAA1551", 
            "RESF1"   ,  "PCBP1-AS1", "MT-ND2"  ,  "IFI6"  ,   "MX1"    ,   "IFIT1"  ,  
            "HLA-DPA1" , "THBS1"    ,  "CCL18"  ,   "PPDPF"  ,   "SFTPB"  ,   "IGFBP7" ,  
            "FOS"     ,  "EGR1"    ,  "NFKBIZ")

library(viridis)
Idents(sub_sce)<-sub_sce$SubCluster
Neutro_dot_p1<-DotPlot(sub_sce, features = unique(markers1),assay = "RNA",
                cols = c("grey","blue"))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

Neutro_dot_p1

endo_data<-Neutro_dot_p1$data
write.table(endo_data,"source_data_NC/Figure_S4D.txt",sep = "\t",row.names = F,col.names = T)






sub_sce<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
markers<-unique(c("POSTN","COL1A1","FAP",
                  "SFRP4","CXCL14","ELN",
                  "APOD","IGF1","THBS4",
                  "MMP3","CXCL8","IL24",
                  "CCL2","CCL11","PTGDS",
                  "PI16","MFAP5","CLU",
                  "SOD2","CXCL2","FOSL1",
                  "CXCL12","MFAP4","ADH1B",
                  "HBB",
                  "IER2","EGR1","FOS",
                  "CD74","HLA-DRA","HLA-DRB1",
                  "ACTG2","DES","MYLK",
                  "MYH11","ADIRF","ACTA2",
                  "RGS5","CCL19","CD36",
                  "MCAM","THY1","MKI67" ))
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=rev(c(   "mCAF_C1-POSTN"  ,   "mCAF_C2-SFRP4",  "iCAF_C1-IGF1"    ,  "iCAF_C2-MMP3"     , "iCAF_C3-CCL2"    ,     "NF_C1-PI16"     ,  
                                                              "NF_C2-SOD2"       , "NF_C3-CXCL12"     , "NF_C4-HBB"       ,  "NF_C5-IER"      ,   "NF_C6-CD74"   ,     "PVL_C1_mPC-ACTG2" ,
                                                              "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36"  ,"PVL_C4_imPC-MCAM", "PVL_C5-prolif")))     
Idents(sub_sce)<-sub_sce$SubCluster
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

library(viridis)
Idents(sub_sce)<-sub_sce$SubCluster
Stromal_dot_p1<-DotPlot(sub_sce, features = unique(markers),assay = "RNA",
                cols = c("grey","blue"))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

Stromal_dot_p1

endo_data<-Stromal_dot_p1$data
endo_data$id <-gsub("PVL","PC",endo_data$id)
write.table(endo_data,"source_data_NC/Figure_S4B.txt",sep = "\t",row.names = F,col.names = T)


total_p1<-Stromal_dot_p1+Epi_dot_p1+Myeloid_dot_p1+Neutro_dot_p1+plot_layout(ncol = 1,nrow = 4,
                                                                             heights = c(4,4,4,3))

pdf("14.Figure/Figure_S4_markers.pdf",width = 12,height = 15)
print(total_p1)
dev.off()

