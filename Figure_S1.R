library(tidyverse)
library(patchwork)
library(Seurat)

setwd("~/data/single_cell/18.pan_Endo/PGF/")



# Figure S1 (A) ----

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



seu_merge<-merge(sub_sce1,sub_sce2)

p1 <- seu_merge@meta.data %>% 
  dplyr::distinct(Cancer.type, Type, Sample_ID) %>% 
  dplyr::group_by(Cancer.type) %>%
  dplyr::mutate(total = n()) %>%
  ungroup() %>%
  dplyr::mutate(Type = forcats::fct_relevel(Type, c("PN", "PT", "mLN", "mBrain")),
                Cancer.type = forcats::fct_reorder(Cancer.type, total, .desc = T)
  ) %>% 
  ggplot(aes(x = Cancer.type, fill = Type)) + 
  geom_bar(width = 0.6, position = "stack") +
  scale_fill_manual('', values = c('#D6604DFF', '#F4A582FF', '#878787FF', '#4D4D4DFF')) +
  labs(x = '', y = 'Number of patients', fill = 'Tissue') +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

total_da1<-table(p1$data$Cancer.type,p1$data$Type)
total_da1<-data.frame(Cancer.type=rownames(total_da1),total_da1)
write.table(total_da1,"source_data_NC/Figure_S1A_top.txt",sep = "\t",row.names = F,col.names = T)


p2 <- seu_merge@meta.data %>% 
  dplyr::group_by(Cancer.type) %>%
  dplyr::mutate(total = n()) %>%
  ungroup() %>%
  dplyr::mutate(Type = forcats::fct_relevel(Type, c("PN", "PT", "mLN", "mBrain")),
                Cancer.type = forcats::fct_reorder(Cancer.type, total, .desc = T)
  ) %>% 
  ggplot(aes(x = Cancer.type, fill = Type)) + 
  geom_bar(width = 0.6, position = "stack") +
  scale_y_continuous(labels = scales::label_number(scale = 1/10000)) +
  scale_fill_manual('', values = c('#D6604DFF', '#F4A582FF', '#878787FF', '#4D4D4DFF')) +
  labs(x = 'Data statistics (scRNA-seq)', y = expression('Number of cells (× 10'^4*')'), fill = 'Tissue') +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )


total_da1<-table(p2$data$Cancer.type,p2$data$Sample.type)
total_da1<-data.frame(Cancer.type=rownames(total_da1),total_da1)
write.table(total_da1,"source_data_NC/Figure_S1A_bottom.txt",sep = "\t",row.names = F,col.names = T)


FigS1A <- p1 / p2 + 
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(theme = theme(legend.position = "top", legend.box = "horizontal"))

# Figure S1 (B) ----
colors <- c("#A6CEE3FF", "#B2DF8AFF", "#FB9A99FF", "#FDBF6FFF", "#CAB2D6FF", "#8DD3C7FF",
            "#FFFFB3FF", "#BEBADAFF", "#FB8072FF", "#80B1D3FF", "#FDB462FF", "#B3DE69FF", "#FCCDE5FF")

seu_neg <- readRDS("3.Cluster/13.Annotation/non_immune.rds")
FigS1B <- Seurat::DimPlot(object = seu_neg, reduction = "tsne", group.by = "Cancer.type",
                         cols = colors, label = F, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

pos_data<-data.frame(Cancer.type=seu_neg$Cancer.type,round(seu_neg@reductions$tsne@cell.embeddings,2))
write.table(pos_data,"source_data_NC/Figure_S1B.txt",sep = "\t",row.names = F,col.names = T)


# Figure S1 (C) ----
seu_pos <- readRDS("3.Cluster/13.Annotation/immune.rds")
FigS1C <- Seurat::DimPlot(object = seu_pos, reduction = "tsne", group.by = "Cancer.type",
                          cols = colors, label = F, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

pos_data<-data.frame(Cancer.type=seu_pos$Cancer.type,round(seu_pos@reductions$tsne@cell.embeddings,2))
write.table(pos_data,"source_data_NC/Figure_S1C.txt",sep = "\t",row.names = F,col.names = T)



# Figure S1 (D) ----
FigS1D_up <- Seurat::FeaturePlot(object = seu_neg, reduction = "tsne", ncol = 5,
                                 features = c('PTPRC','DCN','FAP','RGS5','PECAM1','EPCAM'),
                                 min.cutoff = 0, max.cutoff = 3, order = T, raster = T) &
  {scale_color_gradientn(colours = c('grey90', '#a13037')) +
      labs(x = NULL, y = NULL, color = 'Expression') +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )}

FigS1D_up <- Seurat::FeaturePlot(object = seu_neg, reduction = "tsne", ncol = 5,
                                 features = c('PTPRC','DCN','FAP','RGS5','PECAM1','EPCAM'),
                                 min.cutoff = 0, max.cutoff = 3, order = T, raster = T) 
data1<-data.frame(Cluster=FigS1D_up[[1]]$data$ident,
                  tSNE_1=round(FigS1D_up[[1]]$data$tSNE_1,2),
                  tSNE_2=round(FigS1D_up[[1]]$data$tSNE_2,2),
                  PTPRC=round(FigS1D_up[[1]]$data$PTPRC,2),
                  DCN=round(FigS1D_up[[2]]$data$DCN,2),
                  FAP=round(FigS1D_up[[3]]$data$FAP,2),
                  RGS5=round(FigS1D_up[[4]]$data$RGS5,2),
                  PECAM1=round(FigS1D_up[[5]]$data$PECAM1,2),
                  EPCAM=round(FigS1D_up[[6]]$data$EPCAM,2)
)
write.table(data1,"source_data_NC/Figure_S1D_top.txt",sep = "\t",row.names = F,col.names = T)



FigS1D_down <- Seurat::FeaturePlot(object = seu_pos, reduction = "tsne", ncol = 5,
                                   features = c('PTPRC','CD3D','CD4','FOXP3','CD8A','NKG7',
                                                'CD79A','MZB1','CD1C','FCN1','C1QC','CPA3'),
                                   min.cutoff = 0, max.cutoff = 3, order = T, raster = T) &
  {scale_color_gradientn(colours = c('grey90', '#a13037')) +
      labs(x = NULL, y = NULL, color = 'Expression') +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )}

FigS1D_down <- Seurat::FeaturePlot(object = seu_pos, reduction = "tsne", ncol = 5,
                                   features = c('PTPRC','CD3D','CD4','FOXP3','CD8A','NKG7',
                                                'CD79A','MZB1','CD1C','FCN1','C1QC','CPA3'),
                                   min.cutoff = 0, max.cutoff = 3, order = T, raster = T)
data1<-data.frame(Cluster=FigS1D_down[[1]]$data$ident,
                  tSNE_1=round(FigS1D_down[[1]]$data$tSNE_1,2),
                  tSNE_2=round(FigS1D_down[[1]]$data$tSNE_2,2),
                  PTPRC=round(FigS1D_down[[1]]$data$rna_PTPRC,2),
                  CD3D=round(FigS1D_down[[2]]$data$rna_CD3D,2),
                  CD4=round(FigS1D_down[[3]]$data$CD4,2),
                  FOXP3=round(FigS1D_down[[4]]$data$FOXP3,2),
                  CD8A=round(FigS1D_down[[5]]$data$CD8A,2),
                  NKG7=round(FigS1D_down[[6]]$data$NKG7,2),
                  
                  CD79A=round(FigS1D_down[[7]]$data$CD79A,2),
                  MZB1=round(FigS1D_down[[8]]$data$MZB1,2),
                  CD1C=round(FigS1D_down[[9]]$data$CD1C,2),
                  FCN1=round(FigS1D_down[[10]]$data$FCN1,2),
                  C1QC=round(FigS1D_down[[11]]$data$C1QC,2),
                  CPA3=round(FigS1D_down[[12]]$data$CPA3,2)
                  
)
write.table(data1,"source_data_NC/Figure_S1D_bottom.txt",sep = "\t",row.names = F,col.names = T)

