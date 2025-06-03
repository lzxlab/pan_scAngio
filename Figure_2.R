library(tidyverse)
library(Seurat)

setwd('~/project/scPericyte')

# Figure 2 (A) ----
markers <- c("VEGFA","VEGFB","VEGFC","PGF","ANGPT2","HGF",
             "FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGF11","FGF12",
             "FGF13","FGF14","FGF16","FGF17","FGF18","FGF19","FGF20","FGF21","FGF22","FGF23","MET",
             "FLT1","KDR","FLT4","FGFR1","FGFR2","FGFR3","FGFR4")

seu_neg <- readRDS("3.Cluster/13.Annotation/non_immune.rds")
seu_pos <- readRDS("3.Cluster/13.Annotation/immune.rds")
seu_neutro <- readRDS("3.Cluster/8.SubAnnotation/6.Neutro/sub_sce_annotation.rds")

seu_merge <- merge(subset(seu_neg, features = markers),
                   y = c(subset(seu_pos, features = markers),
                         subset(seu_neutro, features = markers)))
Seurat::DefaultAssay(seu_merge) <- "RNA"
seu_merge <- Seurat::ScaleData(seu_merge)

p_dot <- Seurat::DotPlot(seu_merge, assay = "RNA", group.by = 'TopCluster', features = markers) &
  {Seurat::RotatedAxis() +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
      labs(x = NULL, y = NULL,
           size = "Percentage\nExpressed",
           color = "Average\nExpression")}

p_tree <- p_dot$data %>% 
  dplyr::select(id, features.plot, avg.exp.scaled) %>% 
  tidyr::pivot_wider(id_cols = id, names_from = features.plot, values_from = avg.exp.scaled) %>% 
  tibble::column_to_rownames(var = 'id') %>% 
  {hclust(d = dist(x = .))} %>%
  ggtree::ggtree(., branch.length = 'none')

Fig2A <- p_dot %>% aplot::insert_left(p_tree, width = 0.05)

# Figure 2 (B) ----
plot_ls <- list()

clusters_stroma <- c(
  "mCAF_C1-POSTN", "mCAF_C2-SFRP4", "iCAF_C1-IGF1", "iCAF_C2-MMP3", "iCAF_C3-CCL2", 
  "NF_C1-PI16", "NF_C2-SOD2", "NF_C3-CXCL12", "NF_C4-HBB", "NF_C5-IER", "NF_C6-CD74", 
  "PC_C1_mPC-ACTG2", "PC_C2_mPC-MYH11", "PC_C3_imPC-CD36", "PC_C4_imPC-MCAM", "PC_C5-prolif"
)
colors_stroma <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
                   "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
names(colors_stroma) <- clusters_stroma

seu_stroma <- readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")

plot_ls[['Stroma']] <- Seurat::DimPlot(object = seu_stroma, reduction = "umap", group.by = "SubCluster",
                                       cols = colors_stroma, label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )


clusters_epi <- c(
  "Epi_C1-TFF1", "Epi_C2-prolif", "Epi_C3-PIGR", "Epi_C4-REG1A", "Epi_C5-NDRG1", "Epi_C6-SPRR2A", 
  "Epi_C7-LY6D", "Epi_C8-SAA1", "Epi_C9-TNC", "Epi_C10-SST", "Epi_C11-ATF3", "Epi_C12-EMT", 
  "Epi_C13-IGHG1", "Epi_C14-ISG", "Epi_C15-A2M", "Epi_C16-BCAM"
)
colors_epi <- c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF",
                "#F1CE63FF", "#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF", 
                "#D37295FF", "#FABFD2FF")


names(colors_epi) <- clusters_epi

seu_epi <- readRDS("3.Cluster/5.SubAnnotation/5.Epi/sub_sce_annotation.rds")

plot_ls[['Epi']] <- Seurat::DimPlot(object = seu_epi, reduction = "umap", group.by = "SubCluster",
                                    cols = colors_epi, label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )


clusters_myeloid <- c(
  "DC_C1_pDC-LILRA4", "DC_C2_cDC1-CCL19", "DC_C3_cDC1-CPVL",  "DC_C4_cDC2-FCER1A",
  "Mono_C1-S100A9",  "Mono_C2-RTEN", "Mono_C3-IL1B", 
  "TAM_C1-MARCO", "TAM_C2-F13A1", "TAM_C3-TREM2", "TAM_C4-APOE", "TAM_C5-SPP1", "TAM_C6-CXCL9",
  "TAM_C7-FTL", "TAM_C8-IL32", "TAM_C9-prolif", "TAM_C10-MMP9" 
)
colors_myeloid <- c("#1F77B4FF", "#AEC7E8FF", "#FF7F0EFF", "#FFBB78FF", "#2CA02CFF", "#98DF8AFF", "#D62728FF", 
                    "#FF9896FF", "#9467BDFF", "#C5B0D5FF", "#8C564BFF", "#C49C94FF", "#E377C2FF", "#F7B6D2FF", 
                    "#7F7F7FFF", "#C7C7C7FF", "#BCBD22FF")

names(colors_myeloid) <- clusters_myeloid

seu_myeloid <- readRDS("3.Cluster/5.SubAnnotation/4.myeloid/sub_sce_annotation.rds")

plot_ls[['Myeloid']] <- Seurat::DimPlot(object = seu_myeloid, reduction = "umap", group.by = "SubCluster",
                                        cols = colors_myeloid, label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )


clusters_neutro <- c(
  "Neutro_C1-NIBAN1", "Neutro_C2-CD83", "Neutro_C3-TSPO", "Neutro_C4-CFD", "Neutro_C5-RESF1", "Neutro_C6-IFI6", 
  "Neutro_C7-THBS1", "Neutro_C8-PPDPF", "Neutro_C9-EGR1"
)
colors_neutro <- c("#729ECEFF", "#FF9E4AFF", "#67BF5CFF", "#ED665DFF", "#AD8BC9FF",
                   "#A8786EFF", "#ED97CAFF", "#A2A2A2FF", "#CDCC5DFF")

names(colors_neutro) <- clusters_neutro

seu_neutro <- readRDS("3.Cluster/8.SubAnnotation/6.Neutro/sub_sce_annotation.rds")

plot_ls[['Neutro']] <- Seurat::DimPlot(object = seu_neutro, reduction = "umap", group.by = "SubCluster",
                                       cols = colors_neutrod, label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

Fig2B <- cowplot::plot_grid(plot_ls, nrow = 1)

# Figure 2 (C) ----
Fig2C <- Seurat::DotPlot(seu_merge, assay = "RNA", group.by = 'SubCluster_S', features = markers) &
  {Seurat::RotatedAxis() +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
      labs(x = NULL, y = NULL,
           size = "Percentage\nExpressed",
           color = "Average\nExpression")}


