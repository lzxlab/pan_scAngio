library(tidyverse)
library(Seurat)

setwd('~/project/scPericyte')

# Figure 1 (B) ----
top_clusters <- c("NF","CAF","PC","Endo","Epi","CD4_Tconv","CD4_Treg","CD8","MAIT","NK","B","Plasma","DC","Mono","Macro","Mast")
top_colors <- c("#A6CEE3FF", "#B2DF8AFF", "#FB9A99FF", "#FDBF6FFF", "#CAB2D6FF", "#8DD3C7FF",
            "#FFFFB3FF", "#BEBADAFF", "#FB8072FF", "#80B1D3FF", "#FDB462FF", "#B3DE69FF",
            "#FCCDE5FF", "#BC80BDFF", "#CCEBC5FF", "#FFED6FFF")
names(top_colors) <- top_clusters

seu_neg <- readRDS("3.Cluster/13.Annotation/non_immune.rds")
Fig1B <- Seurat::DimPlot(object = seu_neg, reduction = "tsne", group.by = "TopCluster",
                         cols = top_colors[1:5], label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

# Figure 1 (C) ----
seu_pos <- readRDS("3.Cluster/13.Annotation/immune.rds")
Fig1C <- Seurat::DimPlot(object = seu_pos, reduction = "tsne", group.by = "TopCluster",
                         cols = top_colors[6:16], label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

# Figure 1 (D) ----
top_markers <- c("APOD", "CFD", "DCN", "FAP", "COL10A1", "COL11A1", "RGS5",
                 "ACTA2", "MYH11", "VWF", "PECAM1", "CDH5", "EPCAM", "KRT18",
                 "KRT17", "PTPRC", "CD3D", "CD3G", "CD3E", "CD4", "FOXP3", "CD8A",
                 "CD8B", "SLC4A10", "TRAV1-2", "NCR3", "NKG7", "GNLY", "PRF1", 
                 "CD79A", "CD79B", "CD19", "MS4A1", "MZB1", "JCHAIN", "IGHA1", 
                 "IGHG1", "CD1C", "HLA-DRA", "HLA-DPB1", "FCN1", "CD14", "FCGR3A", 
                 "CD68", "CD163", "APOC1", "C1QC", "CPA3", "TPSB2", "KIT")

seu_merge <- merge(subset(seu_neg, features = top_markers),
                   subset(seu_pos, features = top_markers))
Seurat::DefaultAssay(seu_merge) <- "RNA"
seu_merge <- Seurat::ScaleData(seu_merge)

Fig1D <- Seurat::DotPlot(seu_merge, assay = "RNA", group.by = 'TopCluster', features = top_markers) &
  {Seurat::RotatedAxis() +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
      labs(x = NULL, y = NULL,
           size = "Percentage\nExpressed",
           color = "Average\nExpression")}

# Figure 1 (E) ----
endo_clusters <- c("Artery_EC","Capillary_EC","Venous_EC","Tip_EC","Immature_EC1","Immature_EC2","Other_EC","LEC")
endo_colors <- c("#8DD3C7FF", "#FFFFB3FF", "#BEBADAFF", "#FB8072FF", "#80B1D3FF", "#FDB462FF", "#B3DE69FF", "#FCCDE5FF")
names(endo_colors) <- endo_clusters

seu_endo <- readRDS("3.Cluster/13.Annotation/1.Endo_annotation.rds")
Fig1E <- Seurat::DimPlot(object = seu_endo, reduction = "tsne", group.by = "MidCluster",
                         cols = endo_colors, label = T, raster = T) &
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

# Figure 1 (F) ----
fig_ls <- list()
for (TYPE in c("PN", "PT", "mBrain", "mLN")) {
  seu_temp <- subset(seu_endo, Type == TYPE)
  fig_ls[[TYPE]] <- Seurat::DimPlot(object = seu_temp, reduction = "tsne", group.by = "MidCluster",
                                    cols = endo_colors, label = T, raster = T) &
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
}
Fig1F <- cowplot::plot_grid(fig_ls, nrow = 2)

# Figure 1 (G) ----
endo_markers <- c("CXCL12", "GJA4", "CD36", "CA4", "ACKR1", "SELE", "KDR", 
                  "CXCR4", "APLNR", "SAT1", "KRT17", "ISG15", "PROX1", "CCL21")

Fig1G <- Seurat::DotPlot(seu_endo, assay = "RNA", group.by = 'MidCluster', features = endo_markers) &
  {coord_flip() +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
      labs(x = NULL, y = NULL,
           size = "Percentage\nExpressed",
           color = "Average\nExpression") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))}

# Figure 1 (H)
milo.obj <- readRDS("4.characteristics/7.1.Milo_Endo.milo.obj.RDS")
milo.res <- readRDS("4.characteristics/7.1.Milo_Endo.milo.res.RDS")

Fig1H <- miloR::plotNhoodGraphDA(milo.obj, milo.res, layout = "tsne", alpha = 0.1) +
  scale_fill_gradient2(low = "#3A3A98", mid = "white", high = "#832424", midpoint = 0, name = "log2FC") +
  scale_size(range = c(0.5,4), name = "Nhood size") +
  ggraph::scale_edge_width(range = c(0.2,3), name = "overlap size") +
  scale_x_continuous(name = 'tSNE_1', breaks = NULL) +
  scale_y_continuous(name = 'tSNE_2', breaks = NULL) +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

# Figure 1 (I) ----
Fig1I <- milo.res %>% 
  dplyr::mutate(MidCluster = factor(MidCluster, levels = names(sort(tapply(logFC, MidCluster, median))))) %>% 
  droplevels() %>% 
  {miloR::plotDAbeeswarm(., group.by = "MidCluster", alpha = 1) + 
      geom_jitter(aes(group = MidCluster), size = 1.2, alpha = 0.8, shape = 21) +
      scale_color_gradient2(low = "#3A3A98", mid = "white", high = "#832424", limits = c(-5, 5), oob = scales::squish) +
      geom_boxplot(aes(group = MidCluster), outlier.shape = NA, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic(base_size = 14) + 
      labs(y = "Log2 Fold Change", x = NULL)}

# Figure 1 (J) ----
OR_mtx <- readRDS("4.characteristics/8.1.OR_Endo.or.RDS")
p_mtx <- readRDS("4.characteristics/8.1.OR_Endo.p.RDS")

Fig1J <- ComplexHeatmap::Heatmap(
  matrix = OR_mtx, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (p_mtx[i,j] < 0.001) {
      grid.text("***", x = x, y = y, gp = gpar(col = "black"))
    } else if (p_mtx[i,j] < 0.01) {
      grid.text("**", x = x, y = y, gp = gpar(col = "black"))
    } else if (p_mtx[i,j] < 0.05) {
      grid.text("*", x = x, y = y, gp = gpar(col = "black"))
    } else {
      grid.text("", x = x, y = y, gp = gpar(col = "black"))
    }
  },
  col = c("#7171C5", "white", "#BC5D5D"),
  cluster_rows = TRUE, cluster_columns = TRUE,
  column_names_rot = 45, 
  column_title = "", 
  column_title_gp = gpar(fontface = "bold"),
  width = unit(nrow(OR_mtx)/1.25, "cm"), 
  height = unit(ncol(OR_mtx)/1.5, "cm"), 
  show_heatmap_legend = TRUE, 
  heatmap_legend_param = list(
    title = "Odds Ratio",
    labels_gp = gpar(fontface = "plain")
  )
)
