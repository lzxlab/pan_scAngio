library(tidyverse)
library(patchwork)
library(Seurat)

setwd('~/project/scPericyte')

# Figure S1 (A) ----
seu_merge <- readRDS("3.Cluster/13.Annotation/sc_merge.rds")

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

