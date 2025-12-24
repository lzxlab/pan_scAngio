library(tidyverse)
library(patchwork)
library(Seurat)

setwd("~/data/single_cell/18.pan_Endo/PGF/")



# Figure S3 (A) ----
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




Figure_S3A <- seu_merge@meta.data %>% 
  dplyr::filter(Type %in% c('PT','PN')) %>%
  dplyr::filter(Cancer.type %in% c('COAD','LUAD')) %>% 
  dplyr::filter(Cluster != 'Epi') %>% 
  dplyr::mutate(Sample_ID2 = paste(Cancer.type, Patient_ID, Type, sep = '-')) %>% 
  dplyr::group_by(Cancer.type, Patient_ID) %>%
  dplyr::filter(dplyr::n_distinct(Type) == 2) %>% 
  dplyr::ungroup() %>% 
  droplevels() %>% 
  dplyr::count(Cancer.type, Patient_ID, Type, Cluster, name = "n_cells") %>%
  tidyr::complete(Cancer.type, Cluster, fill = list(n_cells = 0)) %>% 
  dplyr::mutate(n_cells = as.numeric(n_cells)) %>% 
  dplyr::group_by(Cancer.type, Patient_ID, Type) %>%
  dplyr::mutate(prop_cells = n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(Cluster == 'Endo') %>%   
  droplevels() %>% 
  dplyr::select(Cancer.type, Patient_ID, Type, prop_cells) %>% 
  tidyr::pivot_wider(id_cols = c(Cancer.type, Patient_ID), names_from = 'Type', values_from = 'prop_cells') %>% 
  ggplot(aes(x = PN, y = PT)) +
  facet_wrap(. ~ Cancer.type, scales = 'free') +
  geom_point(color = 'grey') +
  geom_smooth(method = "lm", se = F, color = "red") +
  ggpubr::stat_cor(label.x.npc = "left", label.y.npc = 1, vjust = -0.1, method = "spearman") +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = '% among non-Epi in adjacent non-tumor tissues', y = '% among non-Epi\nin tumor tissues') +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    axis.text = element_text(color = 'black')
  )

data1<-Figure_S3A$data
write.table(data1,"source_data_NC/Figure_S3A.txt",sep = "\t",row.names = F,col.names = T)


# Figure S3 (B) ----
FigS3B <- seu_merge@meta.data %>% 
  dplyr::mutate(Type = Type %>% forcats::fct_relevel('PN','PT','mBrain','mLN')) %>% 
  dplyr::filter(Cluster != 'Epi') %>% 
  droplevels() %>% 
  dplyr::mutate(Sample_ID = paste(Cancer.type, Sample_ID, Type, sep = '-')) %>%
  dplyr::count(Cancer.type, Type, Sample_ID, Cluster, name = "n_cells") %>%
  tidyr::complete(Cancer.type, Type, Sample_ID, Cluster, fill = list(n_cells = 0)) %>% 
  dplyr::group_by(Cancer.type, Type, Sample_ID) %>%
  dplyr::filter(sum(n_cells) > 0) %>%
  dplyr::mutate(prop_cells = n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(Cluster == 'Endo') %T>% 
  {wilcox_result <<- stats::wilcox.test(
    prop_cells ~ Type, 
    data = dplyr::filter(., Cluster == 'Endo', Type %in% c('PN', 'PT')))
  } %>% 
  ggplot(aes(x = Type, y = prop_cells, color = Type))+
  geom_boxplot(width = 0.8, position = position_dodge(0), fill = 'transparent', outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.5, alpha = 0.5, width = 0.4, method = "quasirandom") +
  annotate("segment", x = 1, xend = 2, y = 0.95, yend = 0.95, colour = "black") +
  annotate("text", x = 1.5, y = 1, label = sprintf("italic(p)==%s", signif(wilcox_result$p.value, 2)), parse = TRUE) +
  scale_x_discrete(labels = c("Primary\nNormal", "Primary\nTumor", "Lymph Node\nMetastasis", "Brain\nMetastasis")) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = expansion(mult = c(0.1, 0.15))) +
  scale_color_manual(values = c("#BD0026","#FEB24C","#02818A","#67A9CF"), name = "Tissue") +
  labs(x = NULL, y = '% among non-Epi', title = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(lineheight = 0.8, angle = 45, vjust = 1, hjust = 1),
        legend.position = 'none',
        strip.background = element_blank())

data1<-FigS3B$data
write.table(data1,"source_data_NC/Figure_S3B.txt",sep = "\t",row.names = F,col.names = T)



# Figure S3 (C) ----
cancertype_colors <- c("#D0D1E6", "#A6BDDB", "#67A9CF", "#3690C0", "#02818A", "#016C59", "#014636", 
                       "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")

FigS3C <- seu_merge@meta.data %>% 
  dplyr::filter(Type == 'PT') %>%
  dplyr::filter(Cluster != 'Epi') %>% 
  droplevels() %>% 
  dplyr::mutate(Sample_ID = paste(Cancer.type, Sample_ID, sep = '-')) %>%
  dplyr::count(Cancer.type, Sample_ID, Cluster, name = "n_cells") %>%
  tidyr::complete(Cancer.type, Sample_ID, Cluster, fill = list(n_cells = 0)) %>% 
  dplyr::group_by(Cancer.type, Sample_ID) %>%
  dplyr::filter(sum(n_cells) > 0) %>%
  dplyr::mutate(prop_cells = n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(Cluster == 'Endo') %T>% 
  {anova_result <<- aov(prop_cells ~ Cancer.type, data = .)
  f_value <<- summary(anova_result)[[1]]$`F value`[1]
  p_value <<- summary(anova_result)[[1]]$`Pr(>F)`[1]
  } %>% 
  {ggplot(., aes(x = Cancer.type, y = prop_cells, color = Cancer.type))+
      geom_boxplot(width = 0.8, position = position_dodge(0), fill = 'transparent', outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(size = 0.8, alpha = 0.9, width = 0.4, method = "quasirandom") +
      annotate("text", x = 1, y = max(.$prop_cells) * 0.95, parse = TRUE, hjust = 0, 
               label = sprintf("italic(p)==%s*','~italic(F)==%.2f", signif(p_value, 3), f_value)) +
      scale_y_continuous(labels = scales::percent_format(scale = 100)) +
      scale_color_manual(values = cancertype_colors, name = "Cancer Type") +
      labs(x = NULL, y = '% among non-Epi\nin tumor tissues', title = 'All Endothelial cells') +
      theme_classic(base_size = 14) +
      theme(axis.text = element_text(color = 'black'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.title = element_text(vjust = 0.5, hjust = 0.5),
            legend.position = 'none',
            strip.background = element_blank())}

data1<-FigS3C$data
write.table(data1,"source_data_NC/Figure_S3C.txt",sep = "\t",row.names = F,col.names = T)


seu_endo<-readRDS("3.Cluster/13.Annotation/1.Endo_annotation_new.rds")
# Figure S2 (F-K) ----

color_list<-c(RColorBrewer::brewer.pal(9,"PuBuGn")[3:9],RColorBrewer::brewer.pal(9,"YlOrRd")[3:9])

plot_ls <- list()

for (i in c("Artery_EC","Capillary_EC","Venous_EC","Tip_EC","Immature_EC","Other_EC","LEC"))
{
  cell=i
  plot_ls[[i]] <- seu_endo@meta.data %>% 
    dplyr::filter(Type == 'PT') %>%
    dplyr::filter(Cluster == 'Endo') %>% 
    droplevels() %>% 
    dplyr::mutate(Sample_ID = paste(Cancer.type, Sample_ID, sep = '-')) %>%
    dplyr::count(Cancer.type, Sample_ID, MidCluster, name = "n_cells") %>%
    tidyr::complete(Cancer.type, Sample_ID, MidCluster, fill = list(n_cells = 0)) %>% 
    dplyr::group_by(Cancer.type, Sample_ID) %>%
    dplyr::filter(sum(n_cells) > 0) %>%
    dplyr::mutate(prop_cells = n_cells / sum(n_cells)) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(MidCluster == i) %T>% 
    {anova_result <<- aov(prop_cells ~ Cancer.type, data = .)
    f_value <<- summary(anova_result)[[1]]$`F value`[1]
    p_value <<- summary(anova_result)[[1]]$`Pr(>F)`[1]
    } %>% 
    {ggplot(., aes(x = Cancer.type, y = prop_cells, color = Cancer.type))+
        geom_boxplot(width = 0.8, position = position_dodge(0), fill = 'transparent', outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(size = 0.8, alpha = 0.9, width = 0.4, method = "quasirandom") +
        annotate("text", x = 1, y = max(.$prop_cells) * 0.95, parse = TRUE, hjust = 0, 
                 label = sprintf("italic(p)==%s*','~italic(F)==%.2f", signif(p_value, 3), f_value)) +
        scale_y_continuous(labels = scales::percent_format(scale = 100)) +
        scale_color_manual(values = color_list, name = "Cancer Type") +
        labs(x = NULL, y = '% among TECs',title = cell) +
        theme_classic(base_size = 14) +
        theme(axis.text = element_text(color = 'black'),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 13),
              legend.position = 'none',
              strip.background = element_blank())}
  
  
  print(plot_ls[[i]])
}

total_tab<-rbind(plot_ls[[1]]$data,plot_ls[[2]]$data)
total_tab<-rbind(total_tab,plot_ls[[3]]$data)
total_tab<-rbind(total_tab,plot_ls[[4]]$data)
total_tab<-rbind(total_tab,plot_ls[[5]]$data)
total_tab<-rbind(total_tab,plot_ls[[6]]$data)


FigF_K <- cowplot::plot_grid(plot_ls, nrow = 2)

write.table(total_tab,"source_data_NC/Figure_S3D-I.txt",sep = "\t",row.names = F,col.names = T)
