library(tidyverse)
library(patchwork)
library(Seurat)

setwd("~/data/single_cell/18.pan_Endo/PGF/")


# Figure S2 (A) ----
seu_endo <- readRDS("3.Cluster/13.Annotation/1.Endo_annotation.rds")

FigS2A <- Seurat::FeaturePlot(object = seu_endo, reduction = "tsne", ncol = 6,
                              features = c('CXCL12','GJA5','CA4','CD36','ACKR1','SELP',
                                           'PGF','INSR','APLNR','MMP2','KRT17','PROX1'),
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


data1<-data.frame(Cluster=FigS2A[[1]]$data$ident,
                  tSNE_1=FigS2A[[1]]$data$tSNE_1,
                  tSNE_2=FigS2A[[1]]$data$tSNE_2,
                  CXCL12=FigS2A[[1]]$data$CXCL12,
                  GJA5=FigS2A[[2]]$data$GJA5,
                  CA4=FigS2A[[3]]$data$CA4,
                  CD36=FigS2A[[4]]$data$CD36,
                  ACKR1=FigS2A[[5]]$data$ACKR1,
                  SELP=FigS2A[[6]]$data$SELP,
                  
                  PGF=FigS2A[[7]]$data$PGF,
                  INSR=FigS2A[[8]]$data$INSR,
                  APLNR=FigS2A[[9]]$data$APLNR,
                  MMP2=FigS2A[[10]]$data$MMP2,
                  KRT17=FigS2A[[11]]$data$KRT17,
                  PROX1=FigS2A[[12]]$data$PROX1
                  
)

data1[,2:ncol(data1)]<-round(data1[,2:ncol(data1)],2)
write.table(data1,"source_data_NC/Figure_S2A.txt",sep = "\t",row.names = F,col.names = T)



## 堆积柱状图 ----
sc_combined <- qs2::qs_read(stringr::str_glue("{work_dir}/3.Cluster/13.Annotation/sc_combined.qs2"))
meta_All <- sc_combined@meta.data %>% 
  dplyr::filter(Cluster != 'Epi') %>%
  dplyr::mutate(Sample_ID2 = paste(Cancer.type, Sample_ID, Type, sep = '-')) %>% 
  droplevels()

source('~/reference/codes/sc_pipeline/sc_cal.R')
prop_dt1 <- sc_calProportion(
  metadata = meta_All,
  celltype_col = 'Cluster',
  batch_col = 'Sample_ID2',
  group_col = c('Cancer.type'),
  keep_cols = c('Patient_ID','Type')
)

prop_dt1 <- prop_dt1 %>% 
  dplyr::select(-n_cells, -median_proportion) %>% 
  dplyr::mutate(
    Cancer.type = factor(Cancer.type, levels = sort(unique(Cancer.type))),
    Type = forcats::fct_relevel(Type, c("PN", "PT", "mLN", "mBrain")),
    Cluster = forcats::fct_relevel(Cluster, c('Endo', setdiff(all_cluster, c('Epi', 'Endo'))))
  ) %>% 
  tidyr::pivot_wider(names_from = 'Cluster', values_from = 'proportion') %>%
  dplyr::arrange(
    Cancer.type,
    Type,
    desc(Endo), desc(NF), desc(CAF), desc(PC)
  ) %>% 
  tidyr::pivot_longer(
    cols = c(NF:Mast),
    names_to = "MidCluster", values_to = "proportion"
  ) %>%
  dplyr::mutate(
    Sample_ID2 = factor(Sample_ID2, levels = unique(Sample_ID2)),
    MidCluster = forcats::fct_relevel(MidCluster, c('Endo', setdiff(all_cluster, c('Epi', 'Endo'))))
  ) %>% 
  dplyr::relocate(Type, Cancer.type)

group_color <- paletteer::paletteer_d("ggsci::category20c_d3")[1:17]
names(group_color) <- c(levels(prop_dt1$Cancer.type), levels(prop_dt1$Type))

block_df <- prop_dt1 %>% 
  dplyr::select(-c(MidCluster, proportion)) %>% 
  dplyr::rename(Tissue = Type, Cancer = Cancer.type) %>% 
  dplyr::distinct() %>% 
  tidyr::pivot_longer(!c(Sample_ID2, Patient_ID), names_to = 'AnnoType', values_to = 'Anno') %>% 
  dplyr::mutate(AnnoType = AnnoType %>% 
                  forcats::fct_relevel(c('Tissue','Cancer'))) %>% 
  dplyr::mutate(Y_position = 0.1 * as.numeric(AnnoType) + 1)

p1 <- prop_dt1 %>% 
  dplyr::arrange(Cancer.type, Type, Sample_ID2) %>% 
  ggplot() + # 改用fill而不是color
  geom_bar(aes(x = Sample_ID2, y = proportion, fill = MidCluster), stat = "identity", position = "stack") +
  scale_fill_manual(values = c(all_color)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = block_df,
            aes(x = Sample_ID2, y = Y_position, fill = Anno),
            height = 0.08) +
  geom_text(data = block_df %>% dplyr::distinct(AnnoType, Y_position),
            aes(x = Inf, y = Y_position, label = AnnoType),
            hjust = -0.1, size = 11, size.unit = 'pt') +
  scale_fill_manual(values = group_color) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     breaks = 0:4 / 4) +
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_blank(),
    plot.margin = ggplot2::margin(t = 0, r = 300, b = 0, l = 25, unit = "pt"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black'),
    legend.position = 'none'
  ) +
  labs(y = "% among non-Epi") +
  coord_cartesian(ylim = c(0, 1.2), clip = 'off')


sc_Endo <- qs2::qs_read(stringr::str_glue("{work_dir}/3.Cluster/13.Annotation/sc_Endo.qs2"))
meta_Endo1 <- sc_Endo@meta.data %>% 
  # dplyr::filter(Type %in% c('PT','PN')) %>% 
  # dplyr::filter(Cluster == 'Endo') %>% 
  dplyr::mutate(Sample_ID2 = paste(Cancer.type, Sample_ID, Type, sep = '-')) %>% 
  droplevels()

source('~/reference/codes/sc_pipeline/sc_cal.R')
prop_dt2 <- sc_calProportion(
  metadata = meta_Endo1,
  celltype_col = 'MidCluster',
  batch_col = 'Sample_ID2',
  group_col = c('Cancer.type'),
  keep_cols = c('Patient_ID','Type')
)

prop_dt2 <- prop_dt2 %>% 
  dplyr::select(-n_cells, -median_proportion) %>% 
  dplyr::mutate(
    Cancer.type = factor(Cancer.type, levels = sort(unique(Cancer.type))),
    Type = forcats::fct_relevel(Type, c("PN", "PT", "mLN", "mBrain")),
    MidCluster = forcats::fct_relevel(MidCluster, 
                                      c("Artery_EC", "Capillary_EC", "Venous_EC", "Tip_EC", "Immature_EC1", "Immature_EC2", "Other_EC", "LEC"))
  ) %>% 
  # tidyr::pivot_wider(names_from = 'MidCluster', values_from = 'proportion') %>%
  # dplyr::arrange(
  #   Cancer.type,
  #   Type,
  #   desc(Artery_EC), desc(Capillary_EC), desc(Venous_EC), desc(Tip_EC), desc(Immature_EC1), desc(Immature_EC2), desc(Other_EC), desc(LEC)
  # ) %>% 
  # tidyr::pivot_longer(
  #   cols = c(Artery_EC, Capillary_EC, Venous_EC, Tip_EC, Immature_EC1, Immature_EC2, Other_EC, LEC),
  #   names_to = "MidCluster", values_to = "proportion"
  # ) %>%
  dplyr::mutate(
    Sample_ID2 = factor(Sample_ID2, levels = levels(prop_dt1$Sample_ID2)),
    MidCluster = forcats::fct_relevel(MidCluster, 
                                      c("Artery_EC", "Capillary_EC", "Venous_EC", "Tip_EC", "Immature_EC1", "Immature_EC2", "Other_EC", "LEC"))
  ) %>% 
  dplyr::relocate(Type, Cancer.type)

# group_color <- paletteer::paletteer_d("ggsci::category20c_d3")[1:17]
# names(group_color) <- c(levels(prop_dt$Cancer.type), levels(prop_dt$Type))
# 
# block_df <- prop_dt %>% 
#   dplyr::select(-c(MidCluster, proportion)) %>% 
#   dplyr::rename(Tissue = Type, Cancer = Cancer.type) %>% 
#   dplyr::distinct() %>% 
#   tidyr::pivot_longer(!c(Sample_ID2, Patient_ID), names_to = 'AnnoType', values_to = 'Anno') %>% 
#   dplyr::mutate(AnnoType = AnnoType %>% 
#                   forcats::fct_relevel(c('Tissue','Cancer'))) %>% 
#   dplyr::mutate(Y_position = 0.05 * as.numeric(AnnoType) + 1)

p2 <- prop_dt2 %>% 
  dplyr::arrange(Cancer.type, Type, Sample_ID2) %>% 
  ggplot() + # 改用fill而不是color
  geom_bar(aes(x = Sample_ID2, y = proportion, fill = MidCluster), stat = "identity", position = "stack") +
  scale_fill_manual(values = c(endo_color)) +
  # ggnewscale::new_scale_fill() +
  # geom_tile(data = block_df,
  #           aes(x = Sample_ID2, y = Y_position, fill = Anno),
  #           height = 0.04) +
  # geom_text(data = block_df %>% dplyr::distinct(AnnoType, Y_position),
  #           aes(x = Inf, y = Y_position, label = AnnoType), 
  #           hjust = -0.1, size = 11, size.unit = 'pt') +
  # scale_fill_manual(values = group_color) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     breaks = 0:4 / 4) +
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_blank(),
    plot.margin = ggplot2::margin(t = 0, r = 300, b = 0, l = 25, unit = "pt"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black'),
    legend.position = 'none'
  ) +
  labs(y = "% among TECs") +
  coord_cartesian(ylim = c(0, 1), clip = 'off')


library(ComplexHeatmap)

group_color <- paletteer::paletteer_d("ggsci::category20c_d3")[1:17]
names(group_color) <- c(levels(prop_dt1$Cancer.type), levels(prop_dt1$Type))

lgd <- list(
  Legend(title = "Cell Type\n(All Cells)", labels = names(all_color), legend_gp = gpar(fill = unname(all_color))),
  Legend(title = "Cell Type\n(Endothelial)", labels = names(endo_color), legend_gp = gpar(fill = unname(endo_color))),
  Legend(title = "Tissue", labels = levels(prop_dt1$Type), legend_gp = gpar(fill = unname(group_color)[14:17])),
  Legend(title = "Cancer", labels = levels(prop_dt1$Cancer.type), legend_gp = gpar(fill = unname(group_color)[1:13]))
)
lgd_ls <- packLegend(list = lgd, direction = 'horizontal')


pdf(file = stringr::str_glue("{write_dir}/figS2b.Endo_stack_barplot.pdf"),
    width = 16, height = 4.5)
# print(p1/p2)
cowplot::plot_grid(
  p1, p2, ncol = 1,
  align = 'v', axis = 'lr',
  rel_heights = c(1.15, 1)
)
draw(lgd_ls, x = unit(0.99, "npc"), y = unit(0.8, "npc"), just = c("right", "top"))
dev.off()



Figure_S2B_data <- seu_merge@meta.data %>% 
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
  dplyr::ungroup() 

write.table(data1,"source_data_NC/Figure_S2B.txt",sep = "\t",row.names = F,col.names = T)



Figure_S2B_data <- seu_merge@meta.data %>% 
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
  dplyr::ungroup() 

write.table(data1,"source_data_NC/Figure_S2B_top.txt",sep = "\t",row.names = F,col.names = T)



##Stack plot

sub_sce<-readRDS("3.Cluster/13.Annotation/1.Endo_annotation.rds")
sub_sce$Type[which(sub_sce$Type!="PN")]<-"Tumor"
sub_sce$Type[which(sub_sce$Type=="PN")]<-"Normal"

dir.create("3.Cluster/13.plot/1.total/")

combined<-sub_sce

##All cell type
tmp_table<-data.frame(Cluster=combined$MidCluster,Type=combined$Type,Sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type,Sample_ID),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,Sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by=c("Type","Sample_ID"))
sum_table$per<-sum_table$num/sum_table$total_num*100


Cluster<-levels(combined$MidCluster)

g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")
names(g.colSet1)<-Cluster

sum_table1<-sum_table[sum_table$Cluster %in% c("Artery_ECs" , "Capillary_ECs",    "Venous_ECs"),]
sum_table1<-dplyr::summarize(group_by(sum_table1,Type,Sample_ID),per=sum(per))

sum_table1<-sum_table1[order(sum_table1$per,decreasing = T),]
Sample_ID<-unique(c(sum_table1$Sample_ID,sum_table$Sample_ID))
sum_table$Sample_ID<-factor(sum_table$Sample_ID,levels=Sample_ID)

sum_table$Cluster<-factor(sum_table$Cluster,levels = c("Artery_ECs" , "Capillary_ECs",    "Venous_ECs"  , "Tip_ECs"   ,  "Immature_ECs" , "LECs"   ,       "Other_ECs" ))


stack_p2<-ggplot(data=sum_table, aes(x=Sample_ID, y=per, fill=Cluster)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = g.colSet1)+
  labs(x = 'Type', y ="Percentage",title=paste("Distribution of ECs (n = 363)")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,101),breaks = seq(0,100,25),labels = c("0","25%","50%","75%","100%"))+
  theme(axis.text.y = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        axis.text.x = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10,vjust = 0.5,hjust=0.5)
  )+
  facet_grid(~Type,scales = "free_x", space = 'free_x')
print(stack_p2)


write.table(sum_table,"source_data_NC/Figure_S2B_bottm.txt",sep = "\t",row.names = F,col.names = T)




##All cell type
tmp_table<-data.frame(Cluster=combined$MidCluster,Type=combined$Type,Sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
tmp_table$merge_var<-paste(tmp_table$Sample_ID,tmp_table$Cluster)
tmp_table<-tmp_table[!(duplicated(tmp_table$merge_var)),]

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type,Sample_ID),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,Sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by=c("Type","Sample_ID"))
sum_table$per<-sum_table$num/sum_table$total_num*100


Cluster<-levels(combined$MidCluster)

g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")
names(g.colSet1)<-Cluster

sum_table1<-sum_table

sum_table1<-sum_table1[order(sum_table1$total_num,decreasing = T),]
Sample_ID<-unique(c(sum_table1$Sample_ID,sum_table$Sample_ID))
sum_table$Sample_ID<-factor(sum_table$Sample_ID,levels=Sample_ID)

sum_table$Cluster<-factor(sum_table$Cluster,levels = c("Artery_ECs" , "Capillary_ECs",    "Venous_ECs"  , "Tip_ECs"   ,  "Immature_ECs" , "LECs"   ,       "Other_ECs" ))


stack_p3<-ggplot(data=sum_table, aes(x=Sample_ID, y=num, fill=Cluster)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = g.colSet1)+
  labs(x = 'Type', y ="Number of EC lineages",title=paste("")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(axis.text.y = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        axis.text.x = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10,vjust = 0.5,hjust=0.5)
  )+
  facet_grid(~Type,scales = "free_x", space = 'free_x')
print(stack_p3)


write.table(sum_table,"source_data_NC/Figure_S2C.txt",sep = "\t",row.names = F,col.names = T)




sum_table1<-sum_table[!(duplicated(sum_table$Sample_ID)),]

sum_table2<-data.frame(table(sum_table1$total_num))
colnames(sum_table2)<-c("Name","Num")
sum_table3<-data.frame(Name=0,Num=18)

sum_table2<-rbind(sum_table3,sum_table2)
sum_table2$Total_Num<-sum(sum_table2$Num)
sum_table2$per<-round(sum_table2$Num/sum_table2$Total_Num*100,1)
sum_table2$Name<-factor(sum_table2$Name,levels = rev(sum_table2$Name))

cols<-RColorBrewer::brewer.pal(n = 8,name = "Blues")
names(cols)<-sum_table2$Name

library(scales)
pie <- ggplot(sum_table2, aes(x="", y=per, fill=factor(Name)))+
  geom_bar(width=1, stat="identity") + 
  coord_polar("y", start=0) + # 转换为极坐标（饼图核心）
  scale_fill_manual(values = cols)+
  theme(axis.line=element_blank(), plot.title=element_text(hjust=0.5))+
  labs(fill="Name", x=NULL, y=NULL, title="Total EC lineages by samples")+
  geom_text(aes(y = per/3 + c(0, cumsum(per)[-length(per)]), x = 1.3,
                label = percent(per/100)), size=3)

pie

write.table(sum_table2,"source_data_NC/Figure_S2D.txt",sep = "\t",row.names = F,col.names = T)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tSNE")<-Embeddings(object = sub_sce, reduction = "tsne")


Cluster<-levels(combined$MidCluster)

g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")
names(g.colSet1)<-Cluster



library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$tsne@cell.embeddings,col.name = colnames(sub_sce@reductions$tsne@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(MidCluster) %>%
  dplyr::summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )


umap_p2<-ggplot(sub_sce@meta.data ,aes(x=tSNE_1,y=tSNE_2,color=MidCluster))+
  geom_point(aes(color=MidCluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$MidCluster)),
                     labels= c( "Artery_ECs" , "Capillary_ECs",    "Venous_ECs"  , "Tip_ECs"   ,  "Immature_ECs" , "LECs"   ,       "Other_ECs"   ),
                     values = g.colSet1)+
  ggtitle("Clustering of Endothelial cells")+
  #geom_text(aes(label = SubCluster1), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.ticks  = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        legend.position = "right",
        strip.background = element_blank(),
        legend.text = element_text(size = 9))+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=1))+
  facet_wrap(~Cancer.type,ncol=7)


umap_p2

endo_data<-data.frame(Cancer.type=sub_sce$Cancer.type,Cluster=sub_sce$MidCluster,round(sub_sce@reductions$tsne@cell.embeddings,2))
write.table(endo_data,"source_data_NC/Figure_S2E.txt",sep = "\t",row.names = F,col.names = T)


library(patchwork)
total_p1<-stack_p2+stack_p3+plot_layout(ncol = 1,nrow = 2)
total_p<-ggarrange(total_p1,
                   pie+plot_spacer()+plot_layout(ncol = 2,nrow = 1),
                   umap_p2,
                   ncol = 1,nrow = 3,heights = c(2,1,2))
pdf("14.Figure/Figure_S2.pdf",width = 12,height = 12)
print(total_p)
dev.off()
