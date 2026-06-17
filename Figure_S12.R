library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(PCAtools)
library(ggpubr)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(clusterProfiler)
library(msigdbr)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")




library(monocle3)
cds<-readRDS("5.monocle2/7.PC_monocle3/PC_monocle3.rds")

colSet<-c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
          brewer.pal(8, "Set1"),
          brewer.pal(8, "Dark2")[1],
          "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(colSet)<-c(   "mPC"  ,   "imPC"  ,   "prolif"   )

colors <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))

cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, show_trajectory_graph = T,
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  scale_color_gradientn(colors =colors)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5,size = 12))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        legend.title  = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by Monocle3')
cds_p1


if(F){
  sub_sce<-readRDS("3.Cluster/13.Annotation/3.PC_annotation_new.rds")
  sub_sce<-subset(sub_sce,subset=SubCluster %in% c("PC_C3_imPC-CD36" ,"PC_C2_mPC-MYH11", "PC_C4_imPC-MCAM","PC_C5-prolif"))
  sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
  
  
  
  sub_sce$Cluster<-"mPC"
  sub_sce$Cluster[grepl("imPC",sub_sce$SubCluster)]<-"imPC"
  sub_sce$Cluster[grepl("prolif",sub_sce$SubCluster)]<-"prolif"
  sub_sce$Cluster<-factor(sub_sce$Cluster,levels =  c("mPC","imPC","prolif"))
  
  saveRDS(sub_sce,"5.monocle2/7.PC_monocle3/sub_sce.RDS")
}

sub_sce<-readRDS("5.monocle2/7.PC_monocle3/sub_sce.RDS")
Cluster<-unique(sort(sub_sce$Cluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set2"),
               RColorBrewer::brewer.pal(8,"Set1"))
names(g.colSet1)<-Cluster
#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(Cluster) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

umap_p1<-ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=Cluster))+
  geom_point(aes(color=Cluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                     labels= c( "mPC","imPC","prolif"),
                     values = g.colSet1)+
  ggtitle("Clustering of pericytes")+
  geom_text(aes(label = Cluster), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 9))+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=1))


umap_p1


sub_sce<-readRDS("9.CytoTRACE2/CytoTRACE2_obj.RDS")

color1 <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))

CytoTRACE2_p1<-ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2))+
  geom_point(aes(color=CytoTRACE2_Relative),size=1) +
  scale_color_gradientn(colors =rev(color1))+
  
  ggtitle("Clustering of pericytes")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 9))+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=1))


CytoTRACE2_p1




test<-readRDS("5.monocle2/6.PC/monocle2.RDS")
test$Cluster[which(test$SubCluster=="PC_C5-prolif")]<-"prolif"
detach("package:monocle3", unload = TRUE)
library(monocle)
monocle2_time<-plot_cell_trajectory(test,color_by = "Pseudotime",
                                    cell_size =0.4)+
  theme(legend.position = "right")+
  scale_color_gradientn(colors = colors)+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title =  element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by Monocle2')
monocle2_time


monocle2_cluster<-plot_cell_trajectory(test,color_by = "Cluster",
                                       cell_size =0.4)+
  theme(legend.position = "bottom")+
  scale_color_manual(breaks = c("mPC","imPC","prolif"),
                     labels = c("mPC","imPC","prolif"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by cluster')+
  guides(colour = guide_legend(ncol=3,override.aes = list(size=3)))
monocle2_cluster


library(slingshot)
sim<-readRDS("5.monocle2/8.Slingshot/Slingshot.RDS")
plot_data<-data.frame(reducedDims(sim)$PCA)
plot_data$Pseudotime<-sim$slingPseudotime_1
line_data<-SlingshotDataSet(sim)@curves$Lineage1$s
color1 <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

library(ggpubr)
Slingshot_time<-ggplot(data=plot_data,aes(PC_1,PC_2))+
  theme_classic2()+
  geom_point(aes(PC_1,PC_2,fill = Pseudotime),shape=21,stroke=0,color="#FFFFFF00")+
  scale_fill_gradientn(colours =rev(color1))+
  geom_smooth(data=line_data,aes(PC_1,PC_2),color="red3")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by Slingshot')+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))
Slingshot_time



plot_data$Type<-sim$Type
plot_data$Type[which(plot_data$Type=="PN")]<-"Normal"
plot_data$Type[which(plot_data$Type!="Normal")]<-"Tumor"
line_data<-SlingshotDataSet(sim)@curves$Lineage1$s
color1 <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

library(ggpubr)
Slingshot_time1<-ggplot(data=plot_data,aes(PC_1,PC_2))+
  theme_classic2()+
  geom_point(aes(PC_1,PC_2,fill = Pseudotime),shape=21,stroke=0,color="#FFFFFF00",size=0.8)+
  scale_fill_gradientn(colours =rev(color1))+
  geom_smooth(data=line_data,aes(PC_1,PC_2),color="red3")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        legend.position = "none",
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by tissue type')+
  facet_wrap(~Type,ncol=1)+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))
Slingshot_time1



combined <- readRDS("3.Cluster/13.Annotation/3.PC_annotation_new.rds")

sub_sce<-subset(combined,subset=SubCluster %in% c("PC_C3_imPC-CD36" ,"PC_C2_mPC-MYH11", "PC_C4_imPC-MCAM","PC_C5-prolif"))

sub_sce$Cluster<-"mPC"
sub_sce$Cluster[grepl("imPC",sub_sce$SubCluster)]<-"imPC"
sub_sce$Cluster[grepl("prolif",sub_sce$SubCluster)]<-"imPC"
sub_sce<-subset(sub_sce,subset=Cluster %in% c("mPC","imPC"))

combined<-sub_sce
if (T){
  sample_ID<-runif(20000,1,ncol(combined))
  combined$merge_var<-"Drop"
  sample_ID<-read.table("5.monocle2/6.PC/Cell_ID.txt",header = F,row.names = 1)
  sample_ID<-as.numeric(sample_ID[,1])
  combined$merge_var[sample_ID]<-"Keep"
  Idents(combined)<-combined$merge_var
  combined<-subset(combined,idents="Keep")
}


RNA_matrix<-combined@assays$RNA@counts


mymatrix<-as.matrix(RNA_matrix)



m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)
m_t2g$gene_name<-m_t2g$gene_symbol
m_t2g$gs_name<-gsub("HALLMARK_","",m_t2g$gs_name)

type <- unique(m_t2g$gs_name)
type
gs <- list()
for (i in type){
  tmp <- m_t2g$gene_symbol[which(m_t2g$gs_name == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs


library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=10)

counts <- es.dif[,rownames(plot_data)]

pseudotime <- sim$slingPseudotime_1      # 每个细胞的伪时间
cell_Weight <- slingCurveWeights(sim)           # 细胞在各轨迹的权重（分支时用）

#top_counts <- counts[top_genes, ]
library(tradeSeq)
sce <- fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_Weight,
  nknots = 6,  # 【关键】knots数量，基于evaluateK结果调整（如k=6）
  verbose = FALSE
)


# 1. 分析基因与伪时间的整体相关性（associationTest）
Relation <- associationTest(sce)  # p值越小，基因与伪时间关联越强
# 2. 分析轨迹起止点的差异基因（startVsEndTest）
sET <- startVsEndTest(sce)  # 比较“起始细胞”与“终止细胞”的基因表达
# 关键指标：waldStat（值越大，起止点差异越显著）、logFClineage1（起始→终止的log2倍变化）
# 3. 筛选top差异基因（以waldStat排序）
order_sET <- order(sET$waldStat, decreasing = TRUE)
topGeneStart <- names(sce)[order_sET[1]]  # 差异最显著的第一个基因（如示例中的"FTL"）

plotSmoothers(sce, counts, gene = topGeneStart)
plotSmoothers(sce, counts, gene = "MYOGENESIS",)
plotSmoothers(sce, counts, gene = "MYOGENESIS")
plotSmoothers(sce, counts, gene = "NOTCH_SIGNALING")



# 准备数据：整合细胞类型、伪时间、基因表达
coldata <- data.frame(
  celltype = sim@colData$Cluster,
  row.names = colnames(sim)
)


# 匹配sce的细胞顺序
filter_coldata <- data.frame(coldata[colnames(sce), ])
# 添加伪时间信息（第一条轨迹）
filter_coldata$Pseudotime <- sce$crv$pseudotime
# 提取TOP5差异基因的表达量（log2转换，避免数值过大）
top5 <- names(sce)[order_sET[1:5]]
top5_exp <-sce@assays@data$counts %>% t()  # 转置为行=细胞，列=基因
# 合并最终绘图数据
plt_data <- cbind(filter_coldata, top5_exp)

gene<-"ANGIOGENESIS" 
trade_p1 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("Score")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p1


gene<-"NOTCH_SIGNALING" 
trade_p2<- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("Score")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p2


gene<-"MYOGENESIS" 
trade_p3<- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3",method="glm") +     # 加平滑线（无置信区间）
  ylab("Score")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p3



sub_sce<-readRDS("5.monocle2/7.PC_monocle3/sub_sce.RDS")
counts <- sub_sce@assays$RNA@counts[,rownames(plot_data)]


# 准备数据：整合细胞类型、伪时间、基因表达
coldata <- data.frame(
  celltype = sim@colData$Cluster,
  row.names = colnames(sim)
)


# 匹配sce的细胞顺序
filter_coldata <- data.frame(coldata[colnames(sce), ])
# 添加伪时间信息（第一条轨迹）
filter_coldata$Pseudotime <- sce$crv$pseudotime
# 提取TOP5差异基因的表达量（log2转换，避免数值过大）
top5 <- names(sce)[order_sET[1:5]]
top5_exp <- log2(counts[, colnames(sce)] + 0.1) %>% t()  # 转置为行=细胞，列=基因
# 合并最终绘图数据
plt_data <- cbind(filter_coldata, top5_exp)

gene<-"NOTCH3" 
trade_p4 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("log2(Expression+0.1)")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p4


gene<-"HES4" 
trade_p5 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("log2(Expression+0.1)")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p5

gene<-"HEYL" 
trade_p6 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("log2(Expression+0.1)")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p6



gene<-"RBPJ" 
trade_p7 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("log2(Expression+0.1)")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p7




gene<-"HEY1" 
trade_p8 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3") +     # 加平滑线（无置信区间）
  ylab("log2(Expression+0.1)")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p8


library(CellChat)
cellchat<-readRDS("6.cell_chat/4.Cancer.type/PT.RDS")
select_list<-"NOTCH|LGA|ANGPT|PD"
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])

groupSize <- as.numeric(table(cellchat@idents))
levels(cellchat@idents) 
vertex.receiver = seq(1,4) # a numeric vector
cellchatp1<-netVisual_bubble(cellchat, sources.use = c("Epi_BRCA","Epi_CC","Epi_COAD","Epi_ESCC","Epi_GC",
                                                       "Epi_HCC","Epi_HNSCC","Epi_ICC","Epi_LUAD","Epi_OC",
                                                       "Epi_PCA","Epi_PDAC","Epi_PTC"), targets.use = "PC", 
                             pairLR.use = pairLR.use,
                             remove.isolate = FALSE)+
  ggtitle("Cell-cell interactions of tumor cells to pericytes")+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8))



sub_sce<-readRDS("../Non_immune/3.Cluster/5.SubAnnotation/1.Epi/sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$SubCluster

DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

sub_sce<-subset(sub_sce,subset=Type!="PN")
Dot2<-DotPlot(sub_sce,features =c("DLL1","DLL3","DLL4","JAG1","JAG2"),
              group.by = "Cancer.type")+
  theme_bw()+
  coord_flip()+
  ggtitle("Expression of Notch ligands in tumor cells")+
  RotatedAxis()+
  scale_color_gradientn(colors  = c(RColorBrewer::brewer.pal(4,"Reds")))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5))
Dot2








RNA_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.htseq_fpkm-uq.tsv",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix)
gene_probe<-data.table::fread("~/data/TCGA/pancancer_expression/gencode.v22.annotation.gene.probeMap",header = T,stringsAsFactors = F)
phenotype_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.basic_phenotype.tsv",header = T,stringsAsFactors = F)
colnames(gene_probe)[1]<-"Ensembl_ID"
colnames(RNA_matrix)[1]<-"Ensembl_ID"
gene_probe<-gene_probe[,c(1,2)]
RNA_matrix<-merge(gene_probe,RNA_matrix,by="Ensembl_ID")
RNA_matrix<-data.frame(RNA_matrix)
RNA_matrix<-RNA_matrix[which(!(duplicated(RNA_matrix$gene))),]
rownames(RNA_matrix)<-RNA_matrix$gene

RNA_matrix<-RNA_matrix[,3:ncol(RNA_matrix)]
colnames(RNA_matrix)<-gsub("\\.","-",colnames(RNA_matrix))

Gene_name1="PGF"
Gene_name2="ECs"




mymatrix<-as.matrix(RNA_matrix)
imPC_Markers<-c("RGS5","FABP4","NDUFA4L2","STEAP4","HIGD1B","CD36","ARHGDIB","COX4I2","GJA4",
                "NOTCH3","COL4A1","EPS8","RGS16","CAV1","CTSC","CCDC102B","SEPT4",
                "FAM162B","KCNJ8","FABP5","MYO1B","ANGPT2","A2M",
                "COL18A1","IGFBP5","COL4A2",
                "GPX3","MCAM","PGF",
               "PDGFRB","HES4",
               "TDO2","NUDT4","HEYL","HEY1"
)  
mPC_Markers<-c("ACTA2","MYH11","ADIRF","RERGL","TAGLN","PPP1R14A","SORBS2","MYL9","PLN","NDUFA4L2","PTP4A3","BCAM","DSTN","MYLK","EPAS1","TPM2","ACTG2","TPM1","CPM")

mysymbol1<-data.frame(Gene_set="imPC",Gene_symbol=imPC_Markers)
mysymbol2<-data.frame(Gene_set="mPC",Gene_symbol=mPC_Markers)
mysymbol<-rbind(mysymbol1,mysymbol2)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)



type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs


library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=5)


Gene_name1="JAG1"
Gene_name2="imPC"
select_table<-data.frame(sample=colnames(mymatrix),
                         gene1=as.numeric(mymatrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)

project_list<-c(  "TCGA-BLCA" ,
                  "TCGA-LIHC", "TCGA-LGG", "TCGA-MESO", 
                  "TCGA-DLBC",
                  "TCGA-UVM","TCGA-KIRC","TCGA-PCPG" ,"TCGA-THYM","TCGA-TGCT")
select_table<-merge_matrix[which(merge_matrix$project_id %in% project_list),]

#绘图
TCGA_cor_p1 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="grey60",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of JAG1 and imPCs")+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=5)
print(TCGA_cor_p1)










RNA_matrix<-data.table::fread("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/Notch/GSE117083_RAW/fpkm.txt",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix,row.names = 1)
phenotype_matrix<-data.frame(sample=colnames(RNA_matrix),
                             Type=c( "Ctrl" ,"Ctrl" , "Ctrl" , "Rbpj-ko"   ,    "Rbpj-ko",  "Rbpj-ko"   ,   
                                      "Ctrl" ,"Ctrl", "Ctrl" ,"Rbpj-ko"   ,   "Rbpj-ko"   ,   "Rbpj-ko" ),
                             Time=c(rep("Day7",6),rep("Day10",6)))

mymatrix<-as.matrix(RNA_matrix)
imPC_Markers<-c("Rgs5","Fabp4","Steap4","Higd1b","Cd36","Arhgdib","Cox4i2","Gja4",
                "Notch3","Col4a1","Eps8","Rgs16","Cav1","Ctsc","Ccdc102b","Sept4",
                "Fam162b","Kcnj8","Fabp5","Myo1b","Angpt2","A2m",
                "Col18a1","Igfbp5","Col4a2",
                "Gpx3","Mcam","Pgf",
                "Pdgfrb","Hes4",
                "Tdo2","Nudt4","Heyl","Hey1"
                
)  
mPC_Markers<-c("Acta2","Myh11","Adirf","Rergl","Tagln","Ppp1r14a","Sorbs2","Myl9","Pln","Ndufa4l2","Ptp4a3","Bcam","Dstn","Mylk","Epas1","Tpm2","Actg2","Tpm1","Cpm")

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Notch_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_NOTCH_SIGNALING")]

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]


mysymbol1<-data.frame(Gene_set="imPC",Gene_symbol=imPC_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=Angio_sing)
mysymbol3<-data.frame(Gene_set="Notch_sing",Gene_symbol=Notch_sing)
mysymbol<-rbind(mysymbol1,mysymbol2)
mysymbol<-rbind(mysymbol,mysymbol3)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)



type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs


library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=5)




select_table1<-data.frame(sample=colnames(mymatrix),
                         Score=es.dif[1,],
                         Signal="imPC score")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[2,],
                          Signal="Angiogenesis score")
select_table3<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[3,],
                          Signal="Notch signalling")
select_table<-rbind(select_table1,select_table2)
select_table<-rbind(select_table,select_table3)

select_fpkm_matrix<-merge(select_table,phenotype_matrix,by="sample")

select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Type %in% c("Ctrl","Rbpj-ko")),]
select_fpkm_matrix$Time<-factor(select_fpkm_matrix$Time,levels = c("Day7","Day10"))
select_fpkm_matrix$Signal<-factor(select_fpkm_matrix$Signal,levels = c("Notch signalling","Angiogenesis score","imPC score"))
box_p5<-ggplot(data=select_fpkm_matrix,aes(x=Time,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  ggtitle("Mouse brain pericytes (Rbpj-ko vs. Ctrl)")+
  stat_compare_means(method = "t.test",label = "p.format",label.y.npc = 0.95)+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))+
  facet_wrap(~Signal,scale="free_y")
box_p5





RNA_matrix<-data.table::fread("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/Notch/GSE58368_series_matrix.txt/counts.txt",header = T,stringsAsFactors = F)

RNA_matrix<-data.frame(RNA_matrix,row.names = 1)
#RNA_matrix<-RNA_matrix[,1:9]
phenotype_matrix<-data.frame(sample=colnames(RNA_matrix),
                             Type=c( "Notch3(+)"   ,    "Notch3(+)",  "Notch3(+)"   ,   "Notch3(+)"   ,  
                                     "Notch3(-)" ,"Notch3(-)" ,  "Notch3(-)" ,"Notch3(-)","Notch3(-)",
                                     "Notch3(+)"   ,    "Notch3(+)",  "Notch3(+)" , 
                                     "Notch3(-)" ,"Notch3(-)" ,  "Notch3(-)" 
                                   ),
                             PC_source=c( "Brain PCs"   ,    "Brain PCs",  "Brain PCs"   ,   "Brain PCs"   ,  
                                          "Brain PCs" ,"Brain PCs" ,  "Brain PCs" ,"Brain PCs","Brain PCs",
                                          "Aorta PCs"   ,    "Aorta PCs",  "Aorta PCs" , 
                                          "Aorta PCs" ,"Aorta PCs" ,  "Aorta PCs" 
                             ))

mymatrix<-as.matrix(RNA_matrix)
imPC_Markers<-c("Rgs5","Fabp4","Steap4","Higd1b","Cd36","Arhgdib","Cox4i2","Gja4",
                "Notch3","Col4a1","Eps8","Rgs16","Cav1","Ctsc","Ccdc102b","Sept4",
                "Fam162b","Kcnj8","Fabp5","Myo1b","Angpt2","A2m",
                "Col18a1","Igfbp5","Col4a2",
                "Gpx3","Mcam","Pgf",
                "Pdgfrb","Hes4",
                "Tdo2","Nudt4","Heyl","Hey1"
                
)  
mPC_Markers<-c("Acta2","Myh11","Adirf","Rergl","Tagln","Ppp1r14a","Sorbs2","Myl9","Pln","Ndufa4l2","Ptp4a3","Bcam","Dstn","Mylk","Epas1","Tpm2","Actg2","Tpm1","Cpm")

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Notch_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_NOTCH_SIGNALING")]

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]


mysymbol1<-data.frame(Gene_set="imPC",Gene_symbol=imPC_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=Angio_sing)
mysymbol3<-data.frame(Gene_set="Notch_sing",Gene_symbol=Notch_sing)
mysymbol<-rbind(mysymbol1,mysymbol2)
mysymbol<-rbind(mysymbol,mysymbol3)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)



type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs


library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=5)



select_table1<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[1,],
                          Signal="imPC score")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[2,],
                          Signal="Angiogenesis score")
select_table3<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[3,],
                          Signal="Notch signalling")
select_table<-rbind(select_table1,select_table2)
select_table<-rbind(select_table,select_table3)


select_fpkm_matrix<-merge(select_table,phenotype_matrix,by="sample")

select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Type %in% c("Notch3(+)","Notch3(-)")),]
select_fpkm_matrix$Signal<-factor(select_fpkm_matrix$Signal,levels = c("Notch signalling","Angiogenesis score","imPC score"))

box_p1<-ggplot(data=select_fpkm_matrix,aes(x=PC_source,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c( "#224FA2", "#E72C19"))+
  ggtitle("Mouse pericytes (Notch3(+) vs. Notch3(-))")+
  stat_compare_means(method = "t.test",label = "p.format",label.y.npc = 0.95)+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))+
  facet_wrap(~Signal,scale="free_y")
box_p1



total_p2_left<-umap_p1+cds_p1
total_p2_right<- (monocle2_time/Slingshot_time)

total_p2<-ggarrange(total_p2_left,monocle2_cluster,total_p2_right,ncol=3,widths = c(2,1,1.2))


total_p3_left<-Slingshot_time1|(trade_p1/trade_p2)|(trade_p4/trade_p6)+plot_layout(widths = c(1,0.8,0.8))

total_p3<-ggarrange(total_p3_left,cellchatp1,ncol=2,widths = c(1.9,1.5))
total_p4<-ggarrange(Dot2,TCGA_cor_p1,widths = c(1,2))
total_p5<-box_p1+box_p5+plot_layout(ncol = 2)
total_p<-ggarrange(total_p2,total_p3,total_p4,total_p5,ncol = 1,nrow = 4,heights = c(3,3.5,4,3.5))

pdf("14.Figure/4.Figure_4_S12_imPCpseudo.pdf",width = 12,height = 14)
print(total_p)
dev.off()

