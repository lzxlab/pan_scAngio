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
library(SummarizedExperiment)
library(SingleCellExperiment)
library(slingshot)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")

color1 <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)


combined <- readRDS("3.Cluster/13.Annotation/3.PC_annotation_new.rds")
test<-readRDS("5.monocle2/6.PC//monocle2.RDS")
sample_ID<-colnames(test)
combined$merge_var<-"Drop"
combined$merge_var[sample_ID]<-"Keep"
Idents(combined)<-combined$merge_var
combined<-subset(combined,idents="Keep")

combined$Cluster<-"mPC"
combined$Cluster[grepl("CD36",combined$SubCluster)]<-"imPC-CD36"
combined$Cluster[grepl("MCAM",combined$SubCluster)]<-"imPC-MCAM"
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:15)


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
Hypoxia_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1","SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")

gs[["Hypoxia"]]<-Hypoxia_Markers
library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=10)

sim<-readRDS("5.monocle2/8.Slingshot/Slingshot.RDS")
plot_data<-data.frame(reducedDims(sim)$PCA)
plot_data$Pseudotime<-sim$slingPseudotime_1
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

gene<-"Hypoxia" 
trade_p1 <- ggscatter(
  data = plt_data,
  x = "Pseudotime",     # x轴=伪时间（发育顺序）
  y = gene,           # y轴=基因表达量
  color = "Pseudotime",    # 颜色=细胞类型head()
  size = 0.6         # 点大小
) +
  geom_smooth(se = F, color = "red3",span = 0.5) +     # 加平滑线（无置信区间）
  ylab("Score")+
  ggtitle(gene)+
  theme_bw() +                           # 简洁主题
  scale_color_gradientn(colours =rev(color1))+  
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text.x = element_blank(), axis.title.x =element_text(size = 10),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p1


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

gene<-"HIF1A" 
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

gene<-"PGF" 
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
  theme(axis.text = element_text(size = 10), axis.title = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p5


gene<-"ANGPT2" 
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


gene<-"HIGD1B" 
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
  theme(axis.text = element_text(size = 10), axis.title =element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+   # 自定义配色
  theme(legend.position = "none")             # 隐藏图例（避免重复）
trade_p7



##Cluster 
combined <- readRDS("5.monocle2/7.PC_monocle3/sub_sce.RDS")

combined$Cluster<-"mPC"
combined$Cluster[grepl("imPC",combined$SubCluster)]<-"imPC"
combined$Cluster[grepl("prolif",combined$SubCluster)]<-"prolif"
combined$Cluster<-factor(combined$Cluster,levels =  c("mPC","imPC","prolif"))

Gene_name1="Hypoxia"

DefaultAssay(combined)<-"RNA"

RNA_matrix<-combined@assays$RNA@counts


mymatrix<-as.matrix(RNA_matrix)
Hypoxia_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1","SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]


mysymbol1<-data.frame(Gene_set="Hypoxia",Gene_symbol=Hypoxia_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=Angio_sing)
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
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=10)




mt<-es.dif[1,]




mt<-t(es.dif)
combined<-AddMetaData(combined,metadata = mt)


Vln_exp1<-VlnPlot(combined,features = c("Hypoxia","Angiogenesis"),cols = RColorBrewer::brewer.pal(8,"Set2"),
                  pt.size=0,group.by="Cluster",ncol = 2) &
  ylim(c(0.5,1.8)) &
  stat_compare_means(method = "wilcox.test",comparisons  = list(c("mPC","imPC"))) &
  theme(axis.text = element_text(size = 10), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 10),
        #axis.text.x = element_text(size = 10,angle = 45,vjust = 0.5,hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5))
Vln_exp1








RNA_matrix<-data.table::fread("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/Hypoxia/GSE109233_non_normalized.txt/counts.txt",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix,row.names = 1)
phenotype_matrix<-data.frame(sample=colnames(RNA_matrix),
                             Type=c( rep("Hypoxia" ,8),
                                     rep("Normoxia" ,8),
                                     rep("Hypoxia" ,6),
                                     rep("Normoxia" ,6)),
                             Time=c( rep("2h" ,8),
                                     rep("2h" ,8),
                                     rep("6h" ,6),
                                     rep("6h" ,6)),
                             Glucose=c(rep("Glucose (+)",4),rep("Glucose (-)",4),
                                    rep("Glucose (+)",4),rep("Glucose (-)",4),
                                    rep("Glucose (+)",3),rep("Glucose (-)",3),
                                    rep("Glucose (+)",3),rep("Glucose (-)",3)))

mymatrix<-as.matrix(RNA_matrix)
Hypoxia_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1","SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]

Angio_sing<-c("TYMP","VCAN","CD44","FYN","VEGFA","ANGPT2","PGF","TNFAIP6","E2F3","TGFB1",
              "PDGFA","PDGFD","FGF7","SERPINA5",
              "COL3A1"  , "COL5A2"  , "CXCL6"  ,  "FGFR1"  ,  "FSTL1" ,   "ITGAV"  ,  "JAG1",
              "TIMP1"  ,  "TNFRSF21", "VAV2"   ,  "VCAN"   ,  "VEGFA" ,   "VTN", 
              "MMP9","ITGAV","SPP1","PTK2","CCND2","EZH2")
mysymbol1<-data.frame(Gene_set="Hypoxia",Gene_symbol=Hypoxia_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=c(Angio_sing))
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




select_table1<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[1,],
                          Signal="Hypoxia score")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[2,],
                          Signal="Angiogenesis score")
select_table<-rbind(select_table1,select_table2)

select_fpkm_matrix<-merge(select_table,phenotype_matrix,by="sample")

#select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Type %in% c("Ctrl","Rbpj-ko")),]
select_fpkm_matrix$Glucose<-factor(select_fpkm_matrix$Glucose,levels = c("Glucose (+)","Glucose (-)"))
select_fpkm_matrix$Time<-factor(select_fpkm_matrix$Time,levels = c("2h" ,"6h"))
select_fpkm_matrix$Signal<-factor(select_fpkm_matrix$Signal,levels = c("Hypoxia score","Angiogenesis score"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia"))

select_fpkm_matrix1<-select_fpkm_matrix[which(select_fpkm_matrix$Signal %in% "Hypoxia score"),]
box_p1<-ggplot(data=select_fpkm_matrix1,aes(x=Time,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  ggtitle("Hypoxia score in pericytes (GSE109233)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95)+
  ylab("Hypoxia score")+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))
box_p1





gene_name<- unique(c("PGF","HIF1A",
                     "BNIP3","ADM","LDHA"))

data1<-data.frame(mymatrix[gene_name,])
data1$Gene_name<-rownames(data1)
data1<-reshape2::melt(data1,varnams=c("Gene_name"),value.name="gene")
colnames(data1)<-c("gene","sample","Score")
data1$Score<-log2(data1$Score+0.01)
select_fpkm_matrix<-merge(data1,phenotype_matrix,by="sample")


#select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Type %in% c("Ctrl","Rbpj-ko")),]
select_fpkm_matrix$Glucose<-factor(select_fpkm_matrix$Glucose,levels = c("Glucose (+)","Glucose (-)"))
select_fpkm_matrix$Time<-factor(select_fpkm_matrix$Time,levels = c("2h" ,"6h"))
select_fpkm_matrix$gene<-factor(select_fpkm_matrix$gene,levels = c("PGF","HIF1A",
                                                                   "BNIP3","ADM","LDHA"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia"))

select_fpkm_matrix1<-select_fpkm_matrix[which(select_fpkm_matrix$Glucose %in% "Glucose (+)"),]
box_p3<-ggplot(data=select_fpkm_matrix1,aes(x=Time,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  ggtitle("Expression in pericytes (Hypoxia vs Normoxia, GSE109233)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95)+
  theme_classic()+
  ylab("log2(Expression+0.01)")+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))+
  facet_wrap(~gene,scale="free_y",ncol=6)
box_p3





RNA_matrix<-data.table::fread("~/data/lzx_project/19.PC_hypo_RNA/4.matrix/FPKM_symbol.txt",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix,row.names = 1)
phenotype_matrix<-data.frame(sample=colnames(RNA_matrix),
                             Type=c( rep("Hypoxia" ,4),
                                     rep("Normoxia" ,4)))

mymatrix<-as.matrix(RNA_matrix)
Hypoxia_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1",
                    "SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]

Angio_sing<-c("TYMP","VCAN","CD44","FYN","VEGFA","ANGPT2","PGF","TNFAIP6","E2F3","TGFB1",
              "PDGFA","PDGFD","FGF7","SERPINA5",
              "COL3A1"  , "COL5A2"  , "CXCL6"  ,  "FGFR1"  ,  "FSTL1" ,   "ITGAV"  ,  "JAG1",
              "TIMP1"  ,  "TNFRSF21", "VAV2"   ,  "VCAN"   ,  "VEGFA" ,   "VTN", 
              "MMP9","ITGAV","SPP1","PTK2","CCND2","EZH2")
mysymbol1<-data.frame(Gene_set="Hypoxia",Gene_symbol=Hypoxia_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=c(Angio_sing))
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




select_table1<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[1,],
                          Signal="Hypoxia score")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[2,],
                          Signal="Angiogenesis score")
select_table<-rbind(select_table1,select_table2)

select_fpkm_matrix<-merge(select_table,phenotype_matrix,by="sample")


select_fpkm_matrix$Signal<-factor(select_fpkm_matrix$Signal,levels = c("Hypoxia score","Angiogenesis score"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia"))

select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Signal %in% "Hypoxia score"),]
box_p11<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  ggtitle("Human pericytes (Hypoxia vs. Normoxia)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95)+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))
box_p11





select_table1<-data.frame(sample=colnames(mymatrix),
                          Score=log2(mymatrix["ANGPT2",]+0.01),
                          gene="ANGPT2")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=log2(mymatrix["PGF",]+0.01),
                          gene="PGF")
select_table3<-data.frame(sample=colnames(mymatrix),
                          Score=log2(mymatrix["HIF1A",]+0.01),
                          gene="HIF1A")

gene_name<- unique(c("PGF","HIF1A",
                     "BNIP3","ADM","LDHA"))

data1<-data.frame(mymatrix[gene_name,])
data1$Gene_name<-rownames(data1)
data1<-reshape2::melt(data1,varnams=c("Gene_name"),value.name="gene")
colnames(data1)<-c("gene","sample","Score")
data1$Score<-log2(data1$Score+0.01)
select_fpkm_matrix<-merge(data1,phenotype_matrix,by="sample")


select_fpkm_matrix$gene<-factor(select_fpkm_matrix$gene,levels = c("PGF","HIF1A",
                                                                   "BNIP3","ADM","LDHA"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia"))

box_p31<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  ggtitle("Expression in pericytes (Hypoxia vs. Normoxia)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95)+
  theme_classic()+
  ylab("log2(Expression+0.01)")+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))+
  facet_wrap(~gene,scale="free_y",ncol=6)
box_p31




RNA_matrix<-data.table::fread("~/data/lzx_project/25.hypo_PC/4.matrix/FPKM_symbol.txt",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix,row.names = 1)
phenotype_matrix<-data.frame(sample=colnames(RNA_matrix),
                             Type=c( rep("Hypoxia" ,3),
                                     rep("Normoxia" ,3),rep("Re-Oxygen" ,3)))

mymatrix<-as.matrix(RNA_matrix)
Hypoxia_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1",
                    "SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)

Angio_sing<-m_t2g$gene_symbol[which(m_t2g$gs_name=="HALLMARK_ANGIOGENESIS")]

Angio_sing<-c("TYMP","VCAN","CD44","FYN","VEGFA","ANGPT2","PGF","TNFAIP6","E2F3","TGFB1",
              "PDGFA","PDGFD","FGF7","SERPINA5",
              "COL3A1"  , "COL5A2"  , "CXCL6"  ,  "FGFR1"  ,  "FSTL1" ,   "ITGAV"  ,  "JAG1",
              "TIMP1"  ,  "TNFRSF21", "VAV2"   ,  "VCAN"   ,  "VEGFA" ,   "VTN", 
              "MMP9","ITGAV","SPP1","PTK2","CCND2","EZH2")
mysymbol1<-data.frame(Gene_set="Hypoxia",Gene_symbol=Hypoxia_Markers)
mysymbol2<-data.frame(Gene_set="Angiogenesis",Gene_symbol=c(Angio_sing))
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




select_table1<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[1,],
                          Signal="Hypoxia score")
select_table2<-data.frame(sample=colnames(mymatrix),
                          Score=es.dif[2,],
                          Signal="Angiogenesis score")
select_table<-rbind(select_table1,select_table2)

select_fpkm_matrix<-merge(select_table,phenotype_matrix,by="sample")


select_fpkm_matrix$Signal<-factor(select_fpkm_matrix$Signal,levels = c("Hypoxia score","Angiogenesis score"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia","Re-Oxygen"))

select_fpkm_matrix<-select_fpkm_matrix[which(select_fpkm_matrix$Signal %in% "Hypoxia score"),]
box_p41<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399","#003399"))+
  ggtitle("Human pericytes (Hypoxia vs. Normoxia)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95,
                    comparisons = list(c("Normoxia","Hypoxia"),c("Hypoxia","Re-Oxygen")))+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))
box_p41





gene_name<- unique(c("PGF","HIF1A",
                     "BNIP3","ADM","LDHA"))

data1<-data.frame(mymatrix[gene_name,])
data1$Gene_name<-rownames(data1)
data1<-reshape2::melt(data1,varnams=c("Gene_name"),value.name="gene")
colnames(data1)<-c("gene","sample","Score")
data1$Score<-log2(data1$Score+0.01)
select_fpkm_matrix<-merge(data1,phenotype_matrix,by="sample")


select_fpkm_matrix$gene<-factor(select_fpkm_matrix$gene,levels = c("PGF","HIF1A",
                                                                   "BNIP3","ADM","LDHA"))
select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels = c("Normoxia","Hypoxia","Re-Oxygen"))

box_p42<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=Score,fill=Type))+
  #geom_jitter()+
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values = c("#336699", "#993399","#003399"))+
  ggtitle("Expression in pericytes (Hypoxia vs. Normoxia)")+
  stat_compare_means(method = "t.test",label = "p.signif",label.y.npc = 0.95,
                     comparisons = list(c("Normoxia","Hypoxia"),c("Hypoxia","Re-Oxygen")))+
  theme_classic()+
  ylab("log2(Expression+0.01)")+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10))+
  facet_wrap(~gene,scale="free_y",ncol=6)
box_p42


total_p1_left<-trade_p1+(trade_p4/trade_p5)+(trade_p6/trade_p7)+plot_layout(ncol=3,widths = c(2,1,1))
total_p1<-ggarrange(total_p1_left,Vln_exp1,ncol = 2,widths = c(1.5,1))
total_p2<-box_p41+box_p42+plot_layout(ncol=2,widths = c(1,6))
total_p3<-box_p1+box_p3+plot_layout(ncol=2,widths = c(1,6))
total_p<-ggarrange(total_p1,total_p2,total_p3,ncol = 1,nrow = 3,heights = c(3,3.5,3))
pdf("14.Figure/4.Figure_4_S11_imPChypoxia.pdf",width = 12,height = 9.5)
print(total_p)
dev.off()
