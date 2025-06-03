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
library(patchwork)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")



combined<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")



Feature_P1<-FeaturePlot(combined,features = c("ANGPT2","PGF","MCAM"),pt.size = 0.1,raster = F,
                        cols = viridis::viridis(3,option = "D"),ncol=3) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        #axis.text.x = element_text(size = 10),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
Feature_P1

Idents(combined)<-combined$TopCluster
combined<-subset(combined,idents="PVL")

data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("PGF")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Type),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]))


box_p2<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 2.5, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of PGF")+
  xlab("Expression level")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,"PuBu"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlim(c(-2 ,6))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p2



gene_name<- c("ANGPT2")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Type),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]))


box_p3<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 2.5, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of ANGPT2")+
  xlab("Expression level")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,"PuBu"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlim(c(-2 ,6))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p3



##Bevtreatment
data<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/GSE140082_geo.normdata.csv",sep=",",header = T)

colnames(data)[1]<-"ID"

ID_map<-data.table::fread("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/GPL14951-11332.txt")
ID_map<-ID_map[,c("ID","ILMN_Gene")]
data<-merge(ID_map,data,by="ID")

data<-data[which(!(duplicated(data$ILMN_Gene))),]
data<-data.frame(data)
rownames(data)<-data$ILMN_Gene

data<-data[,which(!(grepl("pval",colnames(data))))]
data<-data[,c(-1,-2)]



Gene_name="MCAM+ imPC"
select_gene=Gene_name


##survival_analysis
surv_table<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE140082/Survival.txt",header = T,stringsAsFactors = F,sep="\t")
surv_table<-surv_table[which(surv_table$Treatment=="bevacizumab"),]
surv_table<-na.omit(surv_table)


mymatrix<-as.matrix(data)

Tip_Markers<-c("RGS5","MYL9","EGFL6","ACTA2","MCAM","PGF","ANGPT2","THY1","PDGFRB","COL4A1","ARHGDIB","NOTCH3")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))


##Cut off by OS
merged_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")

library(survminer)
library(survival)
res.cut <- surv_cutpoint(merged_matrix, #数据集
                         time = "OS_Time", #生存状态
                         event = "OS", #生存时间
                         variables = c("gene") #需要计算的数据列名
)
merged_matrix$gene_level<-"Low"
merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))

surv_fit<-survfit(Surv(OS_Time , OS) ~ gene_level,data= merged_matrix)

Bev_gg_surv1<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#224FA2","#E72C19"),
                         title = paste("OS by ",Gene_name,", Bev treatments\n Ovrian cancer (GSE140082)",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "d_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 1470), 
                         break.time.by=365.25,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "OS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("Low","High"),
                         #risk table
                         risk.table = T,# Add risk table
                         risk.table.font=3,
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv1$plot<-Bev_gg_surv1$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_text(size=10),
                                                           legend.position = "none")

Bev_gg_surv1$table<-Bev_gg_surv1$table+theme_classic(base_size = 10)+theme(plot.title = element_text(hjust = 0.5,size=10))
print (Bev_gg_surv1)





surv_fit<-survfit(Surv(PFS_Time , PFS) ~ gene_level,data= merged_matrix)

Bev_gg_surv2<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#224FA2","#E72C19"),
                         title = paste("PFS by ",Gene_name,", Bev treatments\n Ovrian cancer (GSE140082)",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "d_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 1470), 
                         break.time.by=365.25,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "PFS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("Low","High"),
                         #risk table
                         risk.table = T,# Add risk table
                         risk.table.font=3,
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv2$plot<-Bev_gg_surv2$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_text(size=10),
                                                           legend.position = "none")

Bev_gg_surv2$table<-Bev_gg_surv2$table+theme_classic(base_size = 10)+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10))

print (Bev_gg_surv2)






##Bevtreatment


data<-read.table("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE69795/GSE69795_series_matrix.txt",sep="\t",header = T)


colnames(data)[1]<-"ID"
data$ID<-as.character(data$ID)



ID_map<-data.table::fread("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE69795/GPL14951-11332.txt")
library(limma)
ID_map<-ID_map[,c("ID","ILMN_Gene")]
#ID_map$ILMN_Gene<-strsplit2(ID_map$ILMN_Gene," // ")[,2]

ID_map<-data.frame(ID_map)

data<-merge(ID_map,data,by="ID")
data<-data.frame(data)


data<-aggregate(data[,3:ncol(data)],list(data$ILMN_Gene), mean)
rownames(data)<-data$Group.1
data<-data[,-1]

head(data)


Gene_name="MCAM+ imPC"
select_gene=Gene_name


##survival_analysis
surv_table<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE69795/sample_info1.txt",header = T,stringsAsFactors = F,sep="\t")
#surv_table<-surv_table[which(surv_table$Treatment=="bevacizumab"),]
surv_table<-na.omit(surv_table)


mymatrix<-as.matrix(data)

Tip_Markers<-c("RGS5","MYL9","EGFL6","ACTA2","MCAM","PGF","ANGPT2","THY1","PDGFRB","COL4A1","ARHGDIB","NOTCH3")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))


##Cut off by OS
merged_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")

library(survminer)
library(survival)
res.cut <- surv_cutpoint(merged_matrix, #数据集
                         time = "OS", #生存状态
                         event = "OS_status", #生存时间
                         variables = c("gene") #需要计算的数据列名
)
merged_matrix$gene_level<-"Low"
merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))

surv_fit<-survfit(Surv(OS , OS_status) ~ gene_level,data= merged_matrix)

Bev_gg_surv3<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#224FA2","#E72C19"),
                         title = paste("OS by ",Gene_name,", Bev treatments\n Urothelial Cancer (GSE69795)",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "m_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 60), 
                         break.time.by=12,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "OS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("Low","High"),
                         #risk table
                         risk.table = T,# Add risk table
                         risk.table.font=3,
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv3$plot<-Bev_gg_surv3$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_text(size=10),
                                                           legend.position = "none")

Bev_gg_surv3$table<-Bev_gg_surv3$table+theme_classic(base_size = 10)+theme(plot.title = element_text(hjust = 0.5,size=10))
print (Bev_gg_surv3)













##Bevtreatment

data<-read.table("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE72951_series_matrix.txt/GSE72951_series_matrix.txt",sep="\t",header = T)


colnames(data)[1]<-"ID"
data$ID<-as.character(data$ID)



ID_map<-data.table::fread("~/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE72951_series_matrix.txt/GPL14951-11332.txt")
library(limma)
ID_map<-ID_map[,c("ID","ILMN_Gene")]
#ID_map$ILMN_Gene<-strsplit2(ID_map$ILMN_Gene," // ")[,2]

ID_map<-data.frame(ID_map)

data<-merge(ID_map,data,by="ID")
data<-data.frame(data)


data<-aggregate(data[,3:ncol(data)],list(data$ILMN_Gene), mean)
rownames(data)<-data$Group.1
data<-data[,-1]

head(data)


Gene_name="MCAM+ imPC"
select_gene=Gene_name


##survival_analysis
surv_table<-read.table("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/8.RNA_Seq/GSE72951_series_matrix.txt/sample_info.txt",header = T,stringsAsFactors = F,sep="\t")
#surv_table<-surv_table[which(surv_table$Treatment=="bevacizumab"),]
surv_table<-na.omit(surv_table)
colnames(surv_table)[1]<-"Tumor_ID"
colnames(surv_table)[6:7]<-c("OS","OS_status")

mymatrix<-as.matrix(data)

Tip_Markers<-c("RGS5","MYL9","EGFL6","ACTA2","MCAM","PGF","ANGPT2","THY1","PDGFRB","COL4A1","ARHGDIB","NOTCH3")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))


##Cut off by OS
merged_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")

library(survminer)
library(survival)
res.cut <- surv_cutpoint(merged_matrix, #数据集
                         time = "OS", #生存状态
                         event = "OS_status", #生存时间
                         variables = c("gene") #需要计算的数据列名
)
merged_matrix$gene_level<-"Low"
merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))

surv_fit<-survfit(Surv(OS , OS_status) ~ gene_level,data= merged_matrix)

Bev_gg_surv4<-ggsurvplot(surv_fit,
                         conf.int = F,
                         #fun = "cumhaz",
                         linetype =  1, # Change line type by groups
                         size=0.5,
                         censor = F,
                         #surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(),# Change ggplot2 theme
                         
                         palette = c("#224FA2","#E72C19"),
                         title = paste("OS by ",Gene_name,", Bev treatments\n Glioblastoma (GSE72951)",sep=""),
                         #font.family = "Arial",
                         #axis
                         xscale = "m_y",
                         pval = T,
                         surv.scale = "percent",
                         xlim = c(0, 36), 
                         break.time.by=12,
                         xlab = "Time from diagnosis (years)",
                         #ylim = c(0, 0.05),
                         break.y.by=NULL,
                         ylab = "OS rate (%)",
                         #legend
                         legend = c(0.8,0.8),
                         legend.title = Gene_name,
                         legend.labs = c("Low","High"),
                         #risk table
                         risk.table = T,# Add risk table
                         risk.table.font=3,
                         
                         #字体
                         font.tickslab = c(10, "black"),
                         font.x = c(10, "black"),
                         font.y = c(10, "black"),
                         font.main = c(10, "black"),
                         font.legend = c(10, "black"),
)
Bev_gg_surv4$plot<-Bev_gg_surv4$plot+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                                           legend.title = element_text(size=10),
                                                           legend.position = c(0.7,0.7))

Bev_gg_surv4$table<-Bev_gg_surv4$table+theme_classic(base_size = 10)+theme(plot.title = element_text(hjust = 0.5,size=10))
print (Bev_gg_surv4)





library(ggpubr)
library(survival)
library(survminer)
data<-read.table("~/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE115944/GSE115944_series_matrix.txt",sep="\t",header = T)

gene_name<-"Mcam"

colnames(data)[1]<-"ID"

ID_map<-data.table::fread("~/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE115944/GPL20258-1930.txt")
library(limma)
ID_map$ILMN_Gene<-strsplit2(ID_map$gene_assignment," // ")[,2]
ID_map<-ID_map[,c("ID","ILMN_Gene")]
data<-merge(ID_map,data,by="ID")

data<-data[which(!(duplicated(data$ILMN_Gene))),]
data<-data.frame(data)
rownames(data)<-data$ILMN_Gene

#data<-data[,which(!(grepl("pval",colnames(data))))]
data<-data[,c(-1,-2)]



mymatrix<-as.matrix(data)

Tip_Markers<-c("Rgs5","Myl9","Higd1b","Aata2","Mcam","Pgf","Angpt2","Thy1")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=40)


mt<-es.dif[1,]



#Kdr
gene_name<-"MCAM+ imPC"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))

select_fpkm_matrix$Type<-"Control"
select_fpkm_matrix$Type[which(grepl("VEGF",select_fpkm_matrix$Tumor_ID))]<-"Post-AAT"

p1<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "wilcox.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="ssGSVA score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())

p1





#Flt1
gene_name<-"Mcam"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(data[gene_name,]))

select_fpkm_matrix$Type<-"Control"
select_fpkm_matrix$Type[which(grepl("VEGF",select_fpkm_matrix$Tumor_ID))]<-"Post-AAT"

p2<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  ggtitle(gene_name)+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "t.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="Expression level")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())




#Flt1
gene_name<-"Angpt2"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(data[gene_name,]))

select_fpkm_matrix$Type<-"Control"
select_fpkm_matrix$Type[which(grepl("VEGF",select_fpkm_matrix$Tumor_ID))]<-"Post-AAT"

p3<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  ggtitle(gene_name)+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "t.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="Expression level")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())




total_p1_left<-ggarrange(p1,p2,p3,ncol = 3,common.legend = T,legend = "right")












library(ggpubr)
library(survival)
library(survminer)
data<-read.table("~/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE31715/GSE31715_series_matrix.txt",sep="\t",header = T)


colnames(data)[1]<-"ID"
data$ID<-as.character(data$ID)

ID_map<-data.table::fread("~/data/single_cell/18.pan_Endo/Non_immune/8.RNA_Seq/GSE31715/GPL6947-13512.txt")
library(limma)
ID_map<-ID_map[,c("ID","ILMN_Gene")]
data<-merge(ID_map,data,by="ID")
data<-data.frame(data)

data<-data[which(!(duplicated(data$ILMN_Gene))),]
data<-data.frame(data)
rownames(data)<-data$ILMN_Gene

#data<-data[,which(!(grepl("pval",colnames(data))))]
data<-data[,c(-1,-2)]
#data<-data[,which(grepl("Vehicle|Bevacizumab",colnames(data)))]





mymatrix<-as.matrix(data)

Tip_Markers<-c("RGS5","MYL9","HIGD1B","ACTA2","MCAM","PGF","ANGPT2","THY1")



mysymbol<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)


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
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=40)


mt<-es.dif[1,]



#Kdr
gene_name<-"MCAM+ imPC"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(mt))

select_fpkm_matrix$Type<-"Poor"
select_fpkm_matrix$Type[which(grepl("Good",select_fpkm_matrix$Tumor_ID))]<-"Good"

p1<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "wilcox.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="ssGSVA score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())

p1





#Flt1
gene_name<-"MCAM"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(data[gene_name,]))

select_fpkm_matrix$Type<-"Poor"
select_fpkm_matrix$Type[which(grepl("Good",select_fpkm_matrix$Tumor_ID))]<-"Good"

p2<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "t.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="Expression level")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())



#Flt1
gene_name<-"ANGPT2"

select_fpkm_matrix<-data.frame(Tumor_ID=colnames(data),
                               gene=as.numeric(data[gene_name,]))

select_fpkm_matrix$Type<-"Poor"
select_fpkm_matrix$Type[which(grepl("Good",select_fpkm_matrix$Tumor_ID))]<-"Good"

p3<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_jitter()+
  geom_boxplot()+
  scale_color_manual(values = c("#E64B35","#4DBBD5"))+
  ggtitle(gene_name)+
  stat_compare_means(method = "t.test",label.y.npc = 0.9,label = "p.format")+
  theme_classic()+
  labs(x="Response",y="Expression level")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        #axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.background = element_blank())





total_p1_right<-ggarrange(p1,p2,p3,ncol = 3,common.legend = T,legend = "right")









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



combined<-merge(sub_sce1,sub_sce2)

combined$Cluster<-factor(combined$Cluster,levels=c( "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                    "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"))

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p


library(ggcorrplot2)
cor_p1<-ggcorrplot(cor_value, method = "ellipse",type = "lower",p.mat = cor_test_mat, 
                   insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001),col = rev(c("#CA0020" ,"#F4A582", "#F7F7F7", "#92C5DE" ,"#0571B0")))





cell_pair_list<-list(c("PVL","Endo"),
                     c("Endo","CD4_conv"),
                     c("Endo","CD4_Treg"),
                     c("Endo","CD8"))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=1,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste("Correlation of",Cluster1,"and",Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(plot.title = element_text(size = 10,hjust = 0.5),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.line = element_blank(), 
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_blank(), 
          legend.title = element_text(size = 10),
          legend.position = "none",
          legend.text = element_text(size = 9))+ 
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p2<-p1+p2+p3+p4+plot_layout(ncol=2)
cor_p2<-p1


cellchat<-readRDS("6.cell_chat/3.SubCluster/cellchat.RDS")

groupSize <- as.numeric(table(cellchat@idents))

levels(cellchat@idents)
vertex.receiver = seq(1,4)


pathways.show <- "VEGF"
select_list<-"VEGF|ANGPT2"
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])
#pairLR.use$interaction_name<-unique(sort(pairLR.use$interaction_name))
pairLR.use<-pairLR.use[!(grepl("CD99|THBS1|HLA",pairLR.use$interaction_name)),]

pairLR.use<-data.frame(interaction_name=pairLR.use)
netVisual_b1<-netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = "imPC-MCAM",
                               targets.use = c("Artery_ECs"    ,        "Capillary_ECs"  ,     "Immature_ECs"   ,  "LECs"    ,       
                                               "Other_ECs"      ,  "Tip_ECs"   ,      
                                               "Venous_ECs"),
                               remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9,angle = 45,vjust = 1,hjust = 1),
        axis.text.y =element_text(size = 9), 
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8))+
  coord_flip()
netVisual_b1

netVisual_b2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use, targets.use  = "Tip_ECs",
                               sources.use = c("mPC-ACTG2", "mPC-MYH11" ,"imPC-CD36", "imPC-MCAM", "PC-prolif" ),
                               remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "none",
        legend.text = element_text(size = 8))+
  coord_flip()
netVisual_b2


box_p<-box_p2+box_p3+plot_layout(ncol=2,widths = c(1,1))
total_p3<-Bev_gg_surv1$plot+Bev_gg_surv2$plot+Bev_gg_surv3$plot+Bev_gg_surv4$plot+
  Bev_gg_surv1$table+Bev_gg_surv2$table+
  Bev_gg_surv3$table+Bev_gg_surv4$table+
  plot_layout(ncol=4,nrow = 2,heights = c(3,1))
total_p1<-ggarrange(box_p,netVisual_b2+netVisual_b1,widths = c(1,1.5))


total_p<-ggarrange(total_p1,total_p3,nrow = 2,heights = c(2.5,4))
pdf("14.Figure/Figure_5.pdf",width = 12,height = 6.5)
print(total_p)
dev.off()