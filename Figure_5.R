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


select_fpkm_matrix<-data.frame(Cluster=combined$SubCluster,
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

tab1<-data.frame(select_fpkm_matrix,gene_name="PGF")


gene_name<- c("ANGPT2")


select_fpkm_matrix<-data.frame(Cluster=combined$SubCluster,
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

tab2<-data.frame(select_fpkm_matrix,gene_name="ANGPT2")
out_t<-rbind(tab1,tab2)
out_t$Cluster<-gsub("PVL","PC",out_t$Cluster)
write.table(out_t,"source_data_NC/Figure_5A.txt",sep = "\t",row.names = F,col.names = T)



total_out<-data.frame()

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
#merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
#merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,1/3)]<-"Mid"
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,2/3)]<-"High"
merged_matrix<-merged_matrix[which(merged_matrix$gene_level!="Mid"),]
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))


write.table(merged_matrix,"source_data_NC/Figure_5G.txt",sep = "\t",row.names = F,col.names = T)



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



library(CellChat)
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

select_table<-netVisual_b1$data
write.table(select_table,"source_data_NC/Figure_5B_right.txt",sep = "\t",row.names = F,col.names = T)


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

select_table<-netVisual_b2$data
write.table(select_table,"source_data_NC/Figure_5B_left.txt",sep = "\t",row.names = F,col.names = T)





library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/HCT116_PC_antiPGF_volume.csv",sep=",",header = T)

g1<-tab$X20[which(tab$group=="HCT116")]
g2<-tab$X20[which(tab$group=="HCT116+Pericyte")]
g3<-tab$X20[which(tab$group=="HCT116+Pericyte+antiPGF")]
g4<-tab$X20[which(tab$group=="HCT116+Pericyte+antiANGPT2")]
g5<-tab$X20[which(tab$group=="HCT116+Pericyte+antiPGF+antiANGPT2")]
t.test(g1,g2)
t.test(g2,g3)
t.test(g2,g4)
t.test(g2,g5)
color_list<-c("#377EB8","#E41A1C","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")




#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))
tab1$Group<-"HCT116"
tab1$Group[tab1$group=="HCT116+Pericyte"]<-"HCT116+PCs"
tab1$Group[tab1$group=="HCT116+Pericyte+antiPGF"]<-"HCT116+PCs+antiPGF"
tab1$Group[tab1$group=="HCT116+Pericyte+antiANGPT2"]<-"HCT116+PCs+antiANGPT2"
tab1$Group[tab1$group=="HCT116+Pericyte+antiPGF+antiANGPT2"]<-"HCT116+PCs+antiPGF+antiANGPT2"

tab1$Group<-factor(tab1$Group,levels=c( "HCT116","HCT116+PCs" ,"HCT116+PCs+antiPGF","HCT116+PCs+antiANGPT2",
                                        "HCT116+PCs+antiPGF+antiANGPT2"))
library(dplyr)
result_table<-dplyr::summarize(group_by(tab1,Group,Time),
                               mean_var=mean(Value),
                               SD=sd(Value),
                               N_N=n(),  
                               se=SD/sqrt(N_N),  
                               upper_limit=mean_var+SD,  
                               lower_limit=mean_var-SD)



mice_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
  geom_point(data = result_table,aes(x=Time, y=mean_var),pch=15,size=1.5)+
  geom_line(size=0.5)+
  geom_errorbar(data = result_table,aes(ymin = lower_limit, ymax=upper_limit), #误差条表示均值±标准差
                width=0.2, #误差条末端短横线的宽度
                alpha = 0.7,
                size=0.5) +
  scale_color_manual(values = color_list,labels=c("HCT116","HCT116+PCs" ,
                                                            "HCT116+PCs\n(antiPGF)\n",
                                                            "HCT116+PCs\n(antiANGPT2)\n",
                                                            "HCT116+PCs\n(antiPGF+antiANGPT2)"))+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve\n(HCT116 with or without PCs)")+
  theme_classic2()+
  theme(legend.position=c(0.35,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.2),breaks =seq(0,1.2,0.2))
print (mice_line_p1)


write.table(result_table,"source_data_NC/Figure_5F_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(1,3,4)],"source_data_NC/Figure_5F_data.txt",sep = "\t",row.names = F,col.names = T)






##
color_list<-c( "#377EB8","#E41A1C", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/HCT116_PC_antiPGF_weight.csv",sep=",",header=T)
select_table$Type<-"HCT116"
select_table$Type[select_table$group=="HCT116+Pericyte"]<-"HCT116+PCs"
select_table$Type[select_table$group=="HCT116+Pericyte+antiPGF"]<-"HCT116+PCs+antiPGF"
select_table$Type[select_table$group=="HCT116+Pericyte+antiANGPT2"]<-"HCT116+PCs+antiANGPT2"
select_table$Type[select_table$group=="HCT116+Pericyte+antiPGF+antiANGPT2"]<-"HCT116+PCs+antiPGF+antiANGPT2"

select_table$Type<-factor(select_table$Type,levels=c( "HCT116","HCT116+PCs" ,"HCT116+PCs+antiPGF","HCT116+PCs+antiANGPT2",
                                        "HCT116+PCs+antiPGF+antiANGPT2"))
select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
PC_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,size=3,
                     comparisons = list(c("HCT116","HCT116+PCs" ),c("HCT116+PCs+antiPGF","HCT116+PCs" ),
                                        c("HCT116+PCs","HCT116+PCs+antiPGF+antiANGPT2" )
                                        ))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,size=3,
                     comparisons = list(c("HCT116+PCs+antiANGPT2","HCT116+PCs" )))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y.npc = 0.85,size=3,
                     comparisons = list(c("HCT116+PCs+antiANGPT2","HCT116+PCs+antiPGF+antiANGPT2" )))+
  
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =color_list)+
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))

print(PC_weight)

#write.table(select_table[1:2],"source_data_NC/Figure_5F_weight.txt",sep = "\t",row.names = F,col.names = T)





total_p1<-total_p1_left+total_p1_right+plot_layout(ncol=2)


total_p3<-Bev_gg_surv1$plot+Bev_gg_surv2$plot+Bev_gg_surv1$table+Bev_gg_surv2$table+plot_layout(ncol=2,nrow = 2,heights = c(3,1))



total_p4<-ggarrange(box_p,Feature_P1,ncol=1,heights = c(1.5,1))


total_p5<-ggarrange(total_p4,total_p3,ncol=2,nrow = 1,widths = c(1,1.5))

total_p<-ggarrange(total_p5,total_p1,nrow = 2,heights = c(4.5,2.5))
pdf("14.Figure/8.Figure_8.pdf",width = 12,height = 7)
print(total_p)
dev.off()


box_p<-box_p2+box_p3+plot_layout(ncol=2,widths = c(1,1))
total_p3<-Bev_gg_surv1$plot+Bev_gg_surv2$plot+Bev_gg_surv3$plot+Bev_gg_surv4$plot+
  Bev_gg_surv1$table+Bev_gg_surv2$table+
  Bev_gg_surv3$table+Bev_gg_surv4$table+
  plot_layout(ncol=4,nrow = 2,heights = c(3,1))

library(patchwork)
total_p3<-Bev_gg_surv1$plot+Bev_gg_surv2$plot+
  Bev_gg_surv1$table+Bev_gg_surv2$table+
  plot_layout(ncol=2,nrow = 2,heights = c(3,1))
total_p3<-ggarrange(mice_line_p1,PC_weight/plot_spacer(),total_p3,ncol=3,widths = c(0.65,0.45,1))
total_p1<-ggarrange(box_p,netVisual_b2+netVisual_b1,widths = c(1,1.5))


total_p<-ggarrange(total_p1,total_p3,nrow = 2,heights = c(2.5,4))
pdf("14.Figure/5.Figure_5.pdf",width = 12,height = 6.5)
print(total_p)
dev.off()
