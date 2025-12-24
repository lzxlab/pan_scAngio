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
library(GSVA)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")







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
#merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
#merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,1/3)]<-"Mid"
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,2/3)]<-"High"
merged_matrix<-merged_matrix[which(merged_matrix$gene_level!="Mid"),]
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))


write.table(merged_matrix,"source_data_NC/Figure_S12A_urothelial.txt",sep = "\t",row.names = F,col.names = T)


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
#merged_matrix$gene_level[merged_matrix$gene>res.cut$cutpoint$cutpoint]<-"High"
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,1/3)]<-"Mid"
merged_matrix$gene_level[merged_matrix$gene>quantile(merged_matrix$gene,2/3)]<-"High"
merged_matrix<-merged_matrix[which(merged_matrix$gene_level!="Mid"),]
merged_matrix$gene_level<-factor(merged_matrix$gene_level,levels = c("Low","High"))

write.table(merged_matrix[,c(1,7,8,11)],"source_data_NC/Figure_S12A_glioblastoma.txt",sep = "\t",row.names = F,col.names = T)


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





library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/CT26growth_volume.csv",sep=",",header = T)

g1<-tab$X17[which(tab$group=="igg")]
g2<-tab$X17[which(tab$group=="anti-vegfr2")]
g3<-tab$X17[which(tab$group=="anti-mcam")]
g4<-tab$X17[which(tab$group=="combine")]

t.test(g1,g3)
t.test(g1,g2)
t.test(g1,g4)
t.test(g2,g4)

color_list<-c("#E41A1C","#377EB8","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))


##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


tab1$Group<-"IgG"
tab1$Group[tab1$group=="anti-vegfr2"]<-"αVEGFR2"
tab1$Group[tab1$group=="anti-mcam"]<-"MCAM-ADC"
tab1$Group[tab1$group=="combine"]<-"αVEGFR2+MCAM-ADC"
tab1$Group<-factor(tab1$Group,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))

library(dplyr)
result_table<-dplyr::summarize(group_by(tab1,Group,Time),
                               mean_var=mean(Value),
                               SD=sd(Value),
                               N_N=n(),  
                               se=SD/sqrt(N_N),  
                               upper_limit=mean_var+SD,  
                               lower_limit=mean_var-SD)



CT26_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
  geom_point(data = result_table,aes(x=Time, y=mean_var),pch=15,size=1.5)+
  geom_line(size=0.5)+
  geom_errorbar(data = result_table,aes(ymin = lower_limit, ymax=upper_limit), #误差条表示均值±标准差
                width=0.2, #误差条末端短横线的宽度
                alpha = 0.7,
                size=0.5) +
  scale_color_manual(values = color_list)+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve (CT26)")+
  theme_classic2()+
  theme(legend.position=c(0.35,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5),breaks =seq(0,2.5,0.5))
print (CT26_line_p1)



write.table(result_table,"source_data_NC/Figure_S12B_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(5,4,3)],"source_data_NC/Figure_S12B_data.txt",sep = "\t",row.names = F,col.names = T)





##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/CT26_weight.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αV"
select_table$Type[select_table$group=="anti-mcam"]<-"M"
select_table$Type[select_table$group=="combine"]<-"αV+M"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "M" ,"αV", "αV+M"))

##Nodes

select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
CT26_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αV", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 2,size=3,
                     comparisons = list(c("IgG" , "M"),c("M", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 2.2,size=3,
                     comparisons = list(c("IgG" ,"αV+M")))+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3"))+
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights (CT26)")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5),breaks =seq(0,2.4,0.4))

print(CT26_weight)


write.table(select_table[,c(3,2)],"source_data_NC/Figure_S12B_weight.txt",sep = "\t",row.names = F,col.names = T)







library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/LLCgrowth_volume.tsv",sep="\t",header = T)


g1<-tab$X17[which(tab$group=="igg")]
g2<-tab$X17[which(tab$group=="anti-vegfr2")]
g3<-tab$X17[which(tab$group=="anti-mcam")]
g4<-tab$X17[which(tab$group=="combine")]

t.test(g1,g3)
t.test(g1,g2)
t.test(g1,g4)
t.test(g2,g4)

color_list<-c("#E41A1C","#377EB8","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))


##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


tab1$Group<-"IgG"
tab1$Group[tab1$group=="anti-vegfr2"]<-"αVEGFR2"
tab1$Group[tab1$group=="anti-mcam"]<-"MCAM-ADC"
tab1$Group[tab1$group=="combine"]<-"αVEGFR2+MCAM-ADC"
tab1$Group<-factor(tab1$Group,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))

library(dplyr)
result_table<-dplyr::summarize(group_by(tab1,Group,Time),
                               mean_var=mean(Value),
                               SD=sd(Value),
                               N_N=n(),  
                               se=SD/sqrt(N_N),  
                               upper_limit=mean_var+SD,  
                               lower_limit=mean_var-SD)



LLC_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
  geom_point(data = result_table,aes(x=Time, y=mean_var),pch=15,size=1.5)+
  geom_line(size=0.5)+
  geom_errorbar(data = result_table,aes(ymin = lower_limit, ymax=upper_limit), #误差条表示均值±标准差
                width=0.2, #误差条末端短横线的宽度
                alpha = 0.7,
                size=0.5) +
  scale_color_manual(values = color_list)+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve (LLC)")+
  theme_classic2()+
  theme(legend.position=c(0.35,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6),breaks =seq(0,1.6,0.4))
print (LLC_line_p1)



write.table(result_table,"source_data_NC/Figure_S12C_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(5,4,3)],"source_data_NC/Figure_S12C_data.txt",sep = "\t",row.names = F,col.names = T)





##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/LLCgrowth_weight.tsv",sep="\t",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αV"
select_table$Type[select_table$group=="anti-mcam"]<-"M"
select_table$Type[select_table$group=="combine"]<-"αV+M"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "M" ,"αV", "αV+M"))

##Nodes

select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
LLC_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 0.8,size=3,
                     comparisons = list(c("αV", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.2,size=3,
                     comparisons = list(c("IgG" , "M"),c("M", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.4,size=3,
                     comparisons = list(c("IgG" ,"αV+M")))+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3"))+
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights (LLC)")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6),breaks =seq(0,1.6,0.4))

print(LLC_weight)


write.table(select_table[,c(3,2)],"source_data_NC/Figure_S12C_weight.txt",sep = "\t",row.names = F,col.names = T)










library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/RENCAgrowth_volume.tsv",sep="\t",header = T)

g1<-tab$X16[which(tab$group=="igg")]
g2<-tab$X16[which(tab$group=="anti-vegfr2")]
g3<-tab$X16[which(tab$group=="anti-mcam")]
g4<-tab$X16[which(tab$group=="combine")]

t.test(g1,g3)
t.test(g1,g2)
t.test(g1,g4)
t.test(g2,g4)

color_list<-c("#E41A1C","#377EB8","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))


##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


tab1$Group<-"IgG"
tab1$Group[tab1$group=="anti-vegfr2"]<-"αVEGFR2"
tab1$Group[tab1$group=="anti-mcam"]<-"MCAM-ADC"
tab1$Group[tab1$group=="combine"]<-"αVEGFR2+MCAM-ADC"
tab1$Group<-factor(tab1$Group,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))

library(dplyr)
result_table<-dplyr::summarize(group_by(tab1,Group,Time),
                               mean_var=mean(Value),
                               SD=sd(Value),
                               N_N=n(),  
                               se=SD/sqrt(N_N),  
                               upper_limit=mean_var+SD,  
                               lower_limit=mean_var-SD)



RENCA_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
  geom_point(data = result_table,aes(x=Time, y=mean_var),pch=15,size=1.5)+
  geom_line(size=0.5)+
  geom_errorbar(data = result_table,aes(ymin = lower_limit, ymax=upper_limit), #误差条表示均值±标准差
                width=0.2, #误差条末端短横线的宽度
                alpha = 0.7,
                size=0.5) +
  scale_color_manual(values = color_list)+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve (RENCA)")+
  theme_classic2()+
  theme(legend.position=c(0.35,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5),breaks =seq(0,1.5,0.5))
print (RENCA_line_p1)



write.table(result_table,"source_data_NC/Figure_S12D_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(5,4,3)],"source_data_NC/Figure_S12D_data.txt",sep = "\t",row.names = F,col.names = T)





##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/RENCAgrowth_weight.tsv",sep="\t",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αV"
select_table$Type[select_table$group=="anti-mcam"]<-"M"
select_table$Type[select_table$group=="combine"]<-"αV+M"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "M" ,"αV", "αV+M"))

##Nodes

select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
RENCA_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αV", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 2,size=3,
                     comparisons = list(c("IgG" , "M"),c("M", "αV+M")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 2.2,size=3,
                     comparisons = list(c("IgG" ,"αV+M")))+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3"))+
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights (RENCA)")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5),breaks =seq(0,2.4,0.4))

print(RENCA_weight)


write.table(select_table[,c(3,2)],"source_data_NC/Figure_S12D_weight.txt",sep = "\t",row.names = F,col.names = T)





total_p1<-Bev_gg_surv3$plot+Bev_gg_surv4$plot+
  Bev_gg_surv3$table+Bev_gg_surv4$table+
  plot_layout(ncol=2,nrow = 2,heights = c(3,1))


total_p2<-ggarrange(CT26_line_p1+CT26_weight+plot_layout(ncol = 2,widths = c(1.5,1)),ncol = 2,widths = c(2,0.6))
total_p3<-ggarrange(LLC_line_p1+LLC_weight+plot_layout(ncol = 2,widths = c(1.5,1)),ncol = 2,widths = c(2,0.6))
total_p4<-ggarrange(RENCA_line_p1+RENCA_weight+plot_layout(ncol = 2,widths = c(1.5,1)),ncol = 2,widths = c(2,0.6))
total_p<-ggarrange(total_p1,total_p2,total_p3,total_p4,ncol = 1,nrow = 4,heights = c(5,3,3,3))

pdf("14.Figure/12.FigureS12.pdf",width = 9,height = 12)
print(total_p)
dev.off()