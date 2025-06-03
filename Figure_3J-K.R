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




##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/Angiogenesis_plot.csv",sep=",",header=T)
select_table$Type<-as.character(select_table$group)


##Nodes
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
select_table$mg.ml<-select_table$Nb.nodes
  
data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.nodes<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+

  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Number of nodes',title=paste ("Number of nodes")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,1500))

print(Nb.nodes)






##Nb.Junctions
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
select_table$mg.ml<-select_table$Nb.Junctions

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Nb.Junctions<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Number of junctions',title=paste ("Number of junctions")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,500))
  

print(Nb.Junctions)







##Tot.meshes.area
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot.meshes.area[which(select_table$Type=="normoxia")])
select_table$mg.ml<-select_table$Tot.meshes.area/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot.meshes.area<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Area (fold change)',title=paste ("Total meshes area")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))
  

print(Tot.meshes.area)




##Tot.meshes.area
select_table$Type<-factor(select_table$Type,levels=c("normoxia","hypoxia"))
mean_var<-mean(select_table$Tot..lenght[which(select_table$Type=="normoxia")])
select_table$mg.ml<-select_table$Tot..lenght/mean_var

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Tot..lenght<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#4DBBD5", "#E64B35"))+
  labs(x = 'Type', y = 'Length (fold change)',title=paste ("Total length")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,1.4))
  

print(Tot..lenght)



mid_p<-Nb.nodes+Nb.Junctions+Tot..lenght+Tot.meshes.area+plot_layout(ncol=2)
mid_p


library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/pericyte_tumor_curve.csv",sep=",",header = T)

g1<-tab$X19[which(tab$group=="HCT116+Pericyte")]
g2<-tab$X19[which(tab$group=="HCT116")]
color_list<-c("#E41A1C","#377EB8","#984EA3",  "#4DAF4A" , "#FF7F00" ,"#A65628","#F781BF" ,"#999999")


#tab<-tab[,c(1,3:7)]
library(reshape2)
tab1<-melt(tab,value.name = "Value",id.vars = c("group"))
tab1$Time<-as.numeric(gsub("X","",tab1$variable))
tab1$Group<-"HCT116"
tab1$Group[tab1$group=="HCT116+Pericyte"]<-"HCT116+PCs"
tab1$Group<-factor(tab1$Group,levels=c("HCT116+PCs" , "HCT116"))
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
  scale_color_manual(values = color_list)+
  #geom_hline(yintercept = 1,size=0.5)+
  # ylim(c(0,3))+
  ylab("Tumor volume (cm3)")+
  xlab("Time (Days)")+
  ggtitle("Tumor growth curve\n(HCT116 with or without PCs)")+
  theme_classic2()+
  theme(legend.position=c(0.3,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8),breaks =seq(0,0.8,0.2))
print (mice_line_p1)






##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/pericyte_tumoe_weight.csv",sep=",",header=T)
select_table$Type<-"HCT116"
select_table$Type[select_table$group=="HCT116+Pericyte"]<-"HCT116+PCs"
select_table$Type<-factor(select_table$Type,levels=c("HCT116+PCs" , "HCT116"))


##Nodes
select_table$Type<-factor(select_table$Type,levels=c( "HCT116","HCT116+PCs"))
select_table$mg.ml<-select_table$weight

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
PC_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c( "#377EB8", "#E41A1C"))+
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
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8))

print(PC_weight)




select_table<-read.table("Experiments/pericyte_CD31.csv",sep=",",header=T)
select_table$Type<-"HCT116"
select_table$Type[select_table$group=="HCT116+Pericyte"]<-"HCT116+PCs"
select_table$Type<-factor(select_table$Type,levels=c("HCT116+PCs" , "HCT116"))


##Nodes
select_table$Type<-factor(select_table$Type,levels=c( "HCT116","HCT116+PCs"))
select_table$mg.ml<-select_table$X.area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
PC_CD31<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", method = "t.test",label.y.npc = 0.85,size=3)+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c( "#377EB8", "#E41A1C"))+
  labs(x = 'Type', y = 'CD31 density',title=paste ("CD31 density")) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))

print(PC_CD31)




mid_p<-Nb.nodes+Nb.Junctions+Tot..lenght+Tot.meshes.area+plot_layout(ncol=2)
right_p<-PC_weight+PC_CD31+plot_layout(ncol=1)
mid_p


total_p1<-ggarrange(mid_p,mice_line_p1,right_p,ncol=4,widths = c(2,1.5,1.4,0.7))

pdf("14.Figure/Figure_3_part2.pdf",width = 12,height = 3)
print(total_p1)
dev.off()
