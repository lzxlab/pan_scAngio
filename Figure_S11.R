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
library(monocle)
library(ggpubr)
library(patchwork)
library(cowplot)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/18.pan_Endo/PGF/")



##Stack plot

sub_sce<-readRDS("3.Cluster/13.Annotation/3.PC_annotation_new.rds")
sub_sce$Type[which(sub_sce$Type!="PN")]<-"Tumor"
sub_sce$Type[which(sub_sce$Type=="PN")]<-"Normal"

dir.create("3.Cluster/13.plot/1.total/")

combined<-sub_sce

##All cell type
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,Sample_ID=combined$Sample_ID,num=1)
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

sum_table1<-sum_table[sum_table$Cluster %in% c("PC_C4_imPC-MCAM", "PC_C5-prolif"),]
sum_table1<-dplyr::summarize(group_by(sum_table1,Type,Sample_ID),per=sum(per))

sum_table1<-sum_table1[order(sum_table1$per,decreasing = T),]
Sample_ID<-unique(c(sum_table1$Sample_ID,sum_table$Sample_ID))
sum_table$Sample_ID<-factor(sum_table$Sample_ID,levels=Sample_ID)

sum_table$Cluster<-factor(sum_table$Cluster,levels = c("PC_C1_mPC-ACTG2", "PC_C2_mPC-MYH11" ,"PC_C3_imPC-CD36",
                                                       "PC_C4_imPC-MCAM", "PC_C5-prolif"))


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


sum_table$per<-round(sum_table$per,2)
write.table(sum_table,"source_data_NC/Figure_S11A.txt",sep = "\t",row.names = F,col.names = T)



test<-readRDS("5.monocle2/5.CAF/monocle2.RDS")
test$Type[which(test$Type!="PN")]<-"Tumor"
test$Type[which(test$Type=="PN")]<-"Normal"

g.colSet1 <-RColorBrewer::brewer.pal(6,"Set2")[c(1:3)]
trajectory_p1<-plot_cell_trajectory(test,color_by = "Pseudotime",show_branch_points = F,cell_size = 0.3)+
  theme(legend.position = "right")+
  ggtitle("Cell trajectory by pseudotime")
trajectory_p2<-plot_cell_trajectory(test,color_by = "Cluster",show_branch_points = F,cell_size = 0.3)+
  scale_color_manual(values = g.colSet1)+
  theme(legend.position = "none")+
  ggtitle("Cell trajectory by cluster")
trajectory_p3<-plot_cell_trajectory(test,color_by = "Cluster",show_branch_points = F,cell_size = 0.3)+
  scale_color_manual(values = g.colSet1)+
  facet_wrap(~Type,nrow=2)+theme(legend.position = "right")


sum_table<-trajectory_p1$data
sum_table<-data.frame(sum_table[,c("sample_name"  , "data_dim_1" ,  "data_dim_2"  ,          
                                   "sample_state","Sample.type","Cluster")],Pseudotime=pData(test)$Pseudotime)

sum_table$data_dim_1<-round(sum_table$data_dim_1,2)
sum_table$data_dim_2<-round(sum_table$data_dim_2,2)
sum_table$Pseudotime<-round(sum_table$Pseudotime,2)
write.table(sum_table,"source_data_NC/Figure_S11B-C.txt",sep = "\t",row.names = F,col.names = T)



total_p1<-trajectory_p1+trajectory_p2+trajectory_p3+plot_layout(ncol = 3,widths = c(2,2,1))



pdf("14.Figure/Figure_S11A-C",width = 10,height = 3)
print(total_p1)
dev.off()





Enrich_type<-"Reactome"
Type="MCAM+ imPC"
Select_cluster_list<-"PVL_C4_imPC-EGFL6"
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/3.PVL/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
  enrich_table$Cluster<-cluster
  total_table<-rbind(total_table,enrich_table)
  enriched_pathway<-c(enriched_pathway,enrich_table$Description)
}
enriched_pathway<-c("Extracellular matrix organization","Signaling by Rho GTPases","Neutrophil degranulation","ECM proteoglycans","RHO GTPase Effectors","Signaling by Interleukins","Collagen formation","Post-translational protein phosphorylation","Interferon Signaling","Integrin cell surface interactions","Signaling by NOTCH","MAPK family signaling cascades","Signaling by MET")

total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% enriched_pathway),]


S2<- ggplot(total_table1, aes(x= reorder(Description, Count), y=Count,fill=p.adjust)) +
  geom_bar(stat="identity",width=0.8)+
  theme_classic()+
  scale_fill_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  xlab("Pathway")+
  ggtitle(paste(Type," (",Enrich_type,")",sep=""))+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        legend.position = c(0.75,0.3),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))+
  coord_flip()
S2

total_table1$Cluster<-"MCAM+ imPC"
write.table(total_table1,"source_data_NC/Figure_S11D.txt",sep = "\t",row.names = F,col.names = T)



markers<-c("NOTCH3","HEYL","JAG1","PSMB9","NEURL1B","HEY2","PSMB8","CCND1","PSMB10","PSME2","PSME1","HIF1A","PSMA5","STAT1","RBX1")

combined<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")


Idents(combined)<-combined$TopCluster
combined<-subset(combined,idents="PVL")

combined$SubCluster<-factor(combined$SubCluster,levels=rev(c("PVL_C1_mPC-ACTG2", 
                                                             "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36" , "PVL_C4_imPC-MCAM", "PVL_C5-prolif" )))
Idents(combined)<-combined$SubCluster
new.cluster.id<-rev(c(  "mPC-ACTG2" ,
                        "mPC-MYH11" , "imPC-CD36" , 
                        "imPC-MCAM" ,"prolif") )


library(viridis)
names(new.cluster.id)<-levels(combined)
combined<-RenameIdents(combined,new.cluster.id)
combined$SubCluster1<-Idents(combined)

Idents(combined)<-combined$SubCluster1

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Vln_exp1<-VlnPlot(combined,features = markers,pt.size=0,group.by="SubCluster1",stack = T) +
  ggtitle("Gene expression of Notch pathway")+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 0,vjust = 0.5,hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))
Vln_exp1

endo_data<-round(t(data.frame(combined@assays$RNA@scale.data[c("NOTCH3","HEYL","JAG1","PSMB9","NEURL1B","HEY2","PSMB8","CCND1","PSMB10","PSME2","PSME1","HIF1A","PSMA5","STAT1","RBX1"),])),2)
endo_data<-data.frame(Cluster=combined$SubCluster,endo_data)
endo_data$Cluster<-gsub("PVL","PC",endo_data$Cluster)
write.table(endo_data,"source_data_NC/Figure_S11E.txt",sep = "\t",row.names = F,col.names = T)



library(CellChat)
##Cell_chat
cellchat<-readRDS("6.cell_chat/2.All/cellchat.RDS")
#cellchat@idents<-factor(cellchat@idents,levels=c( "B" ,"Plasma"   ,    "CD4_conv"   ,   "CD4_Treg",  "CD8", "MAIT"  ,   "NK"   ,    "DC"  , "Macro"  ,    "Mono"  ,  
#                                                  "Mast" ,    "Endo", "APC","pAD"  ,  "SMC"            )) 

groupSize <- as.numeric(table(cellchat@idents))
cellchat@idents<-factor(cellchat@idents,levels=c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                        "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi")
)
levels(cellchat@idents)
vertex.receiver = seq(1,4)


pathways.show <- "VEGF"
select_list<-"VEGF|PTN|IL|CSF|ANGPT|DLL|JAG|NOTCH"
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])
#pairLR.use$interaction_name<-unique(sort(pairLR.use$interaction_name))
pairLR.use<-pairLR.use[!(grepl("CD99|THBS1|HLA",pairLR.use$interaction_name)),]

pairLR.use<-data.frame(interaction_name=pairLR.use)
netVisual_b1<-netVisual_bubble(cellchat, pairLR.use = pairLR.use, targets.use = 14,sources.use = c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                                                                          "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"),remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "none",
        legend.text = element_text(size = 8))

cell_chat_p2_data<-netVisual_b1$data
cell_chat_p2_data$source<-as.character(cell_chat_p2_data$source)
cell_chat_p2_data$target<-as.character(cell_chat_p2_data$target)
cell_chat_p2_data$source.target<-as.character(cell_chat_p2_data$source.target)
cell_chat_p2_data$source[which(cell_chat_p2_data$source=="PVL")]<-"PC"
cell_chat_p2_data$target[which(cell_chat_p2_data$target=="PVL")]<-"PC"
cell_chat_p2_data$source.target<-gsub("PVL","PC",cell_chat_p2_data$source.target)
#cell_chat_p2_data$group.names<-gsub("PVL","PC",cell_chat_p2_data$group.names)

write.table(cell_chat_p2_data,"source_data_NC/Figure_S11F.txt",sep = "\t",row.names = F,col.names = T)





library(CellChat)
library(patchwork)

cellchat1<-readRDS("6.cell_chat/2.All/PN.RDS")
cellchat2<-readRDS("6.cell_chat/2.All/PT.RDS")



cellchat <- mergeCellChat(list(cellchat1, cellchat2), add.names = c("PN", "PT"))

cellchat@idents<-factor(cellchat@idents,levels=c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                        "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi")
)

groupSize <- as.numeric(table(cellchat1@idents))
levels(cellchat1@idents) 
vertex.receiver = seq(1,4) 




pathways.show <- "VEGF"
select_list<-"VEGF|PTN|IL|CSF|ANGPT|DLL|JAG|NOTCH"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])
pairLR.use<-pairLR.use[!(grepl("CD99|THBS1|HLA",pairLR.use$interaction_name)),]

pairLR.use<-data.frame(interaction_name=pairLR.use)

cell_chat_p2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,  comparison = c(1, 2) , targets.use = 14, sources.use = c(       "B"  ,   "Plasma",  "CD4_conv" ,"CD4_Treg", "CD8"  ,"MAIT",   "NK" ,   
                                                                                                                                   "DC" , "Mono"   , "Macro" ,"Mast" ,"NF", "CAF"   ,   "PVL",     "Endo" ,   "Epi"), remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8))

cell_chat_p2_data<-cell_chat_p2$data
cell_chat_p2_data$source<-as.character(cell_chat_p2_data$source)
cell_chat_p2_data$target<-as.character(cell_chat_p2_data$target)
cell_chat_p2_data$source.target<-as.character(cell_chat_p2_data$source.target)
cell_chat_p2_data$source[which(cell_chat_p2_data$source=="PVL")]<-"PC"
cell_chat_p2_data$target[which(cell_chat_p2_data$target=="PVL")]<-"PC"
cell_chat_p2_data$source.target<-gsub("PVL","PC",cell_chat_p2_data$source.target)
cell_chat_p2_data$group.names<-gsub("PVL","PC",cell_chat_p2_data$group.names)

write.table(cell_chat_p2_data,"source_data_NC/Figure_S11G.txt",sep = "\t",row.names = F,col.names = T)




##Cluster
combined1<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")
Idents(combined1)<-combined1$TopCluster
combined<-subset(combined1,idents="PVL")

combined$SubCluster<-factor(combined$SubCluster,levels=rev(c("PVL_C1_mPC-ACTG2", 
                                                             "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36" , "PVL_C4_imPC-MCAM", "PVL_C5-prolif" )))
Idents(combined)<-combined$SubCluster
new.cluster.id<-rev(c(  "mPC-ACTG2" ,
                        "mPC-MYH11" , "imPC-CD36" , 
                        "imPC-MCAM" ,"prolif") )


library(viridis)
names(new.cluster.id)<-levels(combined)
combined<-RenameIdents(combined,new.cluster.id)
combined$SubCluster1<-Idents(combined)

Idents(combined)<-combined$SubCluster1



Gene_name1="Hypoxia"
Gene_name2="PGF+ Tip cell markers"

DefaultAssay(combined)<-"RNA"

RNA_matrix<-combined@assays$RNA@counts


mymatrix<-as.matrix(RNA_matrix)
Tip_Markers<-c( "ADM","AK4","BNIP3","CA9","CCNG2","ENO1","HK2","LDHA","PFKFB3","PGK1","SLC2A1","VEGFA","PDGFB","PGF","CXCL12","KITLG","ANGPT2")

Endo_Markers<-c("ESM1","PGF","APLN","PXDN","TP53I11","IL32","LXN","LOX","ACTB","ACTG1",
                "CXCR4","ANGPT2","IGF2","PECAM1","VWF")
mysymbol1<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)
mysymbol2<-data.frame(Gene_set="Endo",Gene_symbol=Endo_Markers)
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



select_fpkm_matrix<-data.frame(Sample_ID=combined$Sample_ID,
                               Type=combined$SubCluster1,
                               gene=mt)


select_fpkm_matrix$Type<-factor(select_fpkm_matrix$Type,levels=c("mPC-ACTG2" ,
                                                                 "mPC-MYH11" , "imPC-CD36" , 
                                                                 "imPC-MCAM" ,"prolif"))



box_p1<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_violin(aes(fill=Type),trim=FALSE,color="white") + 
  geom_boxplot(aes(fill=Type),width=0.2,position=position_dodge(0.9),color="black")+
  #scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  ggtitle("Hypoxia score")+
  ylab("Hypoxia score")+
  stat_compare_means(ref.group = 1,label = "p.format",size=3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  ylim(c(0.5,1.8))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))
box_p1

select_fpkm_matrix$gene<-round(select_fpkm_matrix$gene,2)
write.table(select_fpkm_matrix,"source_data_NC/Figure_S11H.txt",sep = "\t",row.names = F,col.names = T)




combined<-readRDS("3.Cluster/13.Annotation/2.Stromal_annotation.rds")


Idents(combined)<-combined$TopCluster
combined<-subset(combined,idents="PVL")

combined$SubCluster<-factor(combined$SubCluster,levels=rev(c("PVL_C1_mPC-ACTG2", 
                                                             "PVL_C2_mPC-MYH11" , "PVL_C3_imPC-CD36" , "PVL_C4_imPC-MCAM", "PVL_C5-prolif" )))
Idents(combined)<-combined$SubCluster
new.cluster.id<-rev(c(  "mPC-ACTG2" ,
                        "mPC-MYH11" , "imPC-CD36" , 
                        "imPC-MCAM" ,"prolif") )


library(viridis)
names(new.cluster.id)<-levels(combined)
combined<-RenameIdents(combined,new.cluster.id)
combined$SubCluster1<-Idents(combined)

Idents(combined)<-combined$SubCluster1

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Idents(combined)<-combined$Type
combined<-subset(combined,idents=c("PT","PN"))

Vln_exp2<-VlnPlot(combined,features = Tip_Markers,pt.size=0,group.by="SubCluster1",stack = T,split.by = "Type",
                  cols = c("#B31B21", "#1465AC")) +
  ggtitle("Gene expression of hypoxia pathway (tumor vs. normal)")+
  stat_compare_means(label = "p.format",label.y.npc = 0.9)+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 0,vjust = 0.5,hjust = 0.5),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))
Vln_exp2


endo_data<-round(t(data.frame(combined@assays$RNA@scale.data[Tip_Markers,])),2)
endo_data<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,endo_data)
endo_data$Cluster<-gsub("PVL","PC",endo_data$Cluster)
write.table(endo_data,"source_data_NC/Figure_S11I.txt",sep = "\t",row.names = F,col.names = T)





total_p1<-ggarrange(S2,Vln_exp1,ncol=2,widths = c(1,1.4))
total_p2<-ggarrange(netVisual_b1,cell_chat_p2,ncol=2,widths = c(1,2))
total_p3<-ggarrange(box_p1,Vln_exp2,ncol=2,widths = c(1,3))

total_p<-ggarrange(total_p1,total_p2,total_p3,ncol=1,nrow = 4,
                   heights = c(3,3,3,3))
pdf("14.Figure/7.Figure_7.pdf",width = 12,height = 12)
print(total_p)
dev.off()

