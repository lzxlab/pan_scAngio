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


library(magrittr)
library(ggsignif)

library(ggsci)
library(RColorBrewer)




Gene_name="VEGFA"
select_gene=Gene_name


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
Tip_Markers<-c("MYL9","GUCY1B3","ACTA2","SOD3","CALD1","TPM2","TAGLN","EGFL6","RGS5","PPP1R14A","SEPT4","HIGD1B","MFGE8","GJA4","FOXS1","C1QTNF1","HEYL","GUCY1A3","COX4I2","KCNMB1","ADAP2","RCSD1","BGN","FSCN1","SPARCL1","CRISPLD2","S100A10","SLC38A11","HEY2","NOTCH3","WFDC1","PTP4A3","TPM1","PERP","GEM","COL4A2","GPR4","CSRP2","KCNA5","C11orf96","MYLK","HSPB6","EFHD1","DKK3","MCAM","S100A16","FRZB","ESAM","RASL12","PLXDC1","EDNRA","FOXF1","A2M")
Endo_Markers<-c("PLVAP","AQP1","PECAM1","RNASE1","CD93","EMCN","SRGN","ADGRL4","FBLN1","GIMAP4","PTPRB","IL3RA","ARHGAP29","PITX2","VWF","ESAM","CDH5","S1PR1","CLEC14A","RAMP3","TIE1","ICAM2","FLT1","CXorf36","CD34","CYYR1","PODXL","FLI1","CCL14","CLDN5","RAMP2","ENG","LDB2","PCDH17","HYAL2","NPDC1","ADCY4","SLCO2A1","TMEM255B","MCTP1","RND1","PLTP","GNG11","CD200","ERG","ACKR1","TGFBR2","COL4A1","COL15A1","EGFL7")
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
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=20)

##survplot
mt<-es.dif[1,]/es.dif[2,]







Gene_name1="PCs"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=es.dif[1,],
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)

select_table<-merge_matrix


#绘图
TCGA_cor_p5 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("ssGSEA score (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of PCs and ECs, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p5)


select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("PCs","ECs")
write.table(select_table,"source_data_NC/Figure_S7B.txt",sep = "\t",row.names = F,col.names = T)




Gene_name1="PGF"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
select_table<-merge_matrix


#绘图
TCGA_cor_p1 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of PGF and ECs, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p1)

select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("PGF","ECs")
write.table(select_table,"source_data_NC/Figure_S9A.txt",sep = "\t",row.names = F,col.names = T)





Gene_name1="ANGPT2"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
select_table<-merge_matrix


#绘图
TCGA_cor_p2 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of ANGPT2 and ECs, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p2)


select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("ANGPT2","ECs")
write.table(select_table,"source_data_NC/Figure_S9B.txt",sep = "\t",row.names = F,col.names = T)






Gene_name1="VEGFA"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)

select_table<-merge_matrix


#绘图
TCGA_cor_p3 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of VEGFA and ECs, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p3)


select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("VEGFA","ECs")
write.table(select_table,"source_data_NC/Figure_S10A.txt",sep = "\t",row.names = F,col.names = T)



Gene_name1="ANGPT2"
Gene_name2="PCs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)

select_table<-merge_matrix


#绘图
TCGA_cor_p6 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of PCs and ANGPT2, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p6)

select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("PGF","PCs")
write.table(select_table,"source_data_NC/Figure_S8B.txt",sep = "\t",row.names = F,col.names = T)


Gene_name1="PGF"
Gene_name2="PCs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)

select_table<-merge_matrix


#绘图
TCGA_cor_p7 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("Expression levels (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of PCs and PGF, TCGA pan-caner dataset (RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = 0.5,colour = "grey"),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p7)


select_table<-select_table[,c("sample","project_id","gene1","gene2")]
colnames(select_table)[c(3,4)]<-c("ANGPT2","PCs")
write.table(select_table,"source_data_NC/Figure_S8A.txt",sep = "\t",row.names = F,col.names = T)


##survplot

Gene_name1="PGF"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix1<-merge_matrix

Gene_name1="VEGFA"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix2<-merge_matrix


project_id_list<-unique(merge_matrix2$project_id)
total_tab<-data.frame()
for (i in 1:length(project_id_list)){
  project_id<-project_id_list[[i]]
  
  
  merged_table<-merge_matrix1[which(merge_matrix1$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  cor.value1<-cor.value
  p.value1<-p.value
  
  merged_table<-merge_matrix2[which(merge_matrix2$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  
  cor.value2<-cor.value
  p.value2<-p.value
  
  tmp_table<-data.frame(Cancer.type=project_id,
                        cor.value=cor.value1,
                        p.value=p.value1,
                        Type="PGF~ECs"
  )
  total_tab<-rbind(total_tab,tmp_table)
  
  tmp_table<-data.frame(Cancer.type=project_id,
                        cor.value=cor.value2,
                        p.value=p.value2,
                        Type="VEGFA~ECs"
  )
  total_tab<-rbind(total_tab,tmp_table)
  
}
total_tab$Type <-factor(total_tab$Type,levels = c("VEGFA~ECs","PGF~ECs"))
box_p<-ggpaired(total_tab, x="Type", y="cor.value", id = "Cancer.type",
                color = "grey70",
              
                add="jitter",line.color = "grey", line.size = 0.5,
                title = "Angiogenic effects",
                #palette=RColorBrewer::brewer.pal(2,"blues")[2:3],
                xlab=" ", 
                ylab="Cor.value", 
                legend.title="Type",show.legend = F) + 
  geom_point(aes(group=Type,color=Cancer.type))+
  theme_classic()+
  #scale_color_manual(values=RColorBrewer::brewer.pal(7,"Set1")[c(1:2)])+
  #geom_text_repel(aes(label=Patient_ID))+
  stat_compare_means(label = "p.format",method = "t.test",
                     paired = T,label.y.npc = 0.9) +#配对t检验
  theme(legend.position = 'right',
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
        axis.text.x = element_text(size = 10,vjust =0.5,hjust = 0.5),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.text  = element_text(size = 10),
        axis.line = element_blank(),
        strip.background = element_blank()
  )+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=4))

print(box_p)



write.table(total_tab,"source_data_NC/Figure_S10B_left.txt",sep = "\t",row.names = F,col.names = T)




Gene_name1="ANGPT2"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[1,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix1<-merge_matrix

Gene_name1="VEGFA"
Gene_name2="ECs"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)
merge_matrix2<-merge_matrix


project_id_list<-unique(merge_matrix2$project_id)
total_tab<-data.frame()
for (i in 1:length(project_id_list)){
  project_id<-project_id_list[[i]]
  
  
  merged_table<-merge_matrix1[which(merge_matrix1$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  cor.value1<-cor.value
  p.value1<-p.value
  
  merged_table<-merge_matrix2[which(merge_matrix2$project_id==project_id),]
  
  cor.value<-cor(merged_table$gene1,merged_table$gene2)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$gene1,merged_table$gene2)$p.value
  p.value<-signif(p.value,digits = 2)
  
  cor.value2<-cor.value
  p.value2<-p.value
  
  tmp_table<-data.frame(Cancer.type=project_id,
                        cor.value=cor.value1,
                        p.value=p.value1,
                        Type="ANGPT2~ECs"
  )
  total_tab<-rbind(total_tab,tmp_table)
  
  tmp_table<-data.frame(Cancer.type=project_id,
                        cor.value=cor.value2,
                        p.value=p.value2,
                        Type="VEGFA~ECs"
  )
  total_tab<-rbind(total_tab,tmp_table)
  
}
total_tab$Type <-factor(total_tab$Type,levels = c("VEGFA~ECs","ANGPT2~ECs"))
box_p1<-ggpaired(total_tab, x="Type", y="cor.value", id = "Cancer.type",
                color = "grey70",
                
                add="jitter",line.color = "grey", line.size = 0.5,
                title = "Angiogenic effects",
                #palette=RColorBrewer::brewer.pal(2,"blues")[2:3],
                xlab=" ", 
                ylab="Cor.value", 
                legend.title="Type",show.legend = F) + 
  geom_point(aes(group=Type,color=Cancer.type))+
  theme_classic()+
  #scale_color_manual(values=RColorBrewer::brewer.pal(7,"Set1")[c(1:2)])+
  #geom_text_repel(aes(label=Patient_ID))+
  stat_compare_means(label = "p.format",method = "t.test",
                     paired = T,label.y.npc = 0.9) +#配对t检验
  theme(legend.position = 'right',
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
        axis.text.x = element_text(size = 10,vjust =0.5,hjust = 0.5),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        strip.text  = element_text(size = 10),
        axis.line = element_blank(),
        strip.background = element_blank()
  )+ 
  guides(colour = guide_legend(override.aes = list(size=3),ncol=4))

print(box_p1)

write.table(total_tab,"source_data_NC/Figure_S10B_right.txt",sep = "\t",row.names = F,col.names = T)




total_p1<-ggarrange(box_p,box_p1,common.legend = T,legend = "right")

pdf("14.Figure/Figure_S5B.pdf",width = 12,height = 8)
print(TCGA_cor_p5)
dev.off()


pdf("14.Figure/Figure_S6.pdf",width = 12,height = 15)
print(TCGA_cor_p7/TCGA_cor_p6)
dev.off()



pdf("14.Figure/Figure_S7.pdf",width = 12,height = 15)
print(TCGA_cor_p2/TCGA_cor_p1)
dev.off()



pdf("14.Figure/Figure_S8.pdf",width = 12,height = 10)
print(TCGA_cor_p3/total_p1+plot_layout(heights = c(7,3)))
dev.off()


