setwd("~/data/single_cell/18.pan_Endo/PGF/")
library(ggplot2)
library(ggpubr)
library(DOSE)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathview)
library(topGO)
library(msigdbr)
library(ggrepel)
library(RColorBrewer)

##NC_vs_KD
##siMARCO_vs_NC
siMARCO_vs_NC<-read.table("~/data/lzx_project/21.PC_CD46/5.DEG/Pos_vs_Neg.edgeR.xls",sep="\t",header = T,row.names = 1)
siMARCO_vs_NC$gene<-rownames(siMARCO_vs_NC)

trans = bitr(siMARCO_vs_NC$gene, fromType="SYMBOL", toType=c("ENTREZID", "GENENAME"), OrgDb="org.Hs.eg.db")
siMARCO_vs_NC<-merge(siMARCO_vs_NC,trans,by.x="gene",by.y="SYMBOL")
siMARCO_vs_NC<-siMARCO_vs_NC[order(siMARCO_vs_NC$logFC,decreasing = T),]

FCgenelist<-siMARCO_vs_NC$logFC
names(FCgenelist)<-siMARCO_vs_NC$ENTREZID
egseKEGG <- gseKEGG(FCgenelist,organism="hsa",
                    pvalueCutoff = 1,
                    use_internal_data = T) # 默认organism="hsa"\


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene,gene_symbol)
m_t2g$gene_name<-m_t2g$gene_symbol
egseKEGG <- GSEA(FCgenelist,TERM2GENE = m_t2g,
                 pvalueCutoff = 1) # 默认organism="hsa"\



egseKEGG[,c('Description','enrichmentScore')]
egseKEGG[1:5,c('Description','enrichmentScore')]


pathway_list<-c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYOGENESIS","HALLMARK_ANGIOGENESIS","HALLMARK_NOTCH_SIGNALING")
up_pathway<-egseKEGG@result[pathway_list,]
up_pathway$Type<-"Up"

pathway_list<-c("HALLMARK_KRAS_SIGNALING_DN","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_UV_RESPONSE_DN","HALLMARK_TGF_BETA_SIGNALING")
down_pathway<-egseKEGG@result[pathway_list,]
down_pathway$Type<-"Down"
total_table<-rbind(up_pathway,down_pathway)
total_table$text_y<- 0
total_table$pos_y<- 1
total_table$pos_y[which(total_table$Type=="Down")]<- 0
total_table$Description<-gsub("HALLMARK_","",total_table$Description)

write.table(total_table,"source_data_NC/Figure_6B.txt",sep = "\t",row.names = F,col.names = T)

GSEA_P1<-ggplot(total_table, aes(x=reorder(Description, NES), y=NES,fill=Type))+
  geom_bar(stat="identity",width=0.8,alpha=0.6)+
  scale_fill_manual(values = c("#336699", "#993399"))+
  geom_text(aes(x=reorder(Description, NES),y=text_y,hjust =pos_y,label = Description),size=2.5)+
  theme_minimal()+
  ggtitle(paste("Hallmarks enrichment (MCAM+ vs. MCAM- PCs)",sep="-")) +
  theme()+
  xlab("Pathways")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+
  coord_flip()
print(GSEA_P1)




gseaScores <- getFromNamespace("gseaScores", "DOSE")

# define function
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}



GSEA_obj<-egseKEGG
pathway<-'HALLMARK_ANGIOGENESIS'
gene_list<-siMARCO_vs_NC$gene


gsdata <- gsInfo(GSEA_obj, geneSetID = pathway)


gsdata1 <- gsdata %>%
  mutate("gene_name" =gene_list ) %>%
  filter(position == 1)


color <- rev(colorRampPalette(c("#336699","white", "#993399"))(10))

# plot
pcurve <- ggplot(gsdata,aes(x = x,y = runningScore,color = runningScore)) +
  geom_hline(yintercept = 0,size = 1,color = 'black',
             lty = 'dashed') +
  geom_line() +
  # geom_segment(data = gsdata1,aes(xend = x,yend = 0)) +
  theme_bw(base_size = 10) +
  scale_color_gradient(low = '#76BA99',high = '#EB4747') +
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.margin = margin(t = .2,r = .2, b = 0,l = .2,unit = "cm")) +
  ylab('Running Enrichment Score')

pcurve


pseg <- ggplot(gsdata,aes(x = x,y = runningScore)) +
  geom_segment(data = gsdata1,
               aes(x = x,xend = x,y = 0,yend = 1),
               color = 'black',show.legend = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 10) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(t = 0,r = .2,b = .2,l = .2,unit = "cm")) +
  xlab('Rank in Ordered Dataset')

pseg


v <- seq(1, sum(gsdata$position), length.out = 9)
inv <- findInterval(rev(cumsum(gsdata$position)), v)
if (min(inv) == 0) inv <- inv + 1

col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

# ymin <- min(p2$data$ymin)
# yy <- max(p2$data$ymax - p2$data$ymin) * .3
ymin <- 0
yy <- 0.3
xmin <- which(!duplicated(inv))
xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
d <- data.frame(ymin = ymin, ymax = yy,
                xmin = xmin,
                xmax = xmax,
                col = col[unique(inv)])

pseg_ht <- pseg + geom_rect(
  aes_(xmin = ~xmin,xmax = ~xmax,
       ymin = ~ymin,ymax = ~ymax,
       fill = ~I(col)),
  data = d,
  alpha = 0.8,
  inherit.aes = FALSE)

pseg_ht


# add gene rank
pseg_ht1 <- pseg_ht + xlab('') +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(t = -.1,r = .2,b = 0,l = .2,unit = "cm"))

prank <- ggplot(gsdata,aes(x = x,y = geneList)) +
  geom_col(width = 1,fill = 'grey80',color = NA) +
  # geom_col(aes(fill = geneList),
  #          width = 1,color = NA,show.legend = F) +
  # scale_fill_gradient2(low = col[1],mid = 'white',high = col[length(col)],midpoint = 0) +
  geom_hline(yintercept = 0,size = 0.8,color = 'black',
             lty = 'dashed') +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(t = -.1,r = .2,b = .2,l = .2,unit = "cm")) +
  coord_cartesian(expand = 0) +
  ylab('Ranked List') +
  xlab('Rank in Ordered Dataset')

prank

write.table(prank$data,"source_data_NC/Figure_6C.txt",sep = "\t",row.names = F,col.names = T)


pall <- aplot::plot_list(gglist = list(pcurve,pseg_ht1,prank),
                         ncol = 1, heights = c(0.5,0.2,0.3))

pall

#geneLabel <- gsdata1 %>% arrange(desc(runningScore)) %>%
#  head(20)
geneLabel <-gsdata1[which(gsdata1$gene_name %in% c("")),]


plabel <- pcurve +
  geom_segment(data = geneLabel,aes(xend = x,yend = 0),
               color = 'red') +
  geom_text_repel(data = geneLabel,
                  aes(label = gene_name),
                  force = 20,
                  max.overlaps = 50,
                  # nudge_y = 0.2,
                  size = 4,
                  fontface = 'italic')

aplot::plot_list(gglist = list(plabel,pseg_ht1,prank),
                 ncol = 1, heights = c(0.5,0.2,0.3))

panother <- ggplot(gsdata,aes(x = x,y = runningScore,color = runningScore)) +
  geom_hline(yintercept = 0,size = 0.8,color = 'black',
             lty = 'dashed') +
  geom_point() +
  geom_line() +
  geom_segment(data = gsdata1,aes(xend = x,yend = 0)) +
  theme_bw(base_size = 10) +
  # scale_color_gradient(low = '#336699',high = '#993399') +
  scale_color_gradient2(low = '#336699',mid = 'white',high = '#993399',midpoint = 0.2) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.margin = margin(t = .2,r = .2, b = 0,l = .2,unit = "cm")) +
  ylab('Running Enrichment Score')

panother


panother_label <-
  panother +
  ggtitle(pathway)+
  geom_text_repel(data = geneLabel,
                  aes(label = gene_name),
                  min.segment.length = 0,
                  force = 20,
                  max.overlaps = 50,
                  # nudge_y = 0.15,
                  size = 4,
                  color = 'black',
                  fontface = 'italic')

panother_label

ht <- ggplot(gsdata,aes(x = x,y = runningScore)) +
  geom_rect(aes_(xmin = ~xmin,xmax = ~xmax,
                 ymin = ~ymin,ymax = ~ymax,
                 fill = ~I(color)),
            data = d,
            alpha = 0.8,
            inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 0,r = .2, b = .2,l = .2,unit = "cm"))

# combine
GSEA_P2<-aplot::plot_list(gglist = list(panother_label,ht),
                 ncol = 1, heights = c(0.9,0.1))


vocano_plot = function(sorted_matrix, Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05,gene_list=NA){
  par(mar = c(5, 6, 5, 5))
  sorted_matrix$PValue[sorted_matrix$PValue<10^-20]<-10^-20
  tab = data.frame(logFC = as.numeric(as.character(sorted_matrix$logFC)), 
                   negLogPval = -log10(as.numeric(as.character(sorted_matrix$PValue))))
  rownames(tab)=rownames(sorted_matrix)
  tab <- na.omit(tab)
  nosigGene = rownames(tab)[(abs(tab$logFC) <= lfc | tab$negLogPval <= -log10(pval))]
  sigGenes_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
  sigGenes_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
  draw_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
  draw_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
  up_count = length(sigGenes_up)
  down_count = length(sigGenes_down)
  nosig_count = length(nosigGene)
  gap = max(tab$logFC)/50
  FCrange=ceiling(max(abs(tab$logFC)))
  tab[sigGenes_up,"SigGenes"]=paste("1.up:",up_count,sep="")
  tab[sigGenes_down,"SigGenes"]=paste("2.down:",down_count,sep="")
  tab[nosigGene,"SigGenes"]=paste("3.noSig:",nosig_count,sep="")
  tab$name=rownames(tab)
  options(stringsAsFactors = FALSE)  ### NOTICE!!!
  DF=data.frame(name=as.character(tab$name),SigGenes=as.factor(tab$SigGenes),logFC=tab$logFC,negLogPval=tab$negLogPval)
  rownames(DF)=rownames(tab)
  #DF <- DF[sort(DF$logFC,index.return=TRUE, decreasing = TRUE)$ix,]
  tophit=DF[c(draw_up[1:10],draw_down[max((length(draw_down)-10), 0):length(draw_down)]),]
  xmax <- ceiling(max(abs(DF$logFC)))
  ymax <- ceiling(max(abs(DF$negLogPval)))*1.1
  
  DF$Gene<-""
  DF$Gene[rownames(DF) %in% gene_list]<-rownames(DF)[rownames(DF) %in% gene_list]
  DF<-DF[order(DF$logFC,decreasing = T),]
  DF<-DF[order(DF$negLogPval,decreasing = T),]
  
  DF$Gene[DF$logFC>lfc][1:3]<-rownames(DF)[DF$logFC>lfc][1:3]
  
  DF<-DF[order(DF$logFC,decreasing = F),]
  DF<-DF[order(DF$negLogPval,decreasing = T),]
  DF$Gene[DF$logFC< -lfc][1:5]<-rownames(DF)[DF$logFC< -lfc][1:5]
  #DF$Gene[which(rownames(DF)=="ERBB2")]<-"ERBB2"
  
  p <- ggplot(DF, aes(x = logFC, y = negLogPval, label=DF$Gene)) +
    geom_point(aes(color = SigGenes))+ xlim(-xmax,xmax) + ylim(0,ymax) +
    geom_text_repel(max.overlaps = 1000,size=2.5)+
    scale_color_manual(values = c("#B31B21", "#1465AC","grey")) +
    theme_bw(base_size = 10) + theme(legend.position = c(0.16,0.7)) +
    xlab("Log2FoldChange")+ylab(paste("-log10","PValue",sep=""))+
    geom_vline(aes(xintercept=-lfc),colour="darkgrey", linetype="dashed")+
    geom_vline(aes(xintercept=lfc),colour="darkgrey", linetype="dashed") +
    geom_hline(aes(yintercept=-log10(pval)),colour="darkgrey", linetype="dashed")+
    ggtitle(paste(Sample_2,Sample_1,sep=paste(rep(" ",12),collapse=""))) +
    #expression("Group 1" %->% "Group 2"),
    annotate("text", x=-xmax*0.9, y=-log10(pval), label= paste("PValue","<",pval,sep=""))+
    annotate("text", x=0, y=-log10(pval), label= "2fold")+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}

libs <- c("limma","ggpubr","ggrepel","reshape2","cluster","splines","ggplot2","gridExtra")
loaded <- sapply(libs, library, character.only = T)
fc = 2
lfc = log2(fc)
pval = 0.05
res<-read.table("~/data/lzx_project/21.PC_CD46/5.DEG/Pos_vs_Neg.edgeR.xls",sep="\t",header = T,row.names = 1)

sorted_ori_de <- res[sort(res$logFC,index.return=TRUE, decreasing = F)$ix,]
sorted_ori_de <- sorted_ori_de[is.finite(sorted_ori_de$logFC),]

vocano_p<-vocano_plot(sorted_ori_de, Sample_1 = 'MCAM+ PCs', Sample_2 = 'MCAM- PCs', 
            lfc = lfc, pval = pval,gene_list=c("VEGFA","RGS5","PGF","ANGPT2","HIGD1B","PDGFRA","MCAM","NOTCH3"))


out_p<-data.frame(gene_name=rownames(sorted_ori_de),sorted_ori_de)
write.table(out_p,"source_data_NC/Figure_6A.txt",sep = "\t",row.names = F,col.names = T)




library(ggplot2)
library(ggpubr)
tab<-read.table("Experiments/E0771growth_volume.csv",sep=",",header = T)

g1<-tab$X17[which(tab$group=="igg")]
g2<-tab$X17[which(tab$group=="anti-vegfr2")]
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



E0771_line_p1 <- ggplot(result_table, aes(x=Time, y=mean_var, color=Group,group=Group))+
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
  ggtitle("Tumor growth curve (E0771)")+
  theme_classic2()+
  theme(legend.position=c(0.35,0.75),
        plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2),breaks =seq(0,2,0.4))
print (E0771_line_p1)



write.table(result_table,"source_data_NC/Figure_6E_curve.txt",sep = "\t",row.names = F,col.names = T)
write.table(tab1[,c(5,4,3)],"source_data_NC/Figure_6E_data.txt",sep = "\t",row.names = F,col.names = T)





##
color_list<-c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3", "#FF7F00" ,"#A65628","#F781BF" ,"#999999")

select_table<-read.table("Experiments/E0771_weight.csv",sep=",",header=T)
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
E0771_weight<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.2,size=3,
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
  labs(x = 'Type', y = 'Tumor weights (g)',title=paste ("Terminal tumor weights (E0771)")) +
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

print(E0771_weight)


write.table(select_table[,c(9,2)],"source_data_NC/Figure_6E_weight.txt",sep = "\t",row.names = F,col.names = T)



select_table<-read.table("Experiments/E0771_red_CD31_plot.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))


##Nodes
select_table$mg.ml<-select_table$X.Area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
E0771_CD31<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αVEGFR2", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3,size=3,
                     comparisons = list(c("IgG" , "MCAM-ADC"),c("MCAM-ADC", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3.5,size=3,
                     comparisons = list(c("IgG" ,"αVEGFR2+MCAM-ADC")))+
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
  labs(x = 'Type', y = 'CD31 density (%)',title=paste ("CD31 density (E0771)")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,4))

print(E0771_CD31)

out_p1<-data.frame(Gene="CD31",select_table[,c(3,1)])



select_table<-read.table("Experiments/E0771_green_MCAM_plot.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))


##Nodes
select_table$mg.ml<-select_table$X.Area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
E0771_MCAM<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αVEGFR2", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3,size=3,
                     comparisons = list(c("IgG" , "MCAM-ADC"),c("MCAM-ADC", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3.5,size=3,
                     comparisons = list(c("IgG" ,"αVEGFR2+MCAM-ADC")))+
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
  labs(x = 'Type', y = 'MCAM density (%)',title=paste ("MCAM density (E0771)")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,4))

print(E0771_MCAM)


out_p2<-data.frame(Gene="MCAM",select_table[,c(3,1)])
out_p<-rbind(out_p1,out_p2)
write.table(out_p,"source_data_NC/Figure_6F.txt",sep = "\t",row.names = F,col.names = T)




select_table<-read.table("Experiments/CT26_red_CD31_plot.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))


##Nodes
select_table$mg.ml<-select_table$X.Area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
CT26_CD31<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αVEGFR2", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3,size=3,
                     comparisons = list(c("IgG" , "MCAM-ADC"),c("MCAM-ADC", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3.5,size=3,
                     comparisons = list(c("IgG" ,"αVEGFR2+MCAM-ADC")))+
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
  labs(x = 'Type', y = 'CD31 density (%)',title=paste ("CD31 density (CT26)")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,4))

print(CT26_CD31)
out_p1<-data.frame(Gene="CD31",select_table[,c(3,1)])


select_table<-read.table("Experiments/CT26_green_MCAM_plot.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))


##Nodes
select_table$mg.ml<-select_table$X.Area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
CT26_MCAM<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αVEGFR2", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3,size=3,
                     comparisons = list(c("IgG" , "MCAM-ADC"),c("MCAM-ADC", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3.5,size=3,
                     comparisons = list(c("IgG" ,"αVEGFR2+MCAM-ADC")))+
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
  labs(x = 'Type', y = 'MCAM density (%)',title=paste ("MCAM density (CT26)")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,4))

print(CT26_MCAM)

out_p2<-data.frame(Gene="MCAM",select_table[,c(3,1)])
out_p<-rbind(out_p1,out_p2)
write.table(out_p,"source_data_NC/Figure_6G.txt",sep = "\t",row.names = F,col.names = T)



select_table<-read.table("Experiments/Angpt2_E0771_Angpt2_red.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))


##Nodes
select_table$mg.ml<-select_table$X.Area

data1<-dplyr::summarize(group_by(select_table,Type),
                        mean.var=mean(mg.ml),
                        sd=sd(mg.ml),
                        lower.ci=mean(mg.ml)-sd(mg.ml),
                        upper.ci=mean(mg.ml)+sd(mg.ml)
)
data1$mg.ml<-data1$mean.var
Angpt2_Angpt2<-ggplot(data=select_table, aes(Type, mg.ml, color = Type))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 1.5,size=3,
                     comparisons = list(c("αVEGFR2", "αVEGFR2+MCAM-ADC"),c("αVEGFR2", "MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3,size=3,
                     comparisons = list(c("IgG" , "MCAM-ADC"),c("MCAM-ADC", "αVEGFR2+MCAM-ADC")))+
  stat_compare_means(data=select_table,label = "p.format", 
                     method = "t.test",label.y = 3.5,size=3,
                     comparisons = list(c("IgG" ,"αVEGFR2+MCAM-ADC")))+
  geom_jitter(width = 0.2 )+
  
  geom_bar(data=data1,aes(x=Type, y=mg.ml,color = Type),stat = "identity",size = 1.1,width = 0.7,fill=NA)+
  geom_errorbar(data=data1,aes(x=Type,ymin = lower.ci, ymax=upper.ci),stat = "identity", #误差条表示均值±标准差
                width=0.10, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="black",
                alpha = 0.7,
                size=0.7) +
  theme_classic()+
  
  scale_color_manual(values =c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3"))+
  labs(x = 'Type', y = 'ANGPT2 density (%)',title=paste ("ANGPT2 density (E0771)")) +
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
  scale_y_continuous(expand = c(0,0),limits = c(0,4))

print(Angpt2_Angpt2)



select_table<-read.table("Experiments/Angpt2_E0771_Angpt2_red.csv",sep=",",header=T)
select_table$Type<-"IgG"
select_table$Type[select_table$group=="anti-vegfr2"]<-"αVEGFR2"
select_table$Type[select_table$group=="anti-mcam"]<-"MCAM-ADC"
select_table$Type[select_table$group=="combine"]<-"αVEGFR2+MCAM-ADC"
select_table$Type<-factor(select_table$Type,levels=c("IgG" , "MCAM-ADC" ,"αVEGFR2", "αVEGFR2+MCAM-ADC"))
colnames(select_table)[1]<-"Angpt2"

select_table2<-read.table("Experiments/Angpt2_E0771_Mcam_green.csv",sep=",",header=T)
colnames(select_table2)[1]<-"Mcam"
select_table<-cbind(select_table,select_table2)
select_table<-data.frame(select_table)

cor_p <- ggplot(select_table, aes(x=Angpt2, y=Mcam,color=Type,group = 1))+
  geom_point(data = select_table,aes(x=Angpt2, y=Mcam),pch=15,position=position_dodge(0),size=2)+
  geom_smooth(size=0.5,method="lm",se=F)+
  stat_cor(label.x.npc = 0.1,label.y.npc = 0.9)+
  scale_color_manual(values =c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3"))+
  #scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
  #ylim(c(0,25))+
  #xlim(c(0,25))+
  ylab("MCAM density")+
  xlab("ANGPT2 density")+
  #expand_limits(y = c(0,100))+
  ggtitle(paste(paste("Correlation of ANGPT2 and MCAM")))+
  #labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
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

print(cor_p)


write.table(select_table[,c(3,1,4)],"source_data_NC/Figure_6H.txt",sep = "\t",row.names = F,col.names = T)



library(patchwork)
total_p1<-ggarrange(vocano_p,GSEA_P1,GSEA_P2,ncol=3)
total_p2<-ggarrange(E0771_line_p1+E0771_weight+plot_layout(ncol = 2,widths = c(1.5,1)),ncol = 2,widths = c(2,1.7))
total_p3<-E0771_CD31+E0771_MCAM+cor_p+Angpt2_Angpt2+CT26_CD31+CT26_MCAM+plot_layout(ncol=4,widths = c(1,1,1.5,1.5))
total_p<-ggarrange(total_p1,total_p2,total_p3,ncol = 1,nrow = 3,heights = c(3.2,3,4))
pdf("14.Figure/6.Figure_6.pdf",width = 14,height = 10.7)
print(total_p)
dev.off()
