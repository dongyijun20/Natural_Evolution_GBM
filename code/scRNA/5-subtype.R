library(NMF)
library(data.table)
library(ggsci)
library(pheatmap)
library(AUCell)
library(GEOquery)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(corrplot)
library(xlsx)

markers<-as.list(read.xlsx("mmc2.xlsx",sheetIndex = 1,startRow = 4))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1.S","G2.M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)
gene_sets <- markers

# pbmc <- readRDS("result/merged_new.rds")
# pbmc <- FindClusters(pbmc, resolution = 1)
# DimPlot(pbmc, label=T)
# pbmc$tumor_ident_fine <- "malignant"
# pbmc$tumor_ident_fine[pbmc$seurat_clusters%in%c(34,20,8,27,22,36,26,31,13,19,12,16,3)] <- "non-malignant"
# saveRDS(pbmc, file = "result/merged_new.rds")
# 
# pdf("plot_new/tumor_ident_fine.pdf")
# DimPlot(pbmc, group.by = "tumor_ident_fine",cols = c("malignant"="red","non-malignant"="blue"))+
#   ggtitle("Malignant Cells Identification",subtitle = paste("Fraction = ", 100*round(sum(pbmc$tumor_ident_fine=="malignant")/ncol(pbmc),2),"%",sep = ""))
# dev.off()

malignant <- readRDS("result/malignant_withXHsample.rds")

library(AUCell)
library(clusterProfiler)
cells_rankings <- AUCell_buildRankings(malignant@assays$RNA@data) 

NES = list(c('S100A10','FOSL2','SPP1','CAV1','ANXA1','VIM','CD44','SERPINH1',
             'LGALS3','CEBPB','ATF5','LGALS1'))
gene_sets[["NES"]] <- NES[[1]]
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
malignant$AUC_AC  <- as.numeric(getAUC(cells_AUC)["AC", ])
malignant$AUC_MES  <- as.numeric(getAUC(cells_AUC)["MES", ])
malignant$AUC_NPC  <- as.numeric(getAUC(cells_AUC)["NPC", ])
malignant$AUC_OPC  <- as.numeric(getAUC(cells_AUC)["OPC", ])
malignant$AUC_NES  <- as.numeric(getAUC(cells_AUC)["NES", ])

pdf("plot_new/MES_NES_FeaturePlot_withXHsample.pdf", width = 20, height = 5)
FeaturePlot(malignant, features = c("AUC_NES", "AUC_MES"), max.cutoff = "q90", min.cutoff = "q10", blend = T)
dev.off()

df<-data.frame(MES = malignant$AUC_MES, sample = malignant$sample)
df %>% group_by(sample) %>% summarise_at(vars (MES), list (name = mean))
comparelist <- list(c("NJ01_1","NJ01_2"),c("NJ02_1","NJ02_2"),c("TT01_1","TT01_2"),c("TT02_1","TT02_2")) 
pdf("plot_new/MES_NES_AUC_withXHsample.pdf")
ggboxplot(df, x = "sample", y = "MES",
          color = "sample", palette = SampleColor)+ ylab("MES score")+
  stat_compare_means(comparisons = comparelist, label = "p.signif")+NoLegend()
df<-data.frame(NES = malignant$AUC_NES, sample = malignant$sample)
df %>% group_by(sample) %>% summarise_at(vars (NES), list (name = mean))
ggboxplot(df, x = "sample", y = "NES",
          color = "sample", palette = SampleColor)+ ylab("NES score")+
  stat_compare_means(comparisons = comparelist, label = "p.signif")+NoLegend()
dev.off()

library(GSEABase)
library(GSVA)

DimPlot(malignant)
table(malignant$seurat_clusters)
exp=AverageExpression(malignant,assays = "RNA", group.by = "seurat_clusters")
counts2=exp[["RNA"]]
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=gene_sets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=T,
                   scale = "none") #以行来标准化，这个功能很不错

subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
names(subtype) <- gsub("g", "", names(subtype))
malignant$subtype <- as.vector(sapply(malignant$seurat_clusters,function(x) subtype[x]))

library(ArchR)
subtype_color<-ArchRPalettes$stallion[1:length(unique(malignant$subtype))]
names(subtype_color) <- sort(unique(malignant$subtype))

pdf("plot_new/malignant_subtype_withXHsample.pdf")
DimPlot(malignant, group.by = "subtype", cols = subtype_color)
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=T,
                   scale = "column") #以行来标准化，这个功能很不错
dev.off()

library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- malignant@meta.data
P1=CellInfo %>% ggplot(aes(x=grade, fill=fct_rev(subtype))) +scale_fill_manual(values = subtype_color)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(subtype , fill=subtype))+geom_bar(stat="count",colour = "black",width = 0.7)+  scale_fill_manual(values = subtype_color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)

pdf("plot_new/proportion_subtype_withXHsample.pdf",width=7,height=6)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

saveRDS(malignant, file = "result/malignant_withXHsample.rds")

