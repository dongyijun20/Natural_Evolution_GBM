myeloid <- readRDS("result/myeloid.rds")
myeloid <- subset(myeloid, TAM_type=="MG")
table(myeloid$batch)

myeloid <- myeloid %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

myeloid <- myeloid %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.5)
myeloid <- myeloid %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(myeloid)
table(myeloid$seurat_clusters)

saveRDS(myeloid, file = "result/MG.rds")

Idents(myeloid) <- "seurat_clusters"
markers <- FindAllMarkers(myeloid, only.pos = TRUE, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
write.table(markers_df,file="result/MG_markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

#markers_df <- read.table("result/myeloid_markers200.txt", header = T)
markers_df <- markers_df[-grep("^MT",markers_df$gene),]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#intersect(markers_df$gene[markers_df$cluster=="4"],data$row[data$regulate=="older"])


pdf("plot_new/MG_heatmap.pdf")
DoHeatmap(subset(myeloid, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
dev.off()
pdf("plot_new/MG_dotplot.pdf", height = 8, width = 5)
DotPlot(myeloid, features = unique(markers_df$gene)) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 10))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

cluster_color<-c('#4b6aa8','#3ca0cf','#c376a7')
DimPlot(myeloid, label=T, cols = cluster_color)
SampleColor <- readRDS("color/SampleColor.rds")
GradeColor <- readRDS("color/GradeColor.rds")

pdf("plot_new/MG_basics.pdf", width = 4, height = 3)
DimPlot(myeloid, label=T, cols = cluster_color)
DimPlot(myeloid, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(myeloid, group.by = "grade", cols = GradeColor)
dev.off()

table(myeloid$sample)
library(forcats)
library(ggplot2)
library(gridExtra)

CellInfo <- myeloid@meta.data
P1=CellInfo %>% ggplot(aes(x=sample, fill=seurat_clusters)) +scale_fill_manual(values = cluster_color)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(seurat_clusters , fill=seurat_clusters))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  scale_fill_manual(values = cluster_color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)

P3=CellInfo %>% ggplot(aes(x=grade, fill=TAM_type)) +scale_fill_manual(values = col.ls)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P3g=ggplotGrob(P3)
P4=ggplot(CellInfo, aes(TAM_type , fill=TAM_type))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  scale_fill_manual(values = col.ls)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P4g=ggplotGrob(P4)
pdf("plot_new/proportion_MG.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
grid.arrange(grobs=list(P3g,P4g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

sce.markers <- read.table("result/MG_markers200.txt", header = T)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
names(gcSample)

## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG", 
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.y = element_text(size=8),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))

