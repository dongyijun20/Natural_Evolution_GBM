### RNA-seq--------------
bulk <- read.csv("data/bulk.csv")
head(bulk)




### single-cell------------
pat <- "WJS-SX"
seu.data <- Read10X(data.dir = paste0("data/N20240174-01-01/0_Cellranger/",pat))
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

### Meta data annotations
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^MT-', col.name = 'percent.mt')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPS', col.name = 'percent.rps')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPL', col.name = 'percent.rpl')
seu_raw$percent.rp <- seu_raw$percent.rps + seu_raw$percent.rpl
seu_raw$barcode <- rownames(seu_raw@meta.data)
seu_raw$barcode_pat <- paste0(rownames(seu_raw@meta.data), '_', pat)

seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^MT-")
x <- PercentageFeatureSet(seu_raw, pattern = "^MT-")
summary(x)

# ### Identify doublets using scrublet
# doublet_rate_tmp <- ncol(seu_raw)*8*1e-6 
# writeMM(seu_raw@assays$RNA@counts, paste0('data/N20240174-01-01/0_Cellranger/', pat, '/matrix_', pat, '_raw.mtx'))
# system(paste('python3 0-scrublet.py', pat, doublet_rate_tmp))
# doublets <- read.table(paste0('data/N20240174-01-01/0_Cellranger/', pat, '/doublets_', pat, '_raw.txt'), header = T)
# seu_raw[['predicted_doublets']] <- doublets$predicted_doublets
# seu_raw[['doublet_scores']] <- doublets$doublet_scores
# #system(paste0('rm data_new/', pat, '/matrix_', pat, '_raw.mtx'))
# #system(paste0('rm data_new/', pat, '/doublets_', pat, '_raw.txt'))

### Seurat workflow
seu_raw <- NormalizeData(seu_raw)
seu_raw <- FindVariableFeatures(seu_raw)
seu_raw <- ScaleData(seu_raw)
seu_raw <- RunPCA(seu_raw)
seu_raw <- RunUMAP(seu_raw, dims = 1:40)
seu_raw <- FindNeighbors(seu_raw, dims = 1:40)
seu_raw <- FindClusters(seu_raw, resolution = 0.1)

### Identify doublets using DF
sweep.list <- paramSweep_v3(seu_raw, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>%
  arrange(desc(BCmetric))
pK <- pK[1, 2]
pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
nExp <- round(doublet_rate_tmp * dim(seu_raw@assays$RNA@counts)[2])
seu_raw <- doubletFinder_v3(seu_raw, PCs = 1:40, pK = pK, nExp = nExp)
seu_raw$doublet <- seu_raw@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)]
seu_raw$DF_score <- seu_raw@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)]

### Filtering steps 
pdf(paste0("plot_new/",pat,"_QC.pdf"))
VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
#table(seu_raw$doublet,seu_raw$predicted_doublets)
DimPlot(seu_raw, cells.highlight = colnames(seu_raw)[seu_raw$doublet=="Doublet"], sizes.highlight = 0.1)
#DimPlot(seu_raw, cells.highlight = colnames(seu_raw)[seu_raw$predicted_doublets==T], sizes.highlight = 0.1)
dev.off()

minFeature <- 200
maxFeature <- 8000
minCount <- 1000
maxCount <- 40000
maxMT <- 20

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature &
                nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & 
                doublet == 'Singlet')
ncol(seu_raw)
ncol(seu)

### Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
seu <- RunUMAP(seu, dims = 1:40)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu)

### Preliminary cell type identification using SingleR
seu_sce <- as.SingleCellExperiment(seu)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels

iced <- DatabaseImmuneCellExpressionData()
pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
pruneScores(pred_iced_main)
seu[['celltype_iced_main']] <- pred_iced_main$pruned.labels

hpca <- HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']] <- pred_hpca_main$pruned.labels

mid <- MonacoImmuneData()
pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
pruneScores(pred_mid_main)
seu[['celltype_mid_main']] <- pred_mid_main$pruned.labels

pdf(paste0("plot_new/",pat,"_basic.pdf"))
DimPlot(seu, label = T)
DimPlot(seu, group.by = "celltype_bped_main", label = T)
DimPlot(seu, group.by = "celltype_iced_main",label = T)
DimPlot(seu, group.by = "celltype_hpca_main",label = T)
DimPlot(seu, group.by = "celltype_mid_main",label = T)
dev.off()

### Save object
seu@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)] <- NULL
seu@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)] <- NULL
saveRDS(seu, file = paste0("result/", pat, ".rds"))


## validation
SX<-readRDS("result/WJS-SX.rds")
HS<-readRDS("result/WJS-HS.rds")
SX$position<-"SX"
HS$position<-"HS"

merge <- merge(SX, c(HS))
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge)
merge <- ScaleData(merge)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge <- CellCycleScoring(merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merge <- RunPCA(merge, npcs = 50)
merge <- merge %>% RunHarmony("position", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution=0.5) %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

SampleColor <- ArchRPalettes$stallion[1:2]
names(SampleColor)<-c("HS","SX")
DimPlot(merge, label = T)
DimPlot(merge, group.by = "position", cols = SampleColor, shuffle = T)
saveRDS(merge,file="result/merge_WJS.rds")

FeaturePlot(merge,c("CD68","CD14","FCGR1A","NUPR1"))

TAM<-subset(merge,seurat_clusters%in%c(9,10))
TAM <- NormalizeData(TAM)
TAM <- FindVariableFeatures(TAM)
TAM <- ScaleData(TAM)
TAM <- RunPCA(TAM, npcs = 50)
TAM <- TAM %>% RunHarmony("position", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution=0.5) %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

DimPlot(TAM, label = T)
DimPlot(TAM, group.by = "position", cols = SampleColor, shuffle = T)

## 先找SX和HS之间的差异基因
Idents(TAM)<-"position"
markers <- FindMarkers(TAM,ident.1 = "HS",ident.2 = "SX",test.use = "wilcox")
# markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
# markers_df <- markers_df[markers_df$p_val_adj<0.01,]
# write.table(markers_df,file="result/myeloid_markers200_grade.txt",quote=F,sep="\t",row.names=F,col.names=T)

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.01)
data$row <- rownames(data)
data$regulate[data$regulate=="Up"] <- "HS"
data$regulate[data$regulate=="Down"] <- "SX"
data$regulate[data$regulate=="Normal"] <- "non-significant"

# plot
pdf("plot_new/volcano_HS_vs_SX.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 0.5, FDR_cut = 0.01, legend_position = "DL",
          colors = c("SX"="navy","HS"="firebrick3","non-significant"="grey"), 
          fills = c("SX"="navy","HS"="firebrick3","non-significant"="grey"), 
          label = "row", label_number = 20, output = FALSE)
dev.off()

markers[rownames(markers)=="NUPR1",]
FeaturePlot(TAM,"NAV3")

library(xlsx)
BMDM_MB_markers <- read.xlsx("reference/13059_2017_1362_MOESM5_ESM.xlsx",sheetIndex = 1)
BMDM_MB_markers <- lapply(1:2, function(x) na.omit(BMDM_MB_markers[,x]))
names(BMDM_MB_markers) <- c("MG","BMDM")
BMDM_MB_markers <- list(BMDM=unique(c(BMDM_MB_markers$BMDM,c("CCR2","CD45RA","CD141","ICAM","CD1C","CD1B",
                                                             "TGFBI","FXYD5","FCGR2B","CLEC12A","CLEC10A","CD207",
                                                             "CD209","CD49D"))),
                        MG=unique(c(BMDM_MB_markers$MG,c("CX3CR1","SALL1","HEXB","P2RY12","TMEM119"))))

### GSVA
Idents(TAM)<-"seurat_clusters"
exp=AverageExpression(TAM,assays = "RNA")
counts2=exp[["RNA"]]
library(GSVA)
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=BMDM_MB_markers, 
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
TAM$TAM_type <- sapply(TAM$seurat_clusters,function(x) subtype[x])
#myeloid$TAM_type <- subtype
DimPlot(TAM,group.by = "TAM_type")

FeaturePlot(TAM,c("IL8","NUPR1","VIM","LDHA","CSTB","PLIN2","S100A10"),min.cutoff = "q5",max.cutoff = "q95",cols = c("lightgrey" ,"#DE1F1F"))
VlnPlot(TAM,c("IL8","NUPR1","VIM","LDHA","CSTB","PLIN2","S100A10"))
VlnPlot(TAM,c("IL8","NUPR1","VIM","LDHA","CSTB","PLIN2","S100A10"),group.by = "position")
saveRDS(TAM,file = "result/TAM_WJS.rds")

## malignant
DimPlot(merge, label = T)
FeaturePlot(merge, features = c("CD8A","CD3E","CDH5","COL1A1","SOX4","AREG"))
malignant <- subset(merge, seurat_clusters%in%c(0:8,11,12,14,15))

malignant <- NormalizeData(malignant)
malignant <- FindVariableFeatures(malignant)
malignant <- ScaleData(malignant)
malignant <- RunPCA(malignant, npcs = 50)
malignant <- malignant %>% RunHarmony("position", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution=0.5) %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

markers<-as.list(read.xlsx("mmc2.xlsx",sheetIndex = 1,startRow = 4))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1.S","G2.M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)
gene_sets <- markers
Idents(malignant)<-"seurat_clusters"
exp=AverageExpression(malignant,assays = "RNA")
counts2=exp[["RNA"]]
library(GSVA)
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=gene_sets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目

pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=F,
                   scale = "column") #以行来标准化，这个功能很不错

subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
malignant$subtype <- sapply(malignant$seurat_clusters,function(x) subtype[x])
malignant
subtype_color<-ArchRPalettes$stallion[1:length(unique(malignant$subtype))]
names(subtype_color) <- sort(unique(malignant$subtype))

pdf("plot_new/malignant_WJS.pdf")
DimPlot(malignant, group.by = "subtype", cols = subtype_color)
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=T,
                   scale = "column") #以行来标准化，这个功能很不错
dev.off()

## crosstalk
merge$celltype <- NA
merge$celltype[match(colnames(malignant),colnames(merge))] <- malignant$subtype
merge$celltype[match(colnames(TAM),colnames(merge))] <- paste(TAM$TAM_type,TAM$seurat_clusters,sep="_")
table(merge$celltype)
GBM <- merge[,!is.na(merge$celltype)]
GBM
merge

#不同分组之间的配对分析
sc.sp=SplitObject(GBM,split.by = "position")
sc.11=GBM[,colnames(sc.sp[["HS"]])]
sc.3=GBM[,colnames(sc.sp[["SX"]])]

library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)

cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

cellchat=cellchat.sc11 
CellChatDB <- CellChatDB.human
cellchat@DB  <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
#################################
cellchat=cellchat.sc3
cellchat@DB  <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
##############################################
source("visualization.R")

cc.list=list(HS=cc.sc11,SX=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
##可视化
##所有细胞群总体观：通讯数量与强度对比
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 

##数量与强度差异网络图
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
##红色是case相对于control上调的，蓝色是下调的

#数量与强度差异热图
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
#case和control对比，红色是上调，蓝色是下调

#保守和特异性信号通路的识别与可视化
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
##左图最下面多个信号通路是case组独有的

##细胞互作数量对比网络图
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,celltype.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,celltype.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )


table(GBM@active.ident)
s.cell=c( "BMDM_1", "MES")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]

netVisual_circle(count1,weight.scale = T,celltype.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(count2,weight.scale = T,celltype.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )

weight.max <- getMaxWeight(cc.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cc.list)) {
  netVisual_circle(cc.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cc.list)[i]))
}

weight.max <- getMaxWeight(cc.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cc.list)) {
  netVisual_circle(cc.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(cc.list)[i]))
}

group.cellType <- c(rep("BMDM", 3), rep("MG", 2), rep("Tumor", 4))
cc.list$HS@meta$celltype <- factor(cc.list$HS@meta$celltype, levels = c("BMDM_1","BMDM_2","BMDM_4", "MG_0","MG_3", "AC","MES","OPC","NPC"))
group.cellType <- factor(group.cellType, levels = c("BMDM", "MG", "Tumor"))
object.list <- lapply(cc.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","weight", "weight.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,8,9),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(1,5,8,9), targets.use = 2,  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,8,9),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in SX", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,8,9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in SX", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "HS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "SX",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "HS",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 2, targets.use = c(1,5,8,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(1,5,8,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(1,5,8,9), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(1,5,8,9), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))

pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("HS", "SX")) # set factor level
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.

saveRDS(cellchat, file = "result/cellchat_comparisonAnalysis_SX_vs_HS.rds")


## score
library(ggpubr)
library(rstatix)

pdf("plot_new/relationship_NESBMDM_Ferro_WJS.pdf")
NES_BMDM = list(c("IL8", "NUPR1","VIM","LDHA","CSTB","PLIN2","S100A10"))

TAM = Seurat::AddModuleScore(TAM, features = NES_BMDM, name='NES_BMDM')
FeaturePlot(TAM, features = 'NES_BMDM1', reduction = 'umap', max.cutoff = "q95", min.cutoff = "q5",cols = c("lightgrey" ,"#DE1F1F"))
df <- TAM@meta.data
stat.test <- df %>% t_test(NES_BMDM1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
ggboxplot(df, x = "seurat_clusters", y = "NES_BMDM1", fill = "seurat_clusters") + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)

Ferroptosis= list(c('ACSL1',	'ACSL3',	'ACSL4',	'ACSL5',	'ACSL6',	'ALOX15',
                    'ATG5',	'ATG7',	'CP',	'CYBB',	'FTH1',	'FTL',	'FTMT',	'GCLC',
                    'GCLM',	'GPX4',	'GSS',	'HMOX1',	'LPCAT3',	'MAP1LC3A',	'MAP1LC3B',
                    'MAP1LC3C',	'NCOA4',	'PCBP1',	'PCBP2',	'PRNP',	'SAT1',	'SAT2',	'SLC11A2',
                    'SLC39A14',	'SLC39A8',	'SLC3A2',	'SLC40A1',	'SLC7A11',	'STEAP3',	'TF',	'TFRC',
                    'TP53',	'VDAC2',	'VDAC3'))
TAM = Seurat::AddModuleScore(TAM, features = Ferroptosis, name='Ferroptosis')
FeaturePlot(TAM, features = 'Ferroptosis1', reduction = 'umap', max.cutoff = "q95", min.cutoff = "q5",cols = c("lightgrey" ,"#DE1F1F"))+ 
  ggtitle("Ferroptosis", subtitle = paste0("correlation with NES_BMDM score = ",round(as.numeric(cor.test(TAM$Ferroptosis1,TAM$NES_BMDM1)$estimate),2)))
FeaturePlot(TAM, features = c('Ferroptosis1',"NES_BMDM1"), reduction = 'umap', max.cutoff = "q95", min.cutoff = "q5",blend = T)

df <- TAM@meta.data
stat.test <- df %>% t_test(Ferroptosis1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
ggboxplot(df, x = "seurat_clusters", y = "Ferroptosis1", fill = "seurat_clusters") + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)
dev.off()
