library(scCATCH)
markers1 <- list("Astrocyte"= c("SPARCL1", "ETNPPL", "CXCL14", "MT1G", "GFAP"),
                 "Microglia (MG)"= c("RNA5SP151", "CCL4L1", "CH25H", "CX3CR1", "P2RY12", "P2RY13"),
                 "Dendritic Cell"= c("ACSL1", "LYZ", "TK1", "HLA-DQA1", "HLA-DPB1","CH25H", "RNA5SP151", "CCL4L1"),
                 "Vascular"= c("APOD", "DCN", "IFITM1", "SPARCL1", "RGS5", "PDGFRB", "NOTCH3"),
                 "Neuron"= c("QDRR", "CCK", "CXCL14", "FABP7", "PTPRZ1"),
                 "T Cell"= c("CD3D", "RNASE1", "CD3E", "GZMK", "PTPRC", "CD45", "CD4", "CD8A"),
                 "B Cell"= c("IGHG1", "IGHG3", "CD79A", "CD19"),
                 "Naïve T Cell"= c("ACSL1", "RNA5SP151"),
                 "Monocyte-derived Macrophage (MDM)"= c("ACSL1", "RNA5SP151", "APOC1", "CD163", "F13A1", "CD14", "CD74"),
                 "Oligodendrocyte"= c("MAG", "QDRR", "PLLP", "APOD", "RNASE1", "OLIG2", "MBP"),
                 "Fibroblast"= c("COL1A1"),
                 "Endothelial Cell"= c("CDH5", "CLDN5", "VWF", "CD34", "PECAM"),
                 "Neutrophil"= c("IL1R2", "CXCR2", "FPR2"))

markers2 <- read.csv("reference/免疫基因表.csv", header = F)
rownames(markers2)<-markers2$V1
markers2<-markers2[,-1]
markers <- list()
for(i in 1:nrow(markers2)) {             # Using for-loop to add columns to list
  markers[[i]] <- markers2[i,][markers2[i,]!=""]
}
names(markers)<-rownames(markers2)
markers[["Microglial cell"]] <- markers1[["Microglial Cell"]]

### human GBM immune atlas
atlas <- read.xlsx("reference/Human_ND_GBM_Full.xlsx", sheetIndex = 1)
atlas <- as.list(atlas)
atlas <- lapply(atlas, na.omit)

markers<-atlas
custom_marker <- data.frame(species = rep("Human",length(unlist(markers))),
                            tissue = rep("Brain",length(unlist(markers))),
                            cancer= rep("Normal",length(unlist(markers))),
                            condition = rep("Normal cell",length(unlist(markers))),
                            subtype1 = NA,
                            subtype2 = NA,
                            subtype3 = NA,
                            celltype = as.character(unlist(sapply(names(markers),function(x) rep(x,length(markers[[x]]))))),
                            gene = as.character(unlist(markers)),
                            resource = NA,
                            pmid = NA, stringsAsFactors = FALSE)

nonmalignant <- readRDS("result/nonmalignant_new.rds")
nonmalignant_subset <- nonmalignant[,sample(colnames(nonmalignant),2000)]
saveRDS(nonmalignant_subset, file="result/nonmalignant_subset.rds")

nonmalignant[["RNA"]] <- JoinLayers(nonmalignant[["RNA"]])
Layers(nonmalignant[["RNA"]])

## run harmony----------
nonmalignant$batch <- nonmalignant$sample

nonmalignant$batch <- sub("(.*)\\_.*","\\1",nonmalignant$sample)

# nonmalignant$batch <- "NJ"
# nonmalignant$batch[nonmalignant$sample=="TT01_1"|nonmalignant$sample=="TT01_2"|nonmalignant$sample=="TT02_1"|nonmalignant$sample=="TT02_2"] <- "TT"

table(nonmalignant$batch)

nonmalignant <- nonmalignant %>% 
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

nonmalignant <- nonmalignant %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 
nonmalignant <- nonmalignant %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(nonmalignant)

library(RColorBrewer)
names<-c("NJ02_1","NJ02_2","TT01_1","TT01_2","TT02_1","TT02_2","NJ01_1","NJ01_2","XH01_1","XH01_2")
SampleColor <- brewer.pal(8, 'Paired')[c(2,1,4,3,6,5,8,7,10,9)]
names(SampleColor)<-names


DimPlot(nonmalignant, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(nonmalignant, group.by = "grade", cols = GradeColor, shuffle = T)
saveRDS(SampleColor, file="color/SampleColor.rds")

GradeColor <- brewer.pal(10, 'Paired')[9:10]
names(GradeColor) <- c("younger","younger")
saveRDS(GradeColor, file="color/GradeColor.rds")

GradeColor <- readRDS("color/GradeColor.rds")
SampleColor <- readRDS("color/SampleColor.rds")

sc_data <- GetAssayData(nonmalignant, slot="data")
sc_cluster <- as.character(nonmalignant@meta.data$seurat_clusters)
obj <- createscCATCH(data = sc_data, cluster = sc_cluster)
obj <- findmarkergene(obj, species = "Human", marker = custom_marker, tissue = c("Brain"), 
                      if_use_custom_marker = T, cancer = F, use_method = "1")
obj <- findcelltype(obj)
table(obj@celltype$cell_type)

celltype = data.frame(ClusterID=obj@celltype$cluster, celltype=obj@celltype$cell_type, stringsAsFactors = FALSE)
nonmalignant@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  nonmalignant@meta.data[which(nonmalignant@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(nonmalignant$celltype)

### cell cycle analysis (some cells may be dividing)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
nonmalignant <- CellCycleScoring(nonmalignant, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(nonmalignant[[]])
pdf("plot_new/nonmalignant_cell_cycle.pdf")
DimPlot(nonmalignant, group.by = "Phase")
dev.off()

macrophage <- c("CCR2","CD45RA","CD141","ICAM","CD1C","CD1B","TGFBI","FXYD5","FCGR2B","CLEC12A","CLEC10A","CD207","CD209","CD49D")
microglial <- c("CX3CR1","SALL1","HEXB","P2RY12","TMEM119")
nonmalignant <- AddModuleScore(nonmalignant, list(macrophage), name = "macrophage")
nonmalignant <- AddModuleScore(nonmalignant, list(microglial), name = "microglial")
FeaturePlot(nonmalignant, features = "macrophage1", max.cutoff = "q99", min.cutoff = "q1")
FeaturePlot(nonmalignant, features = "microglial1", max.cutoff = "q99", min.cutoff = "q1")

pdf("plot_new/nonmalignant_VlnPlot.pdf")
VlnPlot(nonmalignant, features = "macrophage1", group.by = "seurat_clusters")
VlnPlot(nonmalignant, features = "microglial1", group.by = "seurat_clusters")
DotPlot(nonmalignant, features = c(macrophage,microglial)) + RotatedAxis()
dev.off()

### FindMarkers
markers <- FindAllMarkers(nonmalignant, only.pos = TRUE, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
write.table(markers_df,file="result/nonmalignant_markers500_withXHsample.txt",quote=F,sep="\t",row.names=F,col.names=T)

markers_df <- read.table("result/nonmalignant_markers500_withXHsample.txt", header = T)
#write.csv(file="test.csv",rbind(markers_df$cluster[markers_df$cluster==1],markers_df$gene[markers_df$cluster==1]),quote = F,col.names = F,row.names = F)
markers_df = markers_df %>% 
  dplyr::filter(!grepl("^RP", gene)) %>% 
  group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
markers_df$gene[markers_df$cluster==5]

### heatmap------------
# Single cell heatmap of feature expression
Idents(nonmalignant) <- "seurat_clusters"
pdf("plot_new/nonmalignant_heatmap_withXHsample.pdf")
DoHeatmap(subset(nonmalignant, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
dev.off()
pdf("plot_new/nonmalignant_dotplot_withXHsample.pdf",width = 6,height = 10)
DotPlot(nonmalignant, features = unique(markers_df$gene)) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 6))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

## featureplot----------
pdf("plot_new/nonmalignant_featureplot_withXHsample.pdf",width = 10,height = 10)
FeaturePlot(nonmalignant, features=c("CD3E","CD8A", # T cells
                                     "CD68","FCGR1A", # Myeloid cells
                                     "CDH5", # Endothelial cells
                                     "COL1A1", # Fibroblasts
                                     "SOX4", # Tumor cells
                                     "AREG","CD1C" # DC
                                     ), min.cutoff = "q5", max.cutoff = "q95",cols = c("lightgrey" ,"#DE1F1F"))
dev.off()

# ### CellID------------
# library(CelliD)
# library(Seurat)
# library(tidyverse)
# 
# #获取原始表达矩阵
# BaronMatrix <- GetAssayData(nonmalignant, slot="counts")
# #仅考虑编码蛋白质的基因
# data("HgProteinCodingGenes")
# BaronMatrixProt <- BaronMatrix[rownames(BaronMatrix) %in% HgProteinCodingGenes,]
# 
# #这几步类似Seurat的标准流程
# Baron <- CreateSeuratObject(counts = BaronMatrixProt, project = "Baron", min.cells = 5)
# Baron <- NormalizeData(Baron)
# Baron <- ScaleData(Baron, features = rownames(Baron)) 
# Baron <- RunMCA(Baron) #该软件的主要分析函数，将细胞和基因同时降维到一个空间，离细胞近的基因被定义为细胞的特征基因
# 
# #下载参考基因集
# panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
# #根据器官过滤，示例数据就是胰腺的
# panglao_pancreas <- panglao %>% filter(organ == "Brain"|organ == "Connective tissue"|organ == "Immune system")
# #物种包含人
# panglao_pancreas <- panglao_pancreas %>%  filter(str_detect(species,"Hs"))
# #下面两步会得到一个列表，列表的每一个元素是由基因集构成的字符向量，且每个元素被命名为对应的细胞类型的名字
# panglao_pancreas <- panglao_pancreas %>%  
#   group_by(`cell type`) %>%  
#   summarise(geneset = list(`official gene symbol`))
# pancreas_gs <- setNames(panglao_pancreas$geneset, panglao_pancreas$`cell type`)
# 
# HGT_pancreas_gs <- RunCellHGT(Baron, pathways = pancreas_gs, dims = 1:50, n.features = 200)
# pancreas_gs_prediction <- rownames(HGT_pancreas_gs)[apply(HGT_pancreas_gs, 2, which.max)]
# #矩阵的每一列依次：判断是否是最大值，得到一列布尔值，结合矩阵的行名会返回行名中的一个元素（也就是最有可能的细胞类型）。
# #所有列运行完了之后会得到所有细胞最可能的注释
# pancreas_gs_prediction_signif <- ifelse(apply(HGT_pancreas_gs, 2, max)>2, yes = pancreas_gs_prediction, "unassigned")
# #如果`-log10 corrected p-value`的值小于等于2，则认为不显著，注释更正为"unassigned"
# Baron$pancreas_gs_prediction <- pancreas_gs_prediction_signif
# rownames(Baron@meta.data) <- colnames(Baron)
# 
# Baron$batch <- nonmalignant$sample
# Baron <- FindVariableFeatures(Baron, selection.method = "vst", nfeatures = 2000)
# Baron <- ScaleData(Baron)
# Baron <- RunPCA(Baron, npcs = 50, verbose = FALSE)
# Baron <- Baron %>% 
#   RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
# 
# Baron <- Baron %>%
#   FindNeighbors(reduction = "harmony") %>%
#   FindClusters(resolution = 0.5) 
# Baron <- Baron %>%
#   RunUMAP(dims = 1:40, reduction = "harmony")
# 
# cell.type <- sapply(unique(Baron$seurat_clusters), function(x) {
#   k <- table(Baron$pancreas_gs_prediction[Baron$seurat_clusters==x])
#   p <- k["unassigned"]/sum(k)
#   message(p)
#   if(p<0.3|is.na(p)){
#     k <- k[names(k)!="unassigned"]
#     return(names(k)[which.max(k)])
#   } else(return("unassigned"))
# })
# names(cell.type)<-unique(Baron$seurat_clusters)
# 
# Baron$cell.type <- cell.type[Baron$seurat_clusters]
# pdf("plot_new/Cellid.pdf")
# DimPlot(Baron, group.by = "cell.type", pt.size=0.5, label = TRUE,repel = TRUE)
# DimPlot(Baron, group.by = "pancreas_gs_prediction", pt.size=0.5, label = TRUE,repel = TRUE)
# dev.off()

### SingleR--------------
library(celldex)
library(SingleR)
counts <- Read10X("reference/filtered_feature_bc_matrix_HumanNewlyDiagnGBM_Full/filtered_feature_bc_matrix/")
seu_raw <- CreateSeuratObject(counts = counts, project = "Full", min.cells = 1, min.features = 1)
meta <- read.csv(file = "reference/annot_Human_ND_GBM_Full.csv")
table(meta$cluster)
meta$cell <- gsub("_","-",meta$cell)
rownames(meta)<-meta$cell
ref2 <- seu_raw[,meta$cell]
ref2$cluster <- meta$cluster
Idents(ref2)<-"cluster"
ref2 <- subset(ref2, cells = WhichCells(ref2, downsample = 200))
ref2$cluster <- meta[colnames(ref2),]$cluster
table(ref2$cluster)
saveRDS(ref2, file = "reference/GBM_atlas_myeloid_reference.rds")

counts <- Read10X("reference/filtered_feature_bc_matrix_HumanNewlyDiagnGBM_TAM/filtered_feature_bc_matrix/")
seu_raw <- CreateSeuratObject(counts = counts, project = "TAM", min.cells = 1, min.features = 1)
meta <- read.csv(file = "reference/annot_Human_ND_GBM_TAM.csv")
table(meta$cluster)
meta$cell <- gsub("_","-",meta$cell)
rownames(meta)<-meta$cell
ref3 <- seu_raw[,meta$cell]
ref3$cluster <- meta$cluster
Idents(ref3)<-"cluster"
ref3 <- subset(ref3, cells = WhichCells(ref3, downsample = 200))
ref3$cluster <- meta[colnames(ref3),]$cluster
table(ref3$cluster)
saveRDS(ref3, file = "reference/GBM_atlas_TAM_reference.rds")

ref1 <- readRDS("../PD-1/reference/reference_final_subset.rds")
ref2 <- readRDS("reference/GBM_atlas_myeloid_reference.rds")

#ref1: all cells
assay_for_SingleR <- GetAssayData(nonmalignant, slot="data")
cg <- as.SingleCellExperiment(ref1)
predictions <- SingleR(test=assay_for_SingleR, ref=cg, labels=cg$anno_ident)
table(predictions$labels)
celltype<-sapply(unique(as.character(nonmalignant$seurat_clusters)),
                 FUN = function(x) names(table(predictions$labels[nonmalignant$seurat_clusters==x]))[which.max(table(predictions$labels[nonmalignant$seurat_clusters==x]))])
nonmalignant$celltype <- as.vector(sapply(as.character(nonmalignant$seurat_clusters), function(x) celltype[x]))
#nonmalignant$celltype <- predictions$labels
table(nonmalignant$celltype)
col.ls<-ArchRPalettes$stallion[1:length(table(nonmalignant$celltype))]
names(col.ls)<-names(table(nonmalignant$celltype))
DimPlot(nonmalignant, group.by = "celltype", label = T, cols = col.ls)

# ref2: myeloid cells
myeloid <- nonmalignant[,nonmalignant$celltype%in%c("Macrophages","Microglial","pDCs")&nonmalignant$seurat_clusters!=4]
myeloid
myeloid <- myeloid %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 1.2) 
myeloid <- myeloid %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

assay_for_SingleR <- GetAssayData(myeloid, slot="data")
cg <- as.SingleCellExperiment(ref2)
predictions <- SingleR(test=assay_for_SingleR, ref=cg, labels=cg$cluster)
table(predictions$labels)

myeloid$immune_type <- predictions$labels
table(myeloid$immune_type,myeloid$seurat_clusters)

celltype<-sapply(unique(as.character(myeloid$seurat_clusters)),
                 FUN = function(x) names(table(predictions$labels[myeloid$seurat_clusters==x]))[which.max(table(predictions$labels[myeloid$seurat_clusters==x]))])
myeloid$immune_type <- sapply(as.character(myeloid$seurat_clusters), function(x) celltype[x])

saveRDS(myeloid, file = "result/myeloid.rds")

# ref3: TAMs
TAMs <- myeloid[,myeloid$immune_type%in%c("TAM 1","TAM 2","prol. TAM")]
TAMs <- TAMs %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 1.2) 
TAMs <- TAMs %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

assay_for_SingleR <- GetAssayData(TAMs, slot="data")
cg <- as.SingleCellExperiment(ref3)
predictions <- SingleR(test=assay_for_SingleR, ref=cg, labels=cg$cluster)
table(predictions$labels)

TAMs$TAM_type <- predictions$labels
table(TAMs$TAM_type, TAMs$seurat_clusters)

celltype<-sapply(unique(as.character(TAMs$seurat_clusters)),
                 FUN = function(x) names(table(predictions$labels[TAMs$seurat_clusters==x]))[which.max(table(predictions$labels[TAMs$seurat_clusters==x]))])
TAMs$TAM_type <- sapply(as.character(TAMs$seurat_clusters), function(x) celltype[x])

col.ls<-ArchRPalettes$stallion[1:length(table(TAMs$TAM_type))]
names(col.ls)<-names(table(TAMs$TAM_type))

saveRDS(TAMs, file = "result/TAMs.rds")

pdf("plot_new/TAMs.pdf", height = 3, width = 4)
DimPlot(TAMs, group.by = "TAM_type", cols = col.ls)
DimPlot(TAMs)
dev.off()

## 手动注释-----------
nonmalignant$celltype <- "Myeloid cell"
nonmalignant$celltype[nonmalignant$seurat_clusters==4|nonmalignant$seurat_clusters==12] <- "T cell"
nonmalignant$celltype[nonmalignant$seurat_clusters==11] <- "Endothelial cell"
nonmalignant$celltype[nonmalignant$seurat_clusters==9] <- "Fibroblast"
nonmalignant$celltype[nonmalignant$seurat_clusters==8] <- "Dendritic cell"
nonmalignant$celltype[nonmalignant$seurat_clusters==7] <- "Oligodendrocytes"
nonmalignant$celltype[nonmalignant$seurat_clusters==5] <- "Unknown"

CellColor <- ArchRPalettes$stallion[1:length(unique(nonmalignant$celltype))]
names(CellColor) <- unique(nonmalignant$celltype)

pdf("plot_new/nonmalignant_UMAP.pdf", height=4, width=5)
DimPlot(nonmalignant, group.by = "celltype", cols = CellColor)
DimPlot(nonmalignant, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(nonmalignant, group.by = "seurat_clusters", label = T)
dev.off()

saveRDS(nonmalignant, file="result/nonmalignant_withXHsample.rds")

### myeloid----------
nonmalignant<-readRDS("result/nonmalignant_withXHsample.rds")
myeloid <- nonmalignant[,nonmalignant$celltype=="Myeloid cell"]

# choice1
BMDM_MB_markers <- read.xlsx("reference/13059_2017_1362_MOESM5_ESM.xlsx",sheetIndex = 1)
BMDM_MB_markers <- lapply(1:2, function(x) na.omit(BMDM_MB_markers[,x]))
names(BMDM_MB_markers) <- c("MG","BMDM")

# choice2
BMDM_MB_markers <- list(BMDM=unique(c(BMDM_MB_markers$BMDM,c("CCR2","CD45RA","CD141","ICAM","CD1C","CD1B",
                                                             "TGFBI","FXYD5","FCGR2B","CLEC12A","CLEC10A","CD207",
                                                             "CD209","CD49D"))),
                        MG=unique(c(BMDM_MB_markers$MG,c("CX3CR1","SALL1","HEXB","P2RY12","TMEM119"))))
BMDM_MB_markers

# ## AUC
# library(AUCell)
# library(clusterProfiler)
# cells_rankings <- AUCell_buildRankings(myeloid@assays$RNA@data) 
# cells_AUC <- AUCell_calcAUC(BMDM_MB_markers, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
# myeloid$BMDM  <- as.numeric(getAUC(cells_AUC)["BMDM", ])
# myeloid$MG  <- as.numeric(getAUC(cells_AUC)["MG", ])
# FeaturePlot(myeloid, features = c("BMDM", "MG"), max.cutoff = "q90", min.cutoff = "q10", blend = T)

### GSVA
myeloid <- myeloid %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

myeloid <- myeloid %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.35)
myeloid <- myeloid %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(myeloid)
table(myeloid$seurat_clusters)

Idents(myeloid)<-"seurat_clusters"
exp=AverageExpression(myeloid,assays = "RNA")
counts2=exp[["RNA"]]
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=BMDM_MB_markers, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames = T,
                   scale = "none") #以行来标准化，这个功能很不错

subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
myeloid$TAM_type <- as.vector(sapply(myeloid$seurat_clusters,function(x) subtype[x]))
#myeloid$TAM_type <- subtype
DimPlot(myeloid, group.by = "TAM_type")
saveRDS(myeloid, file="result/myeloid_withXHsample.rds")

BMDM <- subset(myeloid, TAM_type=="BMDM")
table(BMDM$batch)

BMDM <- BMDM %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

BMDM <- BMDM %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.35)
BMDM <- BMDM %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(BMDM)
table(BMDM$seurat_clusters)

saveRDS(BMDM, file = "result/BMDM_withXHsample.rds")

## 先找younger和younger之间的差异基因
Idents(BMDM)<-"grade"
markers <- FindMarkers(BMDM,ident.1 = "older",ident.2 = "younger",test.use = "wilcox")
# markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
# markers_df <- markers_df[markers_df$p_val_adj<0.01,]
# write.table(markers_df,file="result/BMDM_markers200_grade.txt",quote=F,sep="\t",row.names=F,col.names=T)

library(ggVolcano)
data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.01)
data$row <- rownames(data)
data$regulate[data$regulate=="Up"] <- "older"
data$regulate[data$regulate=="Down"] <- "younger"
data$regulate[data$regulate=="Normal"] <- "non-significant"

# plot
pdf("plot_new/volcano_BMDM_withXHsample.pdf")
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 0.5, FDR_cut = 0.01, legend_position = "DL",
          colors = c("younger"="navy","older"="firebrick3","non-significant"="grey"), 
          fills = c("younger"="navy","older"="firebrick3","non-significant"="grey"), 
          label = "row", label_number = 10, output = FALSE)
dev.off()


# myeloid <- subset(myeloid, seurat_clusters!=6)
# myeloid <- myeloid %>%
#   FindNeighbors(reduction = "harmony") %>%
#   FindClusters(resolution=0.8)
# myeloid <- myeloid %>%
#   RunUMAP(dims = 1:40, reduction = "harmony")

# pdf("plot_new/myeloid_featureplot_younger.pdf",width = 6,height = 7)
# FeaturePlot(myeloid, features=c("CX3CR1", "TMEM119","APOC2", "TGFBI","LYZ","ANXA1"), min.cutoff = "q5", max.cutoff = "q95",cols = c("lightgrey" ,"#DE1F1F"))
# VlnPlot(myeloid, features=c("CX3CR1", "TMEM119","APOC2", "TGFBI","LYZ","ANXA1"), ncol = 2, pt.size = 0)
# #FeaturePlot(myeloid, features=c("KDM6B"), min.cutoff = "q5", max.cutoff = "q95",cols = c("lightgrey" ,"#DE1F1F"))
# #FeaturePlot(myeloid, features=c("KLF6","NR4A1","BNIP3","GPNMB"), min.cutoff = "q10", max.cutoff = "q90",cols = c("lightgrey" ,"#DE1F1F"))
# dev.off()

# pdf("plot_new/featureplot_BMDM.pdf")
# VlnPlot(BMDM, features=c("CD74","C1QB"), pt.size = 0.01, )#ggtitle("BMDM older vs younger")
# VlnPlot(BMDM, features=c("CD69","CCL3","IL1B"), pt.size = 0.01)#ggtitle("intersection between cluster 4 and older marker")
# VlnPlot(BMDM, features=c("PLIN2","NUPR1"), pt.size = 0.01)#ggtitle("intersection between cluster 4 and older marker")
# dev.off()

Idents(BMDM) <- "seurat_clusters"
markers <- FindAllMarkers(BMDM, only.pos = TRUE, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
write.table(markers_df,file="result/BMDM_markers200_withXHsample.txt",quote=F,sep="\t",row.names=F,col.names=T)

markers_df <- read.table("result/BMDM_markers200_withXHsample.txt", header = T)
markers_df = markers_df %>% 
  dplyr::filter(!grepl("^MT|^RP", gene) & avg_log2FC > 1.25) %>% 
  group_by(cluster) %>%
  arrange(cluster, desc(-log(p_val_adj))) %>%
  distinct(cluster, -log(p_val_adj), .keep_all = TRUE) %>%
  slice_head(n = 10)
#intersect(markers_df$gene[markers_df$cluster=="4"],data$row[data$regulate=="older"])

pdf("plot_new/BMDM_heatmap_withXHsample.pdf")
DoHeatmap(subset(BMDM, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
dev.off()
pdf("plot_new/BMDM_dotplot_withXHsample.pdf", height = 8, width = 5)
DotPlot(BMDM, features = unique(markers_df$gene), scale = T) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 10))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

cluster_color <- RColorBrewer::brewer.pal(10,"Set2")
names(cluster_color) <- unique(BMDM$seurat_clusters)
DimPlot(BMDM, label=T, cols = cluster_color)
FeaturePlot(BMDM, "NUPR1", max.cutoff = "q90", min.cutoff = "q10")

# myeloid$TAM_type<-NA
# myeloid$TAM_type[myeloid$seurat_clusters==2|myeloid$seurat_clusters==1] <- "MHC"
# myeloid$TAM_type[myeloid$seurat_clusters==0|myeloid$seurat_clusters==3] <- "Phago/lipid"
# myeloid$TAM_type[myeloid$seurat_clusters==4|myeloid$seurat_clusters==5] <- "Hypoxia"
# 
# ### GSVA
# gene_sets <- lapply(2:14, function(x) {
#   markers <- readxl::read_xlsx("reference/Human_ND_GBM_TAM.xlsx", sheet = x)
#   gene <- markers$gene[markers$p_val_adj<0.05&markers$avg_logFC>1]
#   c(unique(markers$cluster.1),gene)
# })
# names(gene_sets) <- sapply(gene_sets, function(x) x[1])
# gene_sets <- lapply(gene_sets, function(x) x[-1])
# 
# obj_list <- SplitObject(myeloid,"sample")
# 
# obj_list <- lapply(obj_list, function(myeloid){
#   myeloid <- myeloid %>%
#     FindNeighbors() %>%
#     FindClusters(resolution=0.8)
#   
#   Idents(myeloid)<-"seurat_clusters"
#   exp=AverageExpression(myeloid,assays = "RNA") 
#   counts2=exp[["RNA"]]
#   
#   GSVA_hall <- gsva(expr=as.matrix(counts2), 
#                     gset.idx.list=BMDM_MB_markers, 
#                     mx.diff=T, # 数据为正态分布则T，双峰则F
#                     kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
#                     parallel.sz=1) # 并行线程数目
#   head(GSVA_hall)
#   
#   # pheatmap::pheatmap(GSVA_hall, #热图的数据
#   #                    cluster_rows = T,#行聚类
#   #                    cluster_cols =T,#列聚类，可以看出样本之间的区分度
#   #                    show_colnames=T,
#   #                    scale = "column") #以行来标准化，这个功能很不错
#   
#   subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
#   myeloid$TAM_type <- sapply(myeloid$seurat_clusters,function(x) subtype[x])
#   myeloid
# })
# 
# myeloid <- merge(obj_list[[1]],c(obj_list[[2]],obj_list[[3]],obj_list[[4]],obj_list[[5]],obj_list[[6]],obj_list[[7]],obj_list[[8]]))
# myeloid <- NormalizeData(myeloid)
# myeloid <- FindVariableFeatures(myeloid)
# myeloid <- ScaleData(myeloid)
# myeloid <- RunPCA(myeloid)
# myeloid <- myeloid %>%
#   RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
# myeloid <- myeloid %>%
#   FindNeighbors(reduction = "harmony") %>%
#   FindClusters(resolution=0.5)
# myeloid <- myeloid %>%
#   RunUMAP(dims = 1:40, reduction = "harmony")

TAM_color<-ArchRPalettes$stallion[1:length(table(BMDM$TAM_type))]
names(TAM_color)<-names(table(BMDM$TAM_type))
pdf("plot_new/BMDM_basics_withXHsample.pdf", width = 4, height = 3)
DimPlot(BMDM, label=T, cols = cluster_color)
DimPlot(BMDM, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(BMDM, group.by = "grade", cols = GradeColor)
DimPlot(BMDM, group.by = "TAM_type", cols = TAM_color)
dev.off()

table(myeloid$sample)
library(forcats)
library(ggplot2)
library(gridExtra)

col.ls<-ArchRPalettes$stallion[1:length(table(myeloid$TAM_type))]
names(col.ls)<-names(table(myeloid$TAM_type))

CellInfo <- myeloid@meta.data
P1=CellInfo %>% ggplot(aes(x=sample, fill=TAM_type)) +scale_fill_manual(values = col.ls)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(TAM_type , fill=TAM_type))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  scale_fill_manual(values = col.ls)+
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
pdf("plot_new/proportion_TAM_withXHsample.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
grid.arrange(grobs=list(P3g,P4g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

## Mo-TAM----------------
Mo <- subset(myeloid, TAM_type=="Mo_TAM")
Mo <- Mo %>%
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
Mo <- Mo %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.6)
Mo <- Mo %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

pdf("plot_new/Mo_featureplot.pdf",width = 6,height = 8)
FeaturePlot(Mo, features=c("CCL3","CCL4","IL1B","TNF","GPNMB","LGALS3","ANXA2","VIM"), ncol = 2, min.cutoff = "q10", max.cutoff = "q90",cols = c("lightgrey" ,"#DE1F1F"))
VlnPlot(Mo, features=c("CCL3","CCL4","IL1B","TNF","GPNMB","LGALS3","ANXA2","VIM"), ncol = 2, pt.size = 0)
dev.off()

Idents(Mo)<-"grade"
markers <- FindAllMarkers(Mo, test.use = "wilcox")
library(scRNAtoolVis)
pdf("plot_new/Mo_volcano.pdf", height = 7, width = 5)
jjVolcano(diffData = markers,
          log2FC.cutoff = 0.5, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          #aesCol = c('purple','orange'), #设置点的颜色
          tile.col = cluster_color, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因 
)
dev.off()

markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.01,]
write.table(markers_df,file="result/Mo_markers200.txt",quote=F,sep="\t",row.names=F,col.names=T)

#markers_df <- read.table("result/myeloid_markers200.txt", header = T)
#markers_df <- markers_df[-grep("^MT",markers_df$gene),]
markers_df = markers_df %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("plot_new/Mo_heatmap.pdf")
DoHeatmap(subset(Mo, downsample = 500), features = markers_df$gene, size = 3)+ 
  scale_fill_viridis() + theme(text = element_text(size = 8)) + NoLegend()
dev.off()
pdf("plot_new/Mo_dotplot.pdf", height = 10, width = 6)
DotPlot(Mo, features = unique(markers_df$gene)) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 8))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ #legend
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

cluster_color <- RColorBrewer::brewer.pal(10,"Set2")
names(cluster_color) <- unique(myeloid$seurat_clusters)

Mo$subtype <- "regulatory"
Mo$subtype[Mo$seurat_clusters==2] <- "lipid-associated"
Mo$subtype[Mo$seurat_clusters==0] <- "transitional"
subtype_color <- RColorBrewer::brewer.pal(7,"Set2")[3:5]
names(subtype_color) <- c(unique(Mo$subtype))

pdf("plot_new/Mo_basics.pdf", width = 4, height = 3)
DimPlot(Mo, label=T, cols = cluster_color)
DimPlot(Mo, group.by = "sample", cols = SampleColor, shuffle = T)
DimPlot(Mo, group.by = "grade", cols = GradeColor)
DimPlot(Mo, group.by = "subtype", cols = subtype_color)
dev.off()

myeloid<-readRDS("result/myeloid.rds")
M1<-c("CCL2","CCL3","CCL4","CCL5","CCL8","CCR7","CD74","CSF2",
      "CXCL10","HLA-DRA","HLA-DRB","IFNG","IL1B","IL1R1","IL6","INOS",
      "IRF5","NFKB1","TLR2","TLR4","TNF")
M2<-c("ARG1","CD74","CCL1","CCL17","CCL22","CXCL16","CXCR4","HLA-DRA",
      "HLA-DRB","IL10","IL4","IRF4","MRC1","NFKB1","TGFB1","TNF")






