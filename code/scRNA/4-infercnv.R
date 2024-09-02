library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(infercnv)
pbmc <- readRDS("result/merged_withXHsample.rds")
pbmc <- FindClusters(pbmc, resolution = 2)

table(pbmc$celltype_bped_main)
table(pbmc$celltype_hpca_main)
table(pbmc$celltype_iced_main)
table(pbmc$celltype_mid_main)
immune_bped<-colnames(pbmc)[pbmc$celltype_bped_main%in%c("B-cells","CD4+ T-cells","CD8+ T-cells","DC","Macrophages","Monocytes","NK cells")]
immune_hpca<-colnames(pbmc)[pbmc$celltype_hpca_main%in%c("B_cell","T_cells","DC","Macrophage","Monocyte","NK_cell")]
immune_mid<-colnames(pbmc)[pbmc$celltype_mid_main%in%c("B cells","CD4+ T cells","CD8+ T cells","Dendritic cells","T cells","Monocytes","NK cells")]
immune <- names(table(pbmc$seurat_clusters))[table(pbmc$seurat_clusters[intersect(intersect(immune_bped,immune_hpca),
                                                                       immune_mid)])/table(pbmc$seurat_clusters)>0.9]
# table(pbmc$seurat_clusters)
# 
# subset <- pbmc[,sample(colnames(pbmc),10000)]
# saveRDS(subset, file = "result/merged_withXHsample_subset.rds")

subset <- readRDS("result/merged_withXHsample_subset.rds")
subset <- JoinLayers(subset)
counts_matrix <- LayerData(object = subset, layer = "counts")
annot <- as.data.frame(subset$seurat_clusters)

# Load gene order file
gene_order <- read.table('geneLocate.txt', header = F,row.names = 1)
# Create infercnv object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = annot, 
                                     gene_order_file = gene_order, 
                                     ref_group_names = immune)

# Run inferCNV
infercnv_obj <- infercnv::run(infercnv_obj, 
                              cutoff = 0.1,
                              out_dir = paste0('infercnv/merged_withXHsample_subset'), 
                              cluster_by_groups = T,
                              resume_mode = T,
                              denoise = T, 
                              HMM = T, 
                              output_format = 'pdf',
                              num_threads = 10)

infercnv_obj<-readRDS("infercnv/merged_withXHsample_subset/run.final.infercnv_obj")
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = paste0('infercnv/merged_new_subset/infercnv_re.pdf'),
                   output_format = "pdf") #保存为pdf文件

# Identify malignant cells
subset <- pbmc[,colnames(infercnv_obj@expr.data)]
seu <- add_to_seurat(subset, infercnv_output_path = 'infercnv/merged_withXHsample_subset/')
cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$cnv_avg <- rowMeans(cnvs)
table(seu$cnv_avg,seu$seurat_clusters)

# For the majority of samples cut-off > 0.1 (exceptions below)
sapply(sort(unique(seu$seurat_clusters)), function(x) summary(seu$cnv_avg[seu$seurat_clusters==x]))
seu$malignant <- ifelse(seu$cnv_avg > 0.1, 'malignant', 'non-malignant')
table(seu$malignant,seu$seurat_clusters)

# Add CNV metrics
cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$proportion_scaled_cnv_avg <- rowMeans(cnvs)

cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$proportion_cnv_avg <- rowMeans(cnvs)

cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$has_cnv_avg <- rowMeans(cnvs)

table(seu$malignant)
table(seu$malignant,seu$seurat_clusters)
malignant_cluster <- colnames(table(seu$malignant,seu$seurat_clusters))[table(seu$malignant,seu$seurat_clusters)[1,]>table(seu$malignant,seu$seurat_clusters)[2,]]
malignant_cluster

pbmc$tumor_ident <- "non-malignant"
pbmc$tumor_ident[pbmc$seurat_clusters%in%malignant_cluster]<-"malignant"
table(pbmc$tumor_ident)

saveRDS(pbmc, file = "result/merged_withXHsample.rds")

pdf("plot_new/infercnv_withXHsample.pdf")
DimPlot(pbmc, group.by = "tumor_ident",cols = c("malignant"="red","non-malignant"="blue"))+
  ggtitle("Malignant Cells Identification",subtitle = paste("Fraction = ", 100*round(sum(pbmc$tumor_ident=="malignant")/ncol(pbmc),2),"%",sep = ""))
dev.off()

dat <- read.table("infercnv/merged_withXHsample_subset/expr.infercnv.dat")
cell.group <- read
ht = Heatmap(dat,show_row_names=F,show_column_names=F,cluster_columns=F,
             row_dend_width = unit(0.1,"npc"),
             clustering_method_rows = hclust.method,
             clustering_distance_rows = distant.method,
             split=pbmc$seurat_clusters,column_split=gene.chrom$chrom,column_gap=unit(0.1,"mm"),
             show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "",
                                                                     legend_direction="vertical",nrow=1, legend_height=unit(0.4,"npc")),
             left_annotation = cluster.bar,
             top_annotation = chrom.bar)

draw(ht,heatmap_legend_side = "right", annotation_legend_list=list(lgd), annotation_legend_side="bottom")


# setwd("CopyKat/")
# copykat.test <- copykat(rawmat=exp.rawdata, 
#                         id.type="S", 
#                         cell.line="no", 
#                         ngene.chr=5, 
#                         win.size=25, 
#                         KS.cut=0.15, 
#                         sam.name="integrated", 
#                         distance="euclidean", 
#                         plot.genes = F,
#                         norm_cells = immune,
#                         n.cores = 1)
# setwd("../")
# pred.test <- data.frame(copykat.test$prediction)
# #pred.test <- read.csv("CopyKat/RD002_copykat_prediction.txt",sep="\t")
# CNA.test <- data.frame(copykat.test$CNAmat)
# #CNA.test <- read.csv("CopyKat/RD002_copykat_CNA_results.txt",sep="\t")

# tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
# tumor.cells <- gsub("-",".",tumor.cells)
# tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
# hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads = 4, method = "euclidean"), method = "ward.D2")
# #hcc <- readRDS("CopyKat/RD002_copykat_clustering_results.rds")
# hc.umap <- cutree(hcc,2)
# names(hc.umap)<-gsub("\\.","-",names(hc.umap))
# 
# pbmc@meta.data$copykat.pred <- pred.test$copykat.pred
# pbmc@meta.data$copykat.tumor.pred <- rep("NA", nrow(pbmc@meta.data))
# pbmc@meta.data$copykat.tumor.pred[rownames(pbmc@meta.data) %in% names(hc.umap[hc.umap==1])] <- "Malignant"
# pbmc@meta.data$copykat.tumor.pred[rownames(pbmc@meta.data) %in% names(hc.umap[hc.umap==2])] <- "Nonmalignant"
# # table(pbmc@meta.data$copykat.tumor.pred)
# pbmc<-subset(pbmc,copykat.tumor.pred!="NA")
# levels(pbmc$copykat.tumor.pred)<-c("Malignant","Nonmalignant")
# tumor.fraction <- sapply(unique(pbmc$seurat_clusters), function(x) {
#   sum(pbmc$copykat.tumor.pred[pbmc$seurat_clusters==x]=="Malignant")/sum(pbmc$seurat_clusters==x)
# })
# names(tumor.fraction) <- unique(pbmc$seurat_clusters)
# tumor.cluster <- names(tumor.fraction)[tumor.fraction > 0.8]
# pbmc$tumor.ident <- sapply(pbmc$seurat_clusters, function(x) {
#   if(x%in%tumor.cluster) return("Malignant") else return("Nonmalignant")
# })
# 
# 
# # VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
# # pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 60000 & percent.mt < 60)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# pbmc <- ScaleData(pbmc, features = rownames(pbmc))
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- FindNeighbors(pbmc, dims = 1:30)
# pbmc <- FindClusters(pbmc, resolution = 0.6)
# pbmc <- RunUMAP(pbmc, dims = 1:30)
# pbmc <- RunTSNE(pbmc, dims = 1:30)
# saveRDS(pbmc,file = "RD002.rds")
# 
# pdf("CopyKat/plots_RD002.pdf")
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
# DimPlot(pbmc, label = T)
# DimPlot(pbmc, group.by = "copykat.pred")
# DimPlot(pbmc, group.by = "copykat.tumor.pred")
# DimPlot(pbmc, group.by = "tumor.ident", cols = c("Malignant"="red","Nonmalignant"="blue"))+
#   ggtitle("Malignant Cells Identification",subtitle = paste("Fraction = ", 100*round(sum(pbmc$tumor.ident=="Malignant")/ncol(pbmc),2),"%",sep = ""))
# dev.off()

nonmalignant <- subset(pbmc, tumor_ident=="non-malignant")
nonmalignant <- FindVariableFeatures(nonmalignant, selection.method = "vst", nfeatures = 2000)
nonmalignant <- ScaleData(nonmalignant, features = rownames(nonmalignant))
nonmalignant <- RunPCA(nonmalignant, features = VariableFeatures(object = nonmalignant))
nonmalignant <- FindNeighbors(nonmalignant, dims = 1:40)
nonmalignant <- FindClusters(nonmalignant, resolution = 1)
nonmalignant <- RunUMAP(nonmalignant, dims = 1:40)
saveRDS(nonmalignant,file = "result/nonmalignant_withXHsample.rds")

malignant <- subset(pbmc, tumor_ident=="malignant")
malignant <- FindVariableFeatures(malignant, selection.method = "vst", nfeatures = 2000)
malignant <- ScaleData(malignant, features = rownames(malignant))
malignant <- RunPCA(malignant, features = VariableFeatures(object = malignant))
malignant <- FindNeighbors(malignant, dims = 1:40)
malignant <- FindClusters(malignant, resolution = 1)
malignant <- RunUMAP(malignant, dims = 1:40)
saveRDS(malignant,file = "result/malignant_withXHsample.rds")

pdf("plot_new/divided_withXHsample.pdf")
DimPlot(nonmalignant, label=T)
DimPlot(malignant, label=T)
DimPlot(malignant, group.by = "sample", cols = SampleColor)
DimPlot(malignant, group.by = "grade", cols = GradeColor)
dev.off()

## after determine which is older and which is younger--------

malignant <- readRDS("result/malignant_withXHsample.rds")
younger <- subset(malignant, grade=="younger")
older <- subset(malignant, grade=="older")

younger <- NormalizeData(younger)
younger <- FindVariableFeatures(younger)
younger <- ScaleData(younger)
younger <- RunPCA(younger, npcs = 50)
younger <- RunUMAP(younger, dims = 1:40)
younger <- FindNeighbors(younger, dims = 1:40)
younger <- FindClusters(younger, resolution = 0.5)

pdf("plot_new/younger_malignant_withXHsample.pdf")
DimPlot(younger, label = T)
DimPlot(younger, group.by = "sample", cols = SampleColor)
dev.off()

older <- NormalizeData(older)
older <- FindVariableFeatures(older)
older <- ScaleData(older)
older <- RunPCA(older, npcs = 50)
older <- RunUMAP(older, dims = 1:40)
older <- FindNeighbors(older, dims = 1:40)
older <- FindClusters(older, resolution = 0.5)

pdf("plot_new/older_malignant_withXHsample.pdf")
DimPlot(older, label = T)
DimPlot(older, group.by = "sample", cols = SampleColor)
dev.off()

saveRDS(younger, file = "result/younger_malignant_withXHsample.rds")
saveRDS(older, file = "result/older_malignant_withXHsample.rds")
