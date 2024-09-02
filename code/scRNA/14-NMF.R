library(remotes)
remotes::install_github("carmonalab/GeneNMF") #from Github
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)

# seurat对象的读取，以rds文件为例
data.path = "result/myeloid_withXHsample.rds"
seu <- readRDS(data.path)
# # 提取目标细胞，以肿瘤细胞为例
# seu = subset(seu, 
#              subset = cell_type == "tumor cell")
# 按照样本分开
seu.list <- SplitObject(seu, split.by = "orig.ident")

geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", 
                             k=4:9, L1=c(0,0), 
                             center = T, 
                             scale = T,
                             nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, 
                                        nMP=7,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.3)

## check each metaprogram
library(msigdbr)
library(fgsea)

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "H")
})
top_p

mp.genes <- geneNMF.metaprograms$metaprograms.genes

seu <- AddModuleScore_UCell(seu, features = mp.genes, assay="RNA", ncores=4, name = "")

pdf("plot_new/NMF_myeloid.pdf")
plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
VlnPlot(seu, features=names(mp.genes), group.by = "seurat_clusters",
        pt.size = 0, ncol=4)
VlnPlot(seu, features=names(mp.genes), group.by = "TAM_type",
        pt.size = 0, ncol=4)
FeaturePlot(seu, features = c("MP7","NUPR1"), blend = T)
dev.off()

## signiture scores--------
matrix <- seu@meta.data[,names(mp.genes)]

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                        cell.embeddings = dimred,
                                        assay.used = "RNA",
                                        key = "MP_",
                                        global = FALSE)

set.seed(123)
seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")

FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

DimPlot(seu, reduction = "umap_MP", group.by = "seurat_clusters", label=T) + theme(aspect.ratio = 1,
                                                                                    axis.text = element_blank(),
                                                                                    axis.title = element_blank(),
                                                                                    axis.ticks = element_blank()) + ggtitle("Original cell types") + NoLegend()

