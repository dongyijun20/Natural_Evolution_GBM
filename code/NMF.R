library(remotes)
remotes::install_github("carmonalab/GeneNMF") # R>4.3,换成本地
library(GeneNMF) # 用于NMF
library(Seurat)
library(msigdbr)
library(fgsea)

# seurat对象的读取，以rds文件为例
data.path = "../test.rds"
seu <- readRDS(data.path)
# 提取目标细胞，以肿瘤细胞为例
seu = subset(seu, 
             subset = cell_type == "tumor cell")
# 按照样本分开
seu.list <- SplitObject(seu, split.by = "orig.ident")

geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", 
                             k=4:9, L1=c(0,0),
                             do_centering=TRUE, 
                             nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=10,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.3)

ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))

