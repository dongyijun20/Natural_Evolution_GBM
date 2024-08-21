library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(infercnv)

setwd("~/neurosurgery/NES")
pbmc <- readRDS("result/merged_withXHsample.rds")
subset <- subset(pbmc, sample == "XH01_2")

counts_matrix <- GetAssayData(object = subset, slot = 'counts')
annot <- as.data.frame(subset$tumor_ident)

# Load gene order file
gene_order <- read.table('geneLocate.txt', header = F,row.names = 1)
# Create infercnv object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = annot, 
                                     gene_order_file = gene_order, 
                                     ref_group_names = "non-malignant")

infercnv_obj <- infercnv::run(infercnv_obj, 
                              cutoff = 1,
                              out_dir = paste0('infercnv/XH01_2_malignant'), 
                              cluster_by_groups = F, 
			                        resume = T,
                              analysis_mode="subclusters", 
			                        tumor_subcluster_partition_method = "random_trees",
                              denoise = T, 
                              HMM = T, 
                              output_format = NA,
                              num_threads = 10)

# infercnv_obj <- readRDS("infercnv/older_malignant/run.final.infercnv_obj")
# infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
#                    plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
#                    output_filename = paste0('infercnv/older_malignant/infercnv.pdf'),
#                    output_format = "pdf") #保存为pdf文件


