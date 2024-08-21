library(Seurat)
library(patchwork)
library(monocle)

malignant <- readRDS("result/malignant.rds")
immune.combined <- subset(malignant, sample%in%c("TT02_1","TT02_2"))

## monocle2
DefaultAssay(immune.combined) <- 'RNA'
cds <- as.CellDataSet(immune.combined)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
## ordering by marker gene per cluster
# Idents(immune.combined) <- "celltype_bped_main"
# deg <- FindAllMarkers(immune.combined)
# write.csv(deg, file = "markers.csv", quote = F)

# Idents(immune.combined) <- "sample"
# deg <- FindAllMarkers(immune.combined)
# write.csv(deg, file = "markers_sample.csv", quote = F)
# deg <- read.table("result/myeloid_markers200_grade.txt", header = T)
# sel.gene <- unique(deg$gene[deg$p_val_adj<0.05])

# cds <- detectGenes(cds, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(cds),
#                                     num_cells_expressed >= 10))
# marker_diff <- markerDiffTable(cds[expressed_genes,],
#                                cth,
#                                residualModelFormulaStr = "~Media + num_genes_expressed",
#                                cores = 1)
# candidate_clustering_genes <-
#   row.names(subset(marker_diff, qval < 0.01))
# marker_spec <-
#   calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
# head(selectTopMarkers(marker_spec, 3))
# semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
# HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
# plot_ordering_genes(HSMM)


Idents(immune.combined) <- "grade"
deg <- FindAllMarkers(immune.combined)
dim(deg)
deg = deg %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
sel.gene <- deg$gene
cds <- monocle::setOrderingFilter(cds, sel.gene)

## dimension reduciton
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

## ordering cells
cds <- monocle::orderCells(cds)

saveRDS(cds, file = "result/TT02_cds.rds")

SampleColor <- readRDS("color/SampleColor.rds")
pdf("plot_new/TT02_monocle2.pdf")
monocle::plot_cell_trajectory(cds, color_by = "sample")+ scale_color_manual(values=SampleColor)
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()


#devtools::install_github("rcannood/SCORPIUS")
library(SCORPIUS)

####至此，github官网上expression用的是data的数据，而不是scale data，试试看有什么区别
expression_data <- t(as.matrix(immune.combined@assays$RNA@data))
expression_scale <- t(as.matrix(immune.combined@assays$RNA@scale.data))

group_name =  as.factor(as.character(immune.combined$grade))
table(group_name)
dim(expression_data)
dim(expression_scale)

#######expression_data----
space <- reduce_dimensionality(expression_data, "spearman")
draw_trajectory_plot(space, group_name, contour = TRUE,)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, group_name, traj$path, contour = TRUE)
# warning: setting num_permutations to 10 requires a long time (~30min) to run!
# set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
gimp <- gene_importances(
  expression_data,
  traj$time,
  num_permutations = 0,
  num_threads = 8,
  ntree = 10000,
  ntree_perm = 1000
)
#gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[1:200]
expr_sel <- scale_quantile(expression_data[,gene_sel])

# Draw a time series heatmap
time <- traj$time
draw_trajectory_heatmap(expr_sel, time)

## Also show the progression groupings
draw_trajectory_heatmap(expr_sel, time,
                        progression_group=group_name)

# ######expression_scale----
# set.seed(123)
# space2 <- reduce_dimensionality(expression_scale, "spearman")
# draw_trajectory_plot(space2, group_name, contour = TRUE,)
# traj2 <- infer_trajectory(space2)
# draw_trajectory_plot(space2, group_name, traj2$path, contour = TRUE)
# # warning: setting num_permutations to 10 requires a long time (~30min) to run!
# # set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
# gimp2 <- gene_importances(
#   expression_scale, 
#   traj2$time, 
#   num_permutations = 0, 
#   num_threads = 8, 
#   ntree = 10000,
#   ntree_perm = 1000
# ) 
# #gimp2$qvalue <- p.adjust(gimp2$pvalue, "BH", length(gimp2$pvalue))
# gene_sel2 <- gimp2$gene[1:200]
# expr_sel2 <- scale_quantile(expression_scale[,gene_sel2])
# 
# # Draw a time series heatmap
# time2 <- traj2$time
# # draw_trajectory_heatmap(expr_sel2, time2)
# # 
# # ## Also show the progression groupings
# # draw_trajectory_heatmap(expr_sel2, time2, 
# #                         progression_group=group_name)
# # draw_trajectory_heatmap(expr_sel2, time2, 
# #                         progression_group=group_name,
# #                         show_labels_row=T)
# modules_2 <- extract_modules(scale_quantile(expr_sel2), traj2$time, verbose = F)

pdf("plot_new/SCORPIUS_immune.combined.pdf")
draw_trajectory_plot(space2, group_name, traj2$path, contour = TRUE)
draw_trajectory_heatmap(expr_sel2, time2, 
                        progression_group=group_name,
                        show_labels_row=T, fontsize_row = 4,
                        modules=modules_2)
dev.off()

### monocle3----------------
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(viridis)

immune.combined <- readRDS("BMDM.rds")
myeloid.cds <- as.cell_data_set(immune.combined)
myeloid.cds<-preprocess_cds(myeloid.cds)

# myeloid.cds <- align_cds(myeloid.cds, alignment_group = "batch")
# myeloid.cds <- reduce_dimension(myeloid.cds)
# plot_cells(myeloid.cds, label_groups_by_cluster=FALSE,  color_cells_by = "mac_state")
# plot_cells(myeloid.cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters")
# 
myeloid.cds <- cluster_cells(cds = myeloid.cds, reduction_method = "UMAP")
# myeloid.cds <- learn_graph(myeloid.cds, use_partition = TRUE)
# plot_cells(myeloid.cds)

myeloid.cds <- learn_graph(myeloid.cds, use_partition = F, close_loop = FALSE,
                           learn_graph_control = list(minimal_branch_len=15,euclidean_distance_ratio=5))
myeloid.cds <- order_cells(myeloid.cds, reduction_method = "UMAP")
plot_cells(myeloid.cds)
#hsc <- colnames(immune.combined)[which(immune.combined$seurat_clusters=='0')]
# order cells
#myeloid.cds <- order_cells(myeloid.cds, reduction_method = "UMAP", root_cells = hsc)
#plot_cells(myeloid.cds)

saveRDS(myeloid.cds,file="result/monocle3.rds")

myeloid.cds <- readRDS("result/monocle3.rds")
immune.combined <- AddMetaData(
  object = immune.combined,
  metadata = myeloid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "monocle3"
)

pdf(file="monocle3.pdf", height = 3, width = 3)
# plot trajectories colored by pseudotime
plot_cells(
  cds = myeloid.cds,
  color_cells_by = "mac_state",
  show_trajectory_graph = TRUE,
  label_leaves = TRUE, 
  label_branch_points = TRUE)
plot_cells(
  cds = myeloid.cds,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_leaves = TRUE, 
  label_branch_points = TRUE)
plot_cells(myeloid.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
FeaturePlot(immune.combined, c("monocle3"), pt.size = 0.1) & scale_color_viridis()
dev.off()

S2_BMDM=c("HSPA5","G0S2","IL8","NUPR1","VIM","LDHA","CSTB","PLIN2","S100A10")
ciliated_cds_pr_test_res <- graph_test(myeloid.cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.01))
length(pr_deg_ids)
intersect(pr_deg_ids,S2_BMDM)
rowData(myeloid.cds)$gene_name <- rownames(myeloid.cds)
rowData(myeloid.cds)$gene_short_name <- rowData(myeloid.cds)$gene_name

plot_cells(myeloid.cds, genes=S2_BMDM, 
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

pdf("NUPR1.pdf",width = 4, height = 3.5)
FeaturePlot(immune.combined,"NUPR1", min.cutoff = "q10", max.cutoff = "q90")+scale_color_viridis()
dev.off()

cds_subset <- myeloid.cds[,colData(myeloid.cds)$seurat_clusters%in%c("0","2")&myeloid.cds@clusters$UMAP$clusters==1]
plot_cells(cds_subset)
plot_cells(cds_subset, genes=c("NUPR1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
intersect(pr_deg_ids,S2_BMDM)
saveRDS(pr_deg_ids,file="result/pr_deg_ids.rds")

library(Matrix)
cds_subset<-preprocess_cds(cds_subset, num_dim = 30,norm_method = "none")
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), 
                                cell_group=partitions(cds)[colnames(cds_subset)])
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

# cds_subset <- myeloid.cds[rowData(myeloid.cds)$gene_short_name %in% S2_BMDM,]
# 
# gene_fits <- fit_models(cds_subset, model_formula_str = "~mac_state + batch")
# fit_coefs <- coefficient_table(gene_fits)
# fit_coefs

AFD_lineage_cds <- cds_subset[rowData(cds_subset)$gene_short_name %in% S2_BMDM]
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 2,color_cells_by="seurat_clusters",min_expr=0.1)


