## Figure 3A: bulk RNA-seq location inference-------------
BMDM <- readRDS("result/BMDM_withXHsample_harmony.rds")

older_markers <- FindMarkers(BMDM, ident.1 = "older", ident.2 = "younger", group.by = "grade", logfc.threshold = 0.3, only.pos = T)
markers_list <- readRDS("result/BMDM_markers_filter.rds")
mp.genes <- readRDS("result/mp.genes.rds")
ferroptosis <- readRDS("reference/ferroptosis_markers.rds")
markers<-as.list(read_excel("reference/mmc2.xlsx", skip = 3))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1/S","G2/M")]
BMDM_MG_markers <- readRDS("reference/BMDM_MG_markers.rds")

gene_list <- list(MES = markers$MES,
                  BMDM = BMDM_MG_markers$BMDM,
                  S2 = markers_list$S2, 
                  MP1 = mp.genes$MP1, 
                  Hypoxia = hallmarks$gene_symbol[hallmarks$gs_name=="HALLMARK_HYPOXIA"],
                  Ferroptosis_neg = ferroptosis$negative,
                  Ferroptosis_pos = ferroptosis$positive,
                  Ferroptosis_dual = ferroptosis$dual)

saveRDS(gene_list, file = "reference/gene_list_Fig3A.rds")

## Figure 3A: spatial transcriptomic analysis----------------
library(paletteer)
library(Seurat)
library(msigdbr)
library(viridis)

# samples from our center
GBM_2P <- readRDS("result/46561292-paratumor.rds")
GBM_2T <- readRDS("result/46561292-tumor_core.rds")

# other collected samples
GBM1 <- Load10X_Spatial("collected_data/GBM1/")
GBM1 <- SCTransform(GBM1, assay = "Spatial", verbose = FALSE)
GBM2 <- Load10X_Spatial("collected_data/GBM2/")
GBM2 <- SCTransform(GBM2, assay = "Spatial", verbose = FALSE)
GBM3 <- Load10X_Spatial("collected_data/GBM3/")
GBM3 <- SCTransform(GBM3, assay = "Spatial", verbose = FALSE)
GBM4 <- Load10X_Spatial("collected_data/GBM4/")
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM5_1 <- Load10X_Spatial("collected_data/GBM5-1/")
GBM5_1 <- SCTransform(GBM5_1, assay = "Spatial", verbose = FALSE)
GBM5_2 <- Load10X_Spatial("collected_data/GBM5-2/")
GBM5_2 <- SCTransform(GBM5_2, assay = "Spatial", verbose = FALSE)

GBM_spatial <- readRDS("reference/GBM.spatial.integrated.rds")

obj_list <- list(GBM_2T,GBM_2P,GBM1,GBM2,GBM3,GBM4,GBM5_1,GBM5_2)

obj_list <- lapply(obj_list, function(x) AddModuleScore_UCell(x, gene_list, name = NULL))

saveRDS(obj_list, file = "result/obj_list_ST.rds")

plot_spatial_feature <- function(GBM, features){
  plots <- sapply(features, function(x){
    SpatialFeaturePlot(GBM, x, stroke = 0,
                       pt.size.factor = 2, alpha = c(0.2,0.8),
                       max.cutoff = "q95", min.cutoff = "q5") + 
      scale_fill_gradientn(colors = paletteer_c("viridis::inferno", n = 250))
  })
  return(wrap_plots(plots, nrow = 1))
}

plots <- lapply(obj_list, function(x) plot_spatial_feature(x, names(gene_list)))

pdf("figures/Fig3A.pdf", width = 15, height = 3)
plots
dev.off()

## Figure 3B: statistics of spatial transcriptomics-------------
library(pheatmap)

scores_df <- lapply(obj_list, function(x) x@meta.data[, names(gene_list)]) %>% 
  bind_rows()

correlation_matrix <- cor(scores_df, method = "pearson")

pdf("figures/Fig3B.pdf", height = 5, width = 5)
pheatmap(correlation_matrix, 
         display_numbers = TRUE, 
         clustering_method = "complete", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Correlation of Signatures in ST (n = 8)")
dev.off()

## Figure 3C: scheme figure of serial stereotactic biopsy
## Figure 3D: volcano plot of serial stereotactic biopsy data------------
library(ggVolcano)
library(harmony)
library(readxl)
library(GSVA)

merge_serial <- readRDS("result/merge_WJS.rds")
TAM <- readRDS("result/TAM_WJS_new.rds")
malignant <- readRDS("result/malignant_WJS.rds")

PositionColor <- c("SX"="steelblue","HS"="tomato")

pdf("figures/Sup_Fig3A.pdf", height = 4, width = 5)
DimPlot(TAM, group.by = "position", cols = PositionColor, shuffle = T)
DimPlot(TAM, group.by = "TAM_type")
FeaturePlot(TAM,c("NUPR1"))
dev.off()

# pie chart
prop_data <- TAM@meta.data %>%
  count(position, TAM_subtype) %>%
  group_by(position) %>%
  mutate(proportion = n / sum(n),
         label = paste0(TAM_subtype, " (", round(proportion * 100), "%)"))

# Plot pie charts using facet
pdf("figures/Fig3C.pdf", height = 4, width = 6)
ggplot(prop_data, aes(x = "", y = proportion, fill = TAM_subtype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~position) +
  theme_void() +
  labs(title = "TAM subtype composition in each position") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
  theme(strip.text = element_text(size = 14)) +
  scale_fill_manual(values = c("S2 BMDM"='#ea9994', "other TAMs"='#9cd2ed'))
dev.off()

# define subtypes of malignant cells
markers<-as.list(read_excel("reference/mmc2.xlsx", skip = 3))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1/S","G2/M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)

exp=AverageExpression(malignant, assays = "RNA", group.by = "seurat_clusters")
counts2=exp[["RNA"]]
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=markers, 
                  mx.diff=T, 
                  kcdf="Gaussian", 
                  parallel.sz=4) 
head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall, 
                   cluster_rows = T,
                   cluster_cols =T,
                   show_colnames=T,
                   scale = "none") 

subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
names(subtype) <- gsub("g", "", names(subtype))
malignant$subtype <- as.vector(sapply(malignant$seurat_clusters,function(x) subtype[x]))
table(malignant$subtype)

# define subtype colors
SubtypeColor <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3')
names(SubtypeColor) <- unique(malignant$subtype)

pdf("figures/Sup_Fig3B.pdf", height = 4, width = 5.5)
DimPlot(malignant, group.by = "subtype", cols = SubtypeColor)
dev.off()

# plot the proportion of subtypes
library(forcats)
library(ggplot2)
library(gridExtra)
CellInfo <- malignant@meta.data
P1=CellInfo %>% ggplot(aes(x=position, fill=fct_rev(subtype))) +scale_fill_manual(values = SubtypeColor)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(subtype , fill=subtype))+geom_bar(stat="count",colour = "black",width = 0.7)+  scale_fill_manual(values = SubtypeColor)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)

pdf("figures/Sup_Fig3C.pdf",width=5,height=4)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

Idents(TAM)<-"position"
markers <- FindMarkers(TAM, ident.1 = "HS",ident.2 = "SX",test.use = "wilcox")

data <- add_regulate(markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 1, fdr = 0.05)
data$row <- rownames(data)
data$regulate[data$regulate=="Up"] <- "HS"
data$regulate[data$regulate=="Down"] <- "SX"
data$regulate[data$regulate=="Normal"] <- "non-significant"

# plot
pdf("figures/Fig3D.pdf", height = 5, width = 5)
ggvolcano(data, x = "log2FoldChange", y = "padj", log2FC_cut = 1, FDR_cut = 0.05, legend_position = "DL",
          colors = c(PositionColor,"non-significant"="grey"), 
          fills = c(PositionColor,"non-significant"="grey"), 
          label = "row", label_number = 20, output = FALSE)
dev.off()

## Figure 3E: key markers of S2 BMDM in  stereotactic biopsy data and trajectory analysis----------
# prepare markers
markers_df <- read.table("result/BMDM_markers200_withXHsample.txt", header = T)
markers_list <- markers_df %>%
  filter(p_val_adj < 0.05 & !grepl("^MT|^RP", gene) & avg_log2FC > 1) %>%
  split(.$cluster) %>%  
  lapply(function(x) x$gene)

intersected_genes <- intersect(markers_list$S2, data$row[data$regulate == "center"])
intersected_genes
# [1] "ALDOA"    "GAPDH"    "CSTB"     "TPI1"     "LDHA"     "RGCC"     "SPP1"     "S100A9"   "G0S2"     "TIMP1"   
# [11] "RAB42"    "C15orf48" "S100A8"  

pdf("figures/Sup_Fig3D.pdf", width = 6, height = 3)
VlnPlot(TAM, features = c("NUPR1",intersected_genes), group.by = "position", stack = T) + NoLegend()
dev.off()

# trajetory analysis starts from here
library(monocle)

DefaultAssay(TAM) <- 'RNA'
cds <- as.CellDataSet(TAM)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

sel.gene <- data$row[data$regulate!="non-significant"]
cds <- monocle::setOrderingFilter(cds, sel.gene)

cds <- monocle::reduceDimension(cds, method = 'DDRTree')
cds <- monocle::orderCells(cds)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$position)[,"SX"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 4)
diff_test_res <- diff_test_res[order(diff_test_res$qval),]

sig_gene_names <- intersect(rownames(subset(diff_test_res, qval < 0.05)),intersected_genes)
# [1] "ALDOA"    "GAPDH"    "CSTB"     "TPI1"     "LDHA"     "RGCC"     "SPP1"     "S100A9"   "G0S2"     "TIMP1"   
# [11] "RAB42"    "C15orf48" "S100A8"  

pdf("figures/Sup_Fig3E.pdf", height = 4, width = 4)
monocle::plot_cell_trajectory(cds, color_by = "position") + scale_color_manual(values = PositionColor)
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf("figures/Fig3E.pdf", height = 6, width = 6)
plot_genes_in_pseudotime(cds[sig_gene_names,], color_by = "position", ncol = 3) + scale_color_manual(values = PositionColor)
dev.off()

## Figure 3G: density plot of scores change on pseudotime----------
library(UCell)
library(tidyr)

score_list <- readRDS("result/gene_list_Fig3A.rds")
TAM <- AddModuleScore_UCell(TAM, score_list, name = NULL)
TAM$pseudotime <- cds$Pseudotime

# Gather scores into long format
plot_df <- TAM@meta.data %>%
  dplyr::select(pseudotime, position, MES, BMDM, S2, MP1, Hypoxia, Ferroptosis_neg, Ferroptosis_dual, Ferroptosis_pos) %>%
  pivot_longer(cols = -c(pseudotime, position), names_to = "signature", values_to = "score")
plot_df$signature <- factor(plot_df$signature, levels = c("BMDM", "S2", "MP1", "Hypoxia", "Ferroptosis_neg",  "Ferroptosis_dual", "Ferroptosis_pos"))

# Plot
pdf("figures/Fig3G.pdf", width = 10, height = 4)
ggplot(plot_df, aes(x = pseudotime, y = score)) +
  geom_point(aes(color = position), size = 1) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5, color = "black", size = 1) +
  facet_wrap(~signature, scales = "free_y", ncol = 4) +
  scale_color_manual(values = PositionColor) +
  theme_minimal(base_size = 13) +
  labs(x = "Pseudotime", y = "Signature Score") +
  theme(strip.text = element_text(face = "bold"))
dev.off()

VlnPlot(TAM, c("S2","MP1","Ferroptosis_neg"), group.by = "seurat_clusters")
TAM$TAM_subtype <- "other TAMs"
TAM$TAM_subtype[TAM$seurat_clusters == "1"] <- "S2 BMDM"
table(TAM$TAM_subtype, TAM$TAM_type)

pdf("figures/Sup_Fig3F.pdf", height = 4, width = 4)
DimPlot(TAM, group.by = "seurat_clusters")
VlnPlot(TAM, c("S2","MP1","Ferroptosis_neg"))
DimPlot(TAM, group.by = "TAM_subtype")
dev.off()

## Fig 3H: crosstalk from CellphoneDB----------------
merge_serial$celltype <- NA
merge_serial$celltype[match(colnames(malignant),colnames(merge_serial))] <- malignant$subtype
merge_serial$celltype[merge_serial$celltype %in% c("OPC","NPC","AC")] <- "other tumor cells"
merge_serial$celltype[match(colnames(TAM),colnames(merge_serial))] <- TAM$TAM_subtype
merge_serial_subset <- merge_serial[,!is.na(merge_serial$celltype)]
table(merge_serial_subset$celltype)

merge_serial_subset_HS <- subset(merge_serial_subset, subset = position == "HS")
merge_serial_subset_SX <- subset(merge_serial_subset, subset = position == "SX")

meta_HS <- data.frame(
  Cell = colnames(merge_serial_subset_HS),
  cell_type = merge_serial_subset_HS$celltype
)

# Save meta
write.table(meta_HS[, c("Cell", "cell_type")], "cellphonedb/meta_HS.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Idents(merge_serial_subset_HS) <- "celltype"

# Identify markers for all clusters
all_markers <- FindAllMarkers(
  object = merge_serial_subset_HS,
  only.pos = TRUE,      # Consider only upregulated genes
  min.pct = 0.25,       # Genes expressed in at least 25% of cells in either group
  logfc.threshold = 0.5, # Log-fold-change threshold
)

# Select relevant columns and rename them
degs_data <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::select(cluster, gene, avg_log2FC, p_val, p_val_adj) %>%
  dplyr::rename(cell_type = cluster, logFC = avg_log2FC, P.Value = p_val, adj.P.Val = p_val_adj)

# Save to a tab-delimited file
write.table(degs_data, "cellphonedb/degs_HS.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Save expression
counts_data <- GetAssayData(merge_serial_subset_HS, slot = "counts")
counts_df <- as.data.frame(as.matrix(counts_data))
counts_df <- cbind(Gene = rownames(counts_df), counts_df)
write.table(counts_df, "cellphonedb/counts_HS.txt", sep = "\t", row.names = FALSE, quote = FALSE)

meta_SX <- data.frame(
  Cell = colnames(merge_serial_subset_SX),
  cell_type = merge_serial_subset_SX$celltype
)

# Save meta
write.table(meta_SX[, c("Cell", "cell_type")], "cellphonedb/meta_SX.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Idents(merge_serial_subset_SX) <- "celltype"

# Identify markers for all clusters
all_markers <- FindAllMarkers(
  object = merge_serial_subset_SX,
  only.pos = TRUE,      # Consider only upregulated genes
  min.pct = 0.25,       # Genes expressed in at least 25% of cells in either group
  logfc.threshold = 0.5, # Log-fold-change threshold
)

# Select relevant columns and rename them
degs_data <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::select(cluster, gene, avg_log2FC, p_val, p_val_adj) %>%
  dplyr::rename(cell_type = cluster, logFC = avg_log2FC, P.Value = p_val, adj.P.Val = p_val_adj)

# Save to a tab-delimited file
write.table(degs_data, "cellphonedb/degs_SX.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Save expression
counts_data <- GetAssayData(merge_serial_subset_SX, slot = "counts")
counts_df <- as.data.frame(as.matrix(counts_data))
counts_df <- cbind(Gene = rownames(counts_df), counts_df)
write.table(counts_df, "cellphonedb/counts_SX.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# run cellphonedb in python 3.8

# interpret the results
library(readr)
library(ComplexHeatmap)

# Load CellPhoneDB output
sig_HS <- read_tsv("cellphonedb/outs_HS/degs_analysis_significant_means_04_12_2025_121557.txt")
sig_SX <- read_tsv("cellphonedb/outs_SX/degs_analysis_significant_means_04_12_2025_121559.txt")

get_OR_means <- function(df) {
  keep_cols <- c("MES|S2 BMDM", "MES|other TAMs", "other tumor cells|S2 BMDM", "other tumor cells|other TAMs",
                 "S2 BMDM|MES", "other TAMs|MES", "S2 BMDM|other tumor cells", "other TAMs|other tumor cells")
  
  df %>%
    select(interacting_pair, all_of(keep_cols)) %>%
    filter(if_any(all_of(keep_cols), ~ !is.na(.)))
}

HS_df <- get_OR_means(sig_HS)
SX_df <- get_OR_means(sig_SX)

# 创建矩阵并统一 row 顺序
HS_mat <- HS_df %>%
  column_to_rownames("interacting_pair")

SX_mat <- SX_df %>%
  column_to_rownames("interacting_pair")

all_rows <- union(rownames(HS_mat), rownames(SX_mat)) %>% sort()

HS_mat <- HS_mat[all_rows, , drop = FALSE]
SX_mat <- SX_mat[all_rows, , drop = FALSE]
rownames(HS_mat) <- all_rows
rownames(SX_mat) <- all_rows

HS_mat_num <- as.matrix(sapply(HS_mat, as.numeric))
SX_mat_num <- as.matrix(sapply(SX_mat, as.numeric))

merged_mat <- cbind(HS_mat_num, SX_mat_num)
merged_cols <- c(paste0("HS_",colnames(HS_mat_num)),paste0("SX_",colnames(SX_mat_num)))

max_colnames <- apply(merged_mat, 1, function(x) {
  if (all(is.na(x))) return(NA)
  merged_cols[which.max(x)]
})

keep_rows <- which(max_colnames %in% c("HS_MES|S2 BMDM", "HS_S2 BMDM|MES"))

# 筛选矩阵
HS_mat_filtered <- HS_mat[keep_rows, , drop = FALSE]
SX_mat_filtered <- SX_mat[keep_rows, , drop = FALSE]

mean_expr <- rowMeans(cbind(
  rowMeans(HS_mat_filtered, na.rm = TRUE),
  rowMeans(SX_mat_filtered, na.rm = TRUE)
), na.rm = TRUE)

# 根据平均表达值排序行名（从高到低）
row_order <- names(mean_expr)[order(mean_expr, decreasing = TRUE)]

# 重新排序两个矩阵
HS_mat_filtered <- HS_mat_filtered[row_order, , drop = FALSE]
SX_mat_filtered <- SX_mat_filtered[row_order, , drop = FALSE]

# 颜色映射
all_values <- c(as.matrix(HS_mat_filtered), as.matrix(SX_mat_filtered))
summary(all_values)
col_fun <- colorRamp2(c(0, max(all_values, na.rm = T)), c("#FFF5F0", "firebrick3"))

ht_left <- Heatmap(
  HS_mat_filtered[row_order, , drop = FALSE],
  name = "HS mean",
  col = col_fun,
  na_col = "white", 
  row_order = row_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "HS",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,                  
  column_names_side = "bottom",          
  column_names_gp = gpar(fontsize = 10),  
  rect_gp = gpar(col = "black", lwd = 0.5),
  border = gpar(col = "black", lwd = 2)
)

ht_right <- Heatmap(
  SX_mat_filtered[row_order, , drop = FALSE],
  name = "SX mean",
  col = col_fun,
  na_col = "white",
  row_order = row_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "SX",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 0.5),
  border = gpar(col = "black", lwd = 2)
)

pdf("figures/Fig3H.pdf", height = 4, width = 10)
draw(ht_left + ht_right, heatmap_legend_side = "right")
dev.off()

write.table(plot_df, file = "tables/Fig3E.tsv", quote = F, sep = "\t")
write.table(prop_data, file = "tables/Fig3G.tsv", quote = F, sep = "\t")
write.table(HS_mat_filtered[row_order, , drop = FALSE], file = "tables/Fig3H.tsv", quote = F, sep = "\t")
write.table(SX_mat_filtered[row_order, , drop = FALSE], file = "tables/Fig3I.tsv", quote = F, sep = "\t")


