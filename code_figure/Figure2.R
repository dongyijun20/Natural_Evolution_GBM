## Figure 2A: UMAP of myeloid cells-------
library(GSVA)
merge <- readRDS("result/merged_withXHsample_harmony.rds")
myeloid <- subset(merge, subset = cell_type == "Myeloid")
SampleColor <- readRDS("color/SampleColor.rds")

BMDM_MG_markers <- read_excel("reference/13059_2017_1362_MOESM5_ESM.xlsx")
BMDM_MG_markers <- as.list(BMDM_MG_markers)
names(BMDM_MG_markers) <- c("MG","BMDM")

BMDM_MG_markers <- list(BMDM=unique(c(BMDM_MG_markers$BMDM,c("CCR2","CD1C","CD1B","TGFBI","FXYD5","FCGR2B","CLEC12A","CLEC10A","CD207","CD209"))),
                        MG=unique(c(BMDM_MG_markers$MG,c("CX3CR1","SALL1","HEXB","P2RY12","TMEM119"))))
saveRDS(BMDM_MG_markers, file = "reference/BMDM_MG_markers.rds")

myeloid <- myeloid %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=1)
myeloid <- myeloid %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

Idents(myeloid)<-"seurat_clusters"
exp=AverageExpression(myeloid,assays = "RNA")
counts2=exp[["RNA"]]
GSVA_hall <- gsva(expr=counts2, 
                  gset.idx.list=BMDM_MG_markers, 
                  mx.diff=T, 
                  kcdf="Gaussian",
                  parallel.sz=4) 
head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall,
                   cluster_rows = T,
                   cluster_cols =T,
                   show_colnames = T,
                   scale = "none") 

subtype <- apply(GSVA_hall,2,function(x) rownames(GSVA_hall)[which.max(x)])
myeloid$TAM_type <- as.vector(sapply(myeloid$seurat_clusters,function(x) subtype[x]))
DimPlot(myeloid, group.by = "TAM_type")

# add other labels manually
DimPlot(myeloid, label = T)
myeloid$TAM_type[myeloid$seurat_clusters==9] <- "Proliferating" 
myeloid$TAM_type[myeloid$seurat_clusters==16] <- "Monocyte"
myeloid$TAM_type[myeloid$seurat_clusters==17] <- "DC"

saveRDS(myeloid, file = "result/myeloid_withXHsample_harmony.rds")

TAM_color<-c('#df5734','#6c408e','#537eb7','#d4c2db','#ac6894')
names(TAM_color)<-names(table(myeloid$TAM_type))

p1 <- DimPlot(myeloid, reduction = "umap", group.by = "TAM_type", label = TRUE, cols = TAM_color,
        label.size = 3, repel = TRUE) 
ggsave("figures/Fig2A.svg", plot = p1, height = 4, width = 5.5, dpi = 300)

p2 <- DimPlot(myeloid, reduction = "umap", group.by = "grade", label = TRUE, cols = GradeColor,
        label.size = 3, repel = TRUE) 
ggsave("figures/Sup_Fig2A.svg", plot = p2, height = 4, width = 5.5, dpi = 300)


## Figure 2B: signature scores of myeloid cells in patients with different evolution grade--------

# set the order of x-axis
myeloid$sample <- as.factor(myeloid$sample)

# plot boxplot of BMDM and MG scores
df_BMDM <- data.frame(BMDM = myeloid$BMDM, sample = myeloid$sample)
df_MG <- data.frame(MG = myeloid$MG, sample = myeloid$sample)
comparelist <- list(c("NJ01_1","NJ01_2"),c("NJ02_1","NJ02_2"),c("TT01_1","TT01_2"),c("TT02_1","TT02_2"),c("XH01_1","XH01_2")) 

pdf("figures/Fig2B.pdf", height = 4, width = 5)
ggboxplot(df_BMDM, x = "sample", y = "BMDM",
          color = "sample", palette = SampleColor)+ ylab("BMDM score") +
  stat_compare_means(comparisons = comparelist, label = "p.signif", label.y = max(df_BMDM$BMDM)+0.05) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()
dev.off()

pdf("figures/Sup_Fig2B.pdf", height = 4, width = 5)
ggboxplot(df_MG, x = "sample", y = "MG",
          color = "sample", palette = SampleColor)+ ylab("MG score")+
  stat_compare_means(comparisons = comparelist, method = "wilcox.test", label = "p.signif",  label.y = max(df_MG$MG)+0.05) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()
dev.off()

pdf("figures/Sup_Fig2C.pdf", width = 20, height = 5)
FeaturePlot(myeloid, features = c("BMDM", "MG"), max.cutoff = "q90", min.cutoff = "q10", blend = T) +
  ggtitle("BMDM score and MG score")
dev.off()

## Figure 2C: subtypes of BMDM, UMAP-----------
BMDM <- subset(myeloid, subset = TAM_type == "BMDM")
BMDM <- BMDM %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution=0.2)
BMDM <- BMDM %>%
  RunUMAP(dims = 1:40, reduction = "harmony")
DimPlot(BMDM, group.by = "seurat_clusters")
DimPlot(BMDM, group.by = "grade")

BMDM$subcluster <- paste0("S_", BMDM$seurat_clusters)
DimPlot(BMDM, group.by = "subcluster")

science_myeloid <- readRDS("reference/GSE278456_MyeloidSeurat.RDS")
science_myeloid$celltype_science <- Idents(science_myeloid)
science_myeloid_subset <- subset(science_myeloid, subset = Tumor.Grade=="IV" & !grepl("MCG", celltype_science))
science_myeloid_subset <- science_myeloid_subset[,sample(colnames(science_myeloid_subset),10000)]
DimPlot(science_myeloid_subset)
table(science_myeloid_subset$celltype_science)

anchors <- FindTransferAnchors(reference = science_myeloid_subset, query = BMDM, dims = 1:30, reference.reduction = "pca")
BMDM <- MapQuery(
  anchorset = anchors,
  query = BMDM,
  reference = science_myeloid_subset,
  refdata = list(
    celltype_science = "celltype_science"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)
DimPlot(BMDM, group.by = "predicted.celltype_science", label = T)

table(BMDM$subcluster)
BMDM$subcluster_named <- BMDM$subcluster
BMDM$subcluster_named[BMDM$subcluster_named=="S0"] <- "ELMO1+ MAC"
BMDM$subcluster_named[BMDM$subcluster_named=="S3"] <- "HLA-hi MAC"
BMDM$subcluster_named[BMDM$subcluster_named=="S2"] <- "NUPR1+ MDSC"
BMDM$subcluster_named[BMDM$subcluster_named=="S1"] <- "VCAN+ MDSC"
BMDM$subcluster_named[BMDM$subcluster_named=="S4"] <- "CXCL8+ MDSC"

SubclusterColor <- c('#ea9994','#f2c396','#9cd2ed','#86c7b4','#a992c0')
names(SubclusterColor) <- unique(BMDM$subcluster_named)

pdf("figures/Fig2C.pdf", height = 4, width = 5)
DimPlot(BMDM, reduction = "umap", group.by = "subcluster_named", label = TRUE, cols = SubclusterColor,
        label.size = 3, repel = TRUE) 
dev.off()

saveRDS(BMDM, file = "result/BMDM_withXHsample_harmony.rds")

## Figure 2D: marker genes dotplot of different BMDM clusters----------
Idents(BMDM) <- "subcluster"
markers <- FindAllMarkers(BMDM, only.pos = TRUE, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_df <- markers_df[markers_df$p_val_adj<0.05,]
write.table(markers_df,file="result/BMDM_markers200_withXHsample.txt",quote=F,sep="\t",row.names=F,col.names=T)

markers_df <- read.table("result/BMDM_markers200_withXHsample.txt", header = T)
markers_df <-  markers_df %>% 
  arrange(cluster) %>% 
  dplyr::filter(!grepl("^MT|^RP", gene) & avg_log2FC > 1) 
markers_filter <-  markers_df %>% 
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  slice_head(n = 10)

BMDM$subcluster <- factor(BMDM$subcluster, levels = c("S0", "S1", "S2", "S3", "S4"))
Idents(BMDM) <- "subcluster"

pdf("figures/Fig2D.pdf", height = 8, width = 5)
DotPlot(BMDM, features = unique(markers_filter$gene), scale = T, group.by = "subcluster") + 
  coord_flip() + 
  theme(panel.grid = element_blank(), 
        axis.text.y=element_text(size = 10))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression"))+ 
  scale_colour_gradient2(low = "navy", high = "firebrick3")
dev.off()

MDSC_markers <- list("E-MDSC" = c("GPNMB","SPP1","MT2A","NUPR1","LGALS1","PLIN2","CSTB",
                                  "FTL","MT1X","VIM","MT1G","GAPDH","CTSD","FTH1","LGALS3",
                                  "MIF","BNIP3","HMOX1","ENO1"),
                     "M-MDSC" = c("LYZ","VIM","S100A6","VCAN","S100A10","CXCL2","TIMP1",
                                  "S100A4","EREG","CXCL3","FN1","ANXA1","CD44","LGALS1",
                                  "S100A9","ANXA2","CXCL8","LGALS3","AHNAK"),
                     "pseudotime-E" = c("SLC2A1","HK2","ENO2","GAPDH","SCD","BNIP3","NUPR1",
                                        "FTL","MT1H","SLC2A5","GPNMB","MT1G","MT1X","MT2A","HMOX1",
                                        "APOC1"))
markers <- readRDS("result/gene_list_Fig3A.rds")

pdf("figures/Sup_Fig4F.pdf", height = 6, width = 4)
VlnPlot(science_myeloid_subset, intersect(markers$S2, MDSC_markers$`E-MDSC`), stack = T, flip = T, pt.size = 0) + NoLegend()
dev.off()

## Figure 2E: subtype signature scores of BMDM cells in patients with different evolution grade--------
library(UCell)
library(ggpubr)

markers_df <- read.table("result/BMDM_markers200_withXHsample.txt", header = T)
markers_list <- markers_df %>%
  filter(p_val_adj < 0.05 & !grepl("^MT|^RP", gene) & avg_log2FC > 1) %>%
  split(.$cluster) %>%  
  lapply(function(x) x$gene)
saveRDS(markers_list, file = "result/BMDM_markers_filter.rds")

BMDM <- AddModuleScore_UCell(BMDM, features = markers_list, name = NULL)

# set the order of x-axis
BMDM$sample <- as.factor(BMDM$sample)

# plot boxplot of S2 scores
df_S2 <- data.frame(S2 = BMDM$S2, sample = BMDM$sample)
comparelist <- list(c("NJ01_1","NJ01_2"),c("NJ02_1","NJ02_2"),c("TT01_1","TT01_2"),c("TT02_1","TT02_2"),c("XH01_1","XH01_2")) 

pdf("figures/Fig2E.pdf", height = 4, width = 5)
ggboxplot(df_S2, x = "sample", y = "S2",
          color = "sample", palette = SampleColor)+ ylab("S2 score") +
  stat_compare_means(comparisons = comparelist, method = "wilcox.test", label = "p.signif", label.y = max(df_S2$S2) + 0.05) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()
dev.off()

pdf("figures/Sup_Fig2D.pdf", height = 5, width = 4)
VlnPlot(BMDM, features = c(levels(BMDM$subcluster),"NUPR1"), group.by = "sample", 
        pt.size = 0, stack = T, cols = c(SubclusterColor,"NUPR1"='#bb82b1'), flip = T) + NoLegend()
dev.off()

## correlation in TCGA and CGGA dataset
expression_TCGA <- readRDS("reference/TCGA/TCGA_GBM_expression_data.rds")
fpkm_data_TCGA <- assay(expression_data, "fpkm_unstrand")
rownames(fpkm_data_TCGA) <- rowData(expression_data)$gene_name

fpkm_data_CGGA <- readRDS("reference/CGGA/CGGA_GBM_expression_data.rds")
common_genes <- intersect(rownames(fpkm_data_TCGA), rownames(fpkm_data_CGGA))
expr_tcga_filtered <- fpkm_data_TCGA[common_genes, ]
expr_cgga_filtered <- fpkm_data_CGGA[common_genes, ]
combined_expr <- cbind(expr_tcga_filtered, expr_cgga_filtered)

selected_genes <- intersect(marker_list$S2, rownames(combined_expr))
expr_subset <- combined_expr[selected_genes, ]
cor_matrix <- cor(t(expr_subset), method = "pearson")

library(pheatmap)
breaks_list <- seq(-1, 1, length.out = 100)

pdf("figures/Sup_Fig4G.pdf")
pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(99),  # one less than breaks
         breaks = breaks_list,
         main = "Gene-Gene Correlation (TCGA + CGGA FPKM)")
dev.off()

## Figure 2F: Function of BMDM subclusters--------
library(clusterProfiler)
library(org.Hs.eg.db)

sce.markers <- read.table("result/BMDM_markers200_withXHsample.txt", header = T)
sce.markers <-  sce.markers %>% 
  arrange(cluster) %>% 
  dplyr::filter(!grepl("^MT|^RP", gene)) 

ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

GO_results <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05
)

top_descriptions <- GO_results@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(order_by = GeneRatio, n = 6) %>% 
  ungroup() %>%
  distinct(Description)

GO_results_filter <- GO_results@compareClusterResult %>%
  mutate(Description = str_wrap(Description, width = 150)) %>%
  filter(Description %in% top_descriptions$Description)

GO_results_S2 <-GO_results %>% filter(Cluster=="S2")
GO_results_S2 <- pairwise_termsim(GO_results_S2)

pdf("figures/Sup_Fig2F_2.pdf", height = 4, width = 7)
emapplot(GO_results_S2, cex.params = list(cex_label_category = 2, line = 0.5))
dev.off()

pdf("figures/Sup_Fig2D_2.pdf", height = 8, width = 5)
dotplot(GO_results_filter, showCategory = 6) + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5))
dev.off()

## Figure 2G: NMF analysis of all myeloid cells----------
library(GeneNMF)

seu.list <- SplitObject(myeloid, split.by = "orig.ident")

geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)

set.seed(122)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        weight.explained = 0.7,
                                        nMP = 6,
                                        max.genes=100)

geneNMF.metaprograms$metaprograms.metrics
mp.genes <- geneNMF.metaprograms$metaprograms.genes
grep("NUPR1",mp.genes)
# intersect(mp.genes$MP1,sce.markers$gene[sce.markers$cluster=="S2"])
# [1] "S100A10"  "C15orf48" "S100A9"   "PLIN2"    "S100A6"   "NUPR1"    "RGCC"     "RAB42"    "LGALS3"   "EMP3"     "PLP2"     "LSP1"     "SNHG12"  
# [14] "VIM"      "UPP1"     "LGALS1"   "CSTB"     "FABP5"    "MIF"  

saveRDS(mp.genes, file = "result/mp.genes.rds")

plotMetaPrograms(geneNMF.metaprograms,
                 similarity.cutoff = c(0.1,0.9))

geneNMF.metaprograms$group <- sapply(names(geneNMF.metaprograms$programs.clusters), function(name) {
  if (grepl("_1", name)) {
    "older"
  } else if (grepl("_2", name)) {
    "younger"
  }
})

## custom heatmap
library(pheatmap)
library(viridis)

# Extract data
J <- geneNMF.metaprograms[["programs.similarity"]]
tree <- geneNMF.metaprograms[["programs.tree"]]
cl_members <- geneNMF.metaprograms[["programs.clusters"]]
labs.order <- labels(as.dendrogram(tree))
group_vector <- geneNMF.metaprograms[["group"]]

# Prepare cl_members
cl_names <- names(cl_members)
cl_members[!is.na(cl_members)] <- paste0("MP", cl_members[!is.na(cl_members)])
names(cl_members) <- cl_names

# Recover order of MP clusters
cluster.order <- unique(cl_members[labs.order])
nMP <- length(cluster.order)

# Gaps for heatmap
diffs <- diff(as.integer(as.factor(cl_members)))
gaps <- which(diffs != 0)

# Define annotation column
annotation_col <- data.frame(
  Group = factor(group_vector, levels = c("older", "younger")),
  Metaprogram = factor(cl_members, levels = cluster.order)
)

# Define colors
MetaprogramColor <- setNames(viridis(nMP, option = "D"), levels(annotation_col$Metaprogram))  # Use Viridis for metaprograms

# Define annotation colors
annotation_colors <- list(
  Group = GradeColor,
  Metaprogram = MetaprogramColor
)

# Apply similarity cutoff for the heatmap
similarity.cutoff <- c(0.1, 0.9)
J[J < similarity.cutoff[1]] <- similarity.cutoff[1]
J[J > similarity.cutoff[2]] <- similarity.cutoff[2]

# Define heatmap palette
heatmap_palette <- viridis(100, option = "A", direction = -1)

# Generate the heatmap
pdf("figures/Fig2G.pdf", height = 5.5, width = 6)
pheatmap(
  J,
  scale = "none",
  color = heatmap_palette,
  main = "Clustered Heatmap",
  cluster_rows = tree,
  cluster_cols = tree,
  cutree_rows = nMP,
  cutree_cols = nMP,
  gaps_col = gaps,  
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE
)
dev.off()

seu <- AddModuleScore_UCell(myeloid, features = mp.genes, assay="RNA", name = NULL)

pdf("figures/Sup_Fig2E.pdf", height = 5, width = 4)
VlnPlot(seu, features=names(mp.genes), group.by = "TAM_type", 
        pt.size = 0, stack = T, cols = MetaprogramColor, flip = T) + NoLegend()
dev.off()

## Figure 2H: the relationship between S2, S3 BMDM markers and MP1 metaprograms, the distribution of MP1 across samples----------
BMDM <- AddModuleScore_UCell(BMDM, features = mp.genes, name = NULL)

# plot boxplot of MP scores
df_MP1 <- data.frame(MP1 = BMDM$MP1, sample = BMDM$sample)
comparelist <- list(c("NJ01_1","NJ01_2"),c("NJ02_1","NJ02_2"),c("TT01_1","TT01_2"),c("TT02_1","TT02_2"),c("XH01_1","XH01_2")) 

pdf("figures/Fig2H.pdf", height = 4, width = 5)
ggboxplot(df_MP1, x = "sample", y = "MP1",
          color = "sample", palette = SampleColor)+ ylab("MP1 score") +
  stat_compare_means(comparisons = comparelist, method = "wilcox.test", label = "p.signif", label.y = max(df_MP1$MP1) + 0.05) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()
dev.off()

data_cor <- data.frame(
  MP1 = BMDM$MP1,
  S2 = BMDM$S2,
  subcluster = BMDM$subcluster
)

correlation <- cor(data_cor$MP1, data_cor$S2)
correlation

pdf("figures/Sup_Fig2F.pdf", height = 5, width = 6)
ggplot(data_cor, aes(x = MP1, y = S2, color = subcluster)) +
  geom_point(alpha = 0.8, size = 0.8)  +
  scale_color_manual(values = SubclusterColor) +           
  geom_smooth(method = "lm", color = "lightgrey", se = FALSE) +  
  labs(
    title = "Correlation Plot Between MP1 and S2 Scores",
    subtitle = paste("Correlation Coefficient:", round(correlation, 2)),
    x = "MP1 Score",
    y = "S2 Score"
  ) +
  theme_minimal(base_size = 15)
dev.off()

pdf("figures/Sup_Fig2G.pdf", height = 5, width = 4)
VlnPlot(BMDM, feature = paste0("MP",1:6), group.by = "subcluster", 
        pt.size = 0, stack = T, cols = MetaprogramColor, flip = T) + NoLegend()
dev.off()

## Figure 2I: functional enrichment------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)

mp_list <- lapply(mp.genes, function(x){
  bitr(x, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db', drop = T)[,2]
})
names(mp_list) <- names(mp.genes)

GO_MP1 <- enrichGO(gene = mp_list$MP1, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
GO_MP1 <- pairwise_termsim(GO_MP1)

pdf("figures/Sup_Fig2I.pdf", height = 5, width = 6)
emapplot(GO_MP1, cex.params = list(category_label = 0.5, line = 0.5))
dev.off()

## Figure 2J: relationship with regulators of ferroptosis----------
ferroptosis <- list(
  negative = c("NUPR1", "FTL", "FTH1", "GPX4", "SLC3A2", "NFE2L2", "SAT1", "ZFP36"),
  dual = c("ATF3", "ATF4", "CSTB", "ELAVL1", "HMOX1", "SP1", "TP53"),
  positive = c("ATG5", "ATG7", "BECN1", "STAT3", "ACSL4", "TFRC", "NCOA4", "IREB2")
)

saveRDS(ferroptosis, file = "reference/ferroptosis_markers.rds")

feature_category <- stack(ferroptosis) %>%
  setNames(c("feature", "category"))

category_breaks <- c(0,length(ferroptosis$negative)+0.5,length(ferroptosis$negative)+length(ferroptosis$dual)+0.5)

dotplot <- DotPlot(
  BMDM,
  features = feature_category$feature,
  scale = TRUE, 
  group.by = "subcluster"
)

pdf("figures/Fig2J.pdf", height = 6, width = 5)
dotplot +
  coord_flip() + 
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_blank()
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend("Percent Expression")) +
  scale_colour_gradient2(low = "navy", high = "firebrick3") +
  annotate(
    "segment",
    x = category_breaks, # Adjust for feature groups
    xend = category_breaks,
    y = -Inf,
    yend = Inf,
    color = "black",
    linetype = "dashed"
  ) + annotate(
    "text",
    x = category_breaks+3, # Adjust positions for labels
    y = 6,
    label = c("Negative", "Dual", "Positive"),
    size = 4,
    hjust = 0
  ) + expand_limits(y = 8.5)
dev.off()

## scoring of signatures and their correlation in BMDM---------------
library(pheatmap)

ferroptosis <- reqdRDS("reference/ferroptosis_markers.rds")
mp.genes <- readRDS("result/mp.genes.rds")
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

score_list <- list(Hypoxia = hallmarks$gene_symbol[hallmarks$gs_name=="HALLMARK_HYPOXIA"],
                   Ferroptosis_neg = ferroptosis$negative,
                   Ferroptosis_pos = ferroptosis$positive,
                   Ferroptosis_dual = ferroptosis$dual
)

BMDM <- AddModuleScore_UCell(BMDM, score_list, name = NULL)

scores_df <- BMDM@meta.data[, c("S2","MP1",names(score_list))]

correlation_matrix <- cor(scores_df, method = "pearson")

pdf("figures/Fig2K.pdf", height = 5, width = 5.5)
pheatmap(correlation_matrix, 
         display_numbers = TRUE, 
         clustering_method = "complete", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Correlation of Signatures in scRNA-seq BMDM (n = 10)")
dev.off()

## Figure 2I (modified): emapplot of different gene signatures---------

older_markers <- FindMarkers(BMDM, ident.1 = "older", ident.2 = "younger", group.by = "grade", logfc.threshold = 0.3, only.pos = T)

gene_list <- list(older = rownames(older_markers)[older_markers$p_val_adj<0.05], 
                  S2 = markers_list$S2, 
                  MP1 = mp.genes$MP1, 
                  ferroptosis_neg = ferroptosis$negative)

names_gene_list <- names(gene_list)
gene_list <- lapply(names, function(x) bitr(gene_list[[x]],'SYMBOL','ENTREZID','org.Hs.eg.db', drop = T)[,2])
names(gene_list) <- names_gene_list

compare_results <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",  
  pvalueCutoff = 0.01
)

dotplot(compare_results)

compare_results <- pairwise_termsim(compare_results)

pdf("figures/Sup_Fig2J.pdf", height = 15, width = 18)
emapplot(
  compare_results,
  layout.params = list(layout = "nicely"),
  cex.params = list(
    category_label = 1,
    line = 0.2)
)
dev.off()

pdf("figures/Fig2I.pdf", height = 6, width = 7)
emapplot(
  compare_results, 
  layout.params = list(layout = "nicely"),
  cex.params = list(
    category_label = 0,
    line = 0.2)
)
dev.off()
