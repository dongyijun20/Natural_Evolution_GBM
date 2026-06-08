## Figure 1B-C: Umap of all cells-------
library(RColorBrewer)
library(ArchR)
library(ggunchull)
library(RColorBrewer)
library(Seurat)

# define colors
SampleColor <- brewer.pal(10, 'Paired')[c(2,1,4,3,6,5,8,7,10,9)]
names(SampleColor) <- c("NJ02_1","NJ02_2","TT01_1","TT01_2","TT02_1","TT02_2","NJ01_1","NJ01_2","XH01_1","XH01_2")

GradeColor <- c("younger" = "lightblue",
                "older" = "orange")
saveRDS(GradeColor, file = "color/GradeColor.rds")

# this is the seurat object after QC, normalization, dimension reduction, and harmony
merge <- readRDS("result/merged_withXHsample_harmony.rds")

# delete cells within cluster <= 200 cells, and delete cells that express dual markers (malignant and non-malignant, potentially doublets)
remove_cluster <- names(table(merge$seurat_clusters))[table(merge$seurat_clusters)<=200]
remove_cluster <- c(remove_cluster, 17)
remove_cluster
merge <- merge[,!merge$seurat_clusters %in% remove_cluster]

# adjust resolution of cluster to a suitable value
merge <- FindClusters(merge, resolution = 0.3) # rough classification

# determine cell types by markers
pdf("figures/Sup_Fig1A.pdf", height = 4, width = 5.5)
DimPlot(merge, group.by="seurat_clusters", label = T) + NoLegend()
DimPlot(merge, group.by = "sample", cols = SampleColor)
DimPlot(merge, group.by = "grade", cols = GradeColor)
dev.off()

features_2test <- c("SOX2","SOX4",
                    "MAG","PLLP",
                    "CD3E","PTPRC",
                    "CD14","CD68",
                    "COL1A1","COL1A2",
                    "CD79A",
                    "CDH5")
pdf("figures/Sup_Fig1B.pdf", height = 8, width = 8)
FeaturePlot(merge, feature = features_2test,max.cutoff = "q95", min.cutoff = "q5", ncol = 3, cols = c("grey","firebrick3"), order = T, raster = T)
VlnPlot(merge, feature = features_2test, stack = T, flip = T)
dev.off()

# give labels to cell types by seurat clusters
merge$cell_type[merge$seurat_clusters %in% c(1,2,4,5,7,10,11,12,13,15)] <- "Malignant"
merge$cell_type[merge$seurat_clusters %in% c(0,3,6,14,17)] <- "Myeloid"
merge$cell_type[merge$seurat_clusters %in% c(9)] <- "T Cell"
merge$cell_type[merge$seurat_clusters %in% c(16)] <- "Stromal"
merge$cell_type[merge$seurat_clusters %in% c(8)] <- "Oligodendrocyte"

merge$tumor_ident <- "Non-malignant"
merge$tumor_ident[merge$cell_type=="Malignant"] <- "Malignant"

table(merge$cell_type)
table(merge$tumor_ident)

# genereate consistent cell color
CellColor <- ArchRPalettes$stallion[1:length(unique(merge$cell_type))]
names(CellColor)<-unique(merge$cell_type)

p1 <- DimPlot(merge, group.by = "grade", cols = GradeColor, pt.size = 0.1)
ggsave("figures/Fig1B.svg", plot = p1, height = 4, width = 5, dpi = 300)
p2 <- DimPlot(merge, group.by = "cell_type", cols = CellColor, pt.size = 0.1, shuffle = T)
ggsave("figures/Fig1C.svg", plot = p2, height = 4, width = 5.5, dpi = 300)

## summarize the proportion of cell types in two lesions------
library(dplyr)
library(tidyr)

cell_counts_per_sample <- merge@meta.data %>%
  group_by(sample, grade, cell_type) %>%
  summarise(n = n(), .groups = "drop")

cell_props <- cell_counts_per_sample %>%
  pivot_wider(names_from = cell_type, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  ungroup()
cell_props <- cell_props %>%
  mutate(across(-c(sample, grade, total), ~ . / total))

mean_row <- cell_props %>%
  group_by(grade) %>%
  summarise(across(-sample, mean, .names = "{.col}"), .groups = "drop") %>%
  mutate(sample = paste0("mean_", grade)) %>%
  select(names(cell_props)) 

sd_row <- cell_props %>%
  group_by(grade) %>%
  summarise(across(-sample, sd, .names = "{.col}"), .groups = "drop") %>%
  mutate(sample = paste0("sd_", grade)) %>%
  select(names(cell_props)) 

summary_table <- bind_rows(cell_props, mean_row, sd_row) %>%
  arrange(grade, match(sample, c(unique(cell_props$sample), paste0("mean_", grade), paste0("sd_", grade))))

write.table(summary_table, "tables/TableS1.tsv", quote = F)

## Figure 1D: Umap of tumor cells-------
library(readxl)
library(GSVA)
library(ggsci)

malignant <- subset(merge, subset = tumor_ident == "Malignant")

# load markers of different subtypes
markers<-as.list(read_excel("reference/mmc2.xlsx", skip = 3))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1/S","G2/M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)

# run normalization and dimensional reduction on malignant cells
malignant <- malignant %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>% 
  RunUMAP(dims = 1:40, reduction = "harmony") %>% 
  FindNeighbors(dims = 1:40, reduction = "harmony") %>% 
  FindClusters(resolution = 1) # adjust a suitable resolution

# run GSVA, averaged by seurat clusters
exp=AverageExpression(malignant,assays = "RNA", group.by = "seurat_clusters")
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

p3 <- DimPlot(malignant, group.by = "subtype", cols = SubtypeColor)
ggsave("figures/Fig1D.svg", plot = p3, height = 4, width = 5, dpi = 300)

## Figure 1E: proportion of tumor subtypes-------
library(msigdbr)
library(scrabble)
library(gridExtra)
library(ggrastr)

CellInfo=malignant@meta.data
Idents(malignant)="subtype"

Counts <- as.matrix(GetAssayData(object = malignant, assay = "RNA", slot = "counts"))

#get scores
subtypelist=lapply(names(markers),function(x) markers[[x]][markers[[x]]%in%rownames(Counts)])
names(subtypelist)<-names(markers)

#there appears to be a type-o somewhere in scrabble; need to run this function to get it to run  - taken from the source code
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}

sc <- score(Counts,
         groups=subtypelist,
         binmat = NULL,
         bins = NULL,
         controls = NULL,
         bin.control = F,
         center = T,
         nbin = 30,
         n = 100,
         replace = T)

h=hierarchy (sc, quadrants = c("AC","OPC","MES","NPC"), log.scale = T)

for( i in 1: length(names(SampleColor))){
  name=names(SampleColor)[i]
  color=SampleColor[i]
  CellInfo$SampleColor[CellInfo$sample== name] <- color
}
SampleColor_vector <- as.character(CellInfo$SampleColor)
names(SampleColor_vector) <- row.names(CellInfo)
malignant <- AddMetaData(object = malignant, metadata = SampleColor_vector ,col.name = 'SampleColor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

groups=malignant@meta.data[,c("grade","subtype")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
head(matrix)
row.names(matrix)=matrix$Row.names
x=matrix$grade
y=matrix$subtype
col <- unique(x)
# names(col) <- unique(y)
matrix=matrix[,-1]
head(matrix)

title="All Glioma Samples"
p0 <- ggplot(matrix, aes(x = X, y =Y, color = grade))+geom_point_rast(size=0.3)+ 
  scale_color_manual(values=GradeColor,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour=GradeColor)+theme(legend.position = "none")+ 
  theme(legend.text=element_text(size=15),legend.title =element_text(size=15,face="bold"))+ 
  theme(axis.title.x = element_text(hjust = 0, vjust=-2, colour="black",size=10,face="bold"))+ 
  theme(axis.title.y = element_text(hjust = 0, vjust=4, colour="black",size=10,face="bold"))+  
  scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
  theme(panel.background =element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  ggtitle(title)+theme(plot.title = element_text(size=10,face="bold"))+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  geom_hline(yintercept=0, color = "black", size=0.5)+
  geom_vline(xintercept=0, color = "black", size=0.5)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  + 
  annotate("text",x = xmin+2.5, y = ymax-1,label = "MES-Like",color="white",fontface="bold",size=5)+ 
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=5)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=5)+ 
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=5)  
P0=ggplotGrob(p0)
P0

Final <- list()
Final[[1]] <- P0
for (i in 1:length(names(SampleColor))) {
  ClusterMD=malignant@meta.data[malignant@meta.data$SampleColor==SampleColor[i],]
  groups=ClusterMD[,c("SampleColor","subtype")]
  title=paste0(names(SampleColor)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$SampleColor[is.na(matrix$SampleColor)] <- "gray"
  matrix$subtype=as.character(matrix$subtype)
  matrix$subtype[is.na(matrix$subtype)] <- "Other"
  row.names(matrix)=matrix$Row.names
  matrix=matrix%>%arrange(subtype)
  x=matrix$SampleColor
  y=matrix$subtype
  col= unique(x)
  #names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(SampleColor)))+geom_point_rast(size=0.1)+geom_point_rast(data = subset(matrix, subtype !="Other"),size=0.5)+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=10,face="bold") )+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  + 
    annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=3)+ 
    annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=3)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
    annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=3)+ 
    annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=3)  
  Final[[i+1]] = ggplotGrob(P)
}
numofplots= length(Final)
layout_matrix <- rbind(
  c(1, 8, 2, 4, 6, 10),
  c(1, 9, 3, 5, 7, 11)
)

pdf("figures/Fig1E.pdf", height = 4, width = 12)
grid.arrange(grobs=Final, widths = c(2,1,1,1,1,1),layout_matrix = layout_matrix) 
dev.off()

## Figure 1F: GSEA of differential genes between older and younger tumor cells--------
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

Idents(malignant) <- "grade"
markers_malignant <- FindMarkers(malignant, ident.1 = "older",ident.2 = "younger", test.use = "wilcox")

markers_malignant <- markers_malignant %>%
  dplyr::mutate(gene = rownames(markers_malignant),
                geneID = bitr(rownames(markers_malignant), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = F)[-2251,2]
  )

gene4enrich <- markers_malignant %>% 
    dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.3) %>%
    na.omit() %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(avg_log2FC, geneID)
  
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- hallmarks[, c("gs_name", "entrez_gene")]

# Perform GSEA
gsea_results <- GSEA(gene4enrich, TERM2GENE = hallmark_list, pvalueCutoff = 0.05)

dotplot(gsea_results)
gsea_results <- pairwise_termsim(gsea_results)
emapplot(gsea_results)

pdf("figures/Fig1F.pdf", width = 5.5, height = 5)
dotplot(gsea_results, font.size = 10)
dev.off()

go_result <- enrichGO(gene = names(gene4enrich), OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
go_result <- pairwise_termsim(go_result)
kegg_result <- enrichKEGG(gene = names(gene4enrich), organism = "hsa", pvalueCutoff = 0.05)
kegg_result <- pairwise_termsim(kegg_result)

pdf("figures/Sup_Fig1D.pdf", width = 6, height = 4)
dotplot(go_result, showCategory = 10, font.size = 10)
emapplot(go_result, cex_label_category = 0.5)
dotplot(kegg_result, showCategory = 10, font.size = 10)
emapplot(kegg_result, cex_label_category = 0.5)
dev.off()

## Figure 1G: MES and NES signature scores of different samples--------
library(AUCell)
library(ggpubr)
library(patchwork)
library(tidyr)

gene_sets <- markers
NES = list(c('S100A10','FOSL2','SPP1','CAV1','ANXA1','VIM','CD44','SERPINH1',
             'LGALS3','CEBPB','ATF5','LGALS1'))
gene_sets[["NES"]] <- NES[[1]]
cells_rankings <- AUCell_buildRankings(malignant[["RNA"]]$data)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

malignant$AUC_AC  <- as.numeric(getAUC(cells_AUC)["AC", ])
malignant$AUC_MES  <- as.numeric(getAUC(cells_AUC)["MES", ])
malignant$AUC_NPC  <- as.numeric(getAUC(cells_AUC)["NPC", ])
malignant$AUC_OPC  <- as.numeric(getAUC(cells_AUC)["OPC", ])
malignant$AUC_NES  <- as.numeric(getAUC(cells_AUC)["NES", ])

# For plot: 10 boxplots per panel (cell-level); sample order and color from SampleColor
sample_levels <- c("NJ02_1","NJ02_2","TT01_1","TT01_2","TT02_1","TT02_2","NJ01_1","NJ01_2","XH01_1","XH01_2")
df_plot <- data.frame(
  sample = factor(as.character(malignant$sample), levels = sample_levels),
  MES = malignant$AUC_MES,
  NES = malignant$AUC_NES
)

# Patient-level stats for caption (LMM only)
df_scores <- malignant@meta.data %>%
  dplyr::mutate(
    sample = as.character(sample),
    patient = sub("_.*", "", sample),
    lesion = sub(".*_", "", sample)
  ) %>%
  dplyr::group_by(patient, lesion) %>%
  dplyr::summarise(
    MES = mean(AUC_MES, na.rm = TRUE),
    NES = mean(AUC_NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lesion = factor(lesion, levels = c("1", "2")))

df_wide <- df_scores %>%
  tidyr::pivot_wider(id_cols = patient, names_from = lesion, values_from = c(MES, NES))

p_t_mes <- t.test(df_wide$MES_1, df_wide$MES_2, paired = TRUE)$p.value
p_wilcox_mes <- suppressWarnings(wilcox.test(df_wide$MES_1, df_wide$MES_2, paired = TRUE)$p.value)

p_t_nes <- t.test(df_wide$NES_1, df_wide$NES_2, paired = TRUE)$p.value
p_wilcox_nes <- suppressWarnings(wilcox.test(df_wide$NES_1, df_wide$NES_2, paired = TRUE)$p.value)

if (!requireNamespace("lmerTest", quietly = TRUE)) {
  stop("Package 'lmerTest' is required for LMM-based p-value calculation.")
}
meta_ll <- malignant@meta.data
meta_ll$patient <- sub("_.*", "", as.character(meta_ll$sample))
meta_ll$lesion <- factor(sub(".*_", "", as.character(meta_ll$sample)), levels = c("1", "2"))

lmm_mes <- lmerTest::lmer(AUC_MES ~ lesion + (1|patient), data = meta_ll)
lmm_nes <- lmerTest::lmer(AUC_NES ~ lesion + (1|patient), data = meta_ll)

p_lmm_mes <- summary(lmm_mes)$coefficients["lesion2", "Pr(>|t|)"]
p_lmm_nes <- summary(lmm_nes)$coefficients["lesion2", "Pr(>|t|)"]

caption_mes <- sprintf("n=%d pairs\nLMM %s\nPaired t %s | Wilcoxon %s", 
                       n_pairs, fmt_p(p_lmm_mes), fmt_p(p_t_mes), fmt_p(p_wilcox_mes))

caption_nes <- sprintf("n=%d pairs\nLMM %s\nPaired t %s | Wilcoxon %s", 
                       n_pairs, fmt_p(p_lmm_nes), fmt_p(p_t_nes), fmt_p(p_wilcox_nes))

p_mes <- ggpubr::ggboxplot(df_plot, x = "sample", y = "MES", color = "sample", fill = "sample", palette = SampleColor) +
  labs(x = NULL, y = "MES score", subtitle = caption_mes) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none",
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

p_nes <- ggpubr::ggboxplot(df_plot, x = "sample", y = "NES", color = "sample", fill = "sample", palette = SampleColor) +
  labs(x = NULL, y = "NES score", subtitle = caption_nes) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none",
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

p_combined <- p_mes + p_nes +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))

pdf("figures/Fig1G_review.pdf", height = 4, width = 10)
print(p_combined)
dev.off()

pdf("figures/Sup_Fig1C.pdf", width = 20, height = 5)
FeaturePlot(malignant, features = c("AUC_NES", "AUC_MES"), max.cutoff = "q90", min.cutoff = "q10", blend = T)
dev.off()

## Figure Sup_Fig1E: stem-like possibilities in tumor cells----------
library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)
library(lme4) # 必须加载 lme4 以使用 GLMM

# ====================================================================
# 1 & 2. 依然使用 UCell 评分 (消除文库大小假阳性)
# ====================================================================
# 假设 sig_list 已经按照前面的逻辑建好了
malignant <- AddModuleScore_UCell(malignant, features = sig_list, name = "_UCell")
meta <- malignant@meta.data

# ====================================================================
# 3. 动态卡点 (Top 25%) 与严格互斥分类
# ====================================================================
top_pct <- 0.25 
stem_cut <- quantile(meta$StemAssoc_UCell, 1 - top_pct)
diff_cut <- quantile(meta$Diff_UCell, 1 - top_pct)

meta$is_stem <- meta$StemAssoc_UCell > stem_cut
meta$is_diff <- meta$Diff_UCell > diff_cut

# 定义互斥的逻辑列 (用于后续 GLMM 回归)
meta$is_stem_only <- meta$is_stem & !meta$is_diff
meta$is_diff_only <- meta$is_diff & !meta$is_stem
meta$is_both      <- meta$is_stem & meta$is_diff

# 计算用于画图的比例数据
df_plot <- meta %>%
  mutate(patient = sub("(.*)\\_.*","\\1",sample)) %>% 
  dplyr::group_by(patient, grade) %>%
  summarise(
    total_cells = n(),
    Stem_associated = sum(is_stem_only) / total_cells * 100,
    Differentiated  = sum(is_diff_only) / total_cells * 100,
    Stem_and_Diff   = sum(is_both) / total_cells * 100,
    .groups = "drop"
  )

# 保留成对样本
patients_with_pairs <- df_plot %>%
  dplyr::group_by(patient) %>%
  filter(n() == 2) %>%
  pull(patient) %>%
  unique()

df_plot <- df_plot %>% filter(patient %in% patients_with_pairs)
df_plot$grade <- factor(df_plot$grade, levels = c("older", "younger"))

# ====================================================================
# 4. 统计学核武器：GLMM (二项分布混合模型) 计算真正的 P 值
# ====================================================================
# 确保 grade 是因子，并设置 older 为对照组
meta_ll <- meta %>% mutate(patient = sub("(.*)\\_.*","\\1",sample)) %>%
  filter(patient %in% patients_with_pairs)
meta_ll$grade <- factor(meta_ll$grade, levels = c("older", "younger"))

# 使用 glmer (binomial) 拟合模型 (加入 bobyqa 优化器防报错)
glmer_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

m_stem <- glmer(is_stem_only ~ grade + (1|patient), data = meta_ll, family = binomial, control = glmer_control)
m_diff <- glmer(is_diff_only ~ grade + (1|patient), data = meta_ll, family = binomial, control = glmer_control)
m_both <- glmer(is_both      ~ grade + (1|patient), data = meta_ll, family = binomial, control = glmer_control)

# 提取 P 值
p_stem <- summary(m_stem)$coefficients["gradeyounger", "Pr(>|z|)"]
p_diff <- summary(m_diff)$coefficients["gradeyounger", "Pr(>|z|)"]
p_both <- summary(m_both)$coefficients["gradeyounger", "Pr(>|z|)"]

# 复用之前定义好的格式化函数防零截断
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p == 0) return("p < 2.2e-16")
  if (p < 0.001) return(sprintf("p = %.2e", p))
  sprintf("p = %.3g", p)
}

# ====================================================================
# 5. 绘图 (隐藏 T-test，直接将 GLMM P值写在图顶端)
# ====================================================================
plot_colors <- c("older" = "#E69F00", "younger" = "#0072B2")

p1 <- ggpaired(df_plot, x = "grade", y = "Stem_associated", id = "patient", 
               color = "grade", fill = "white", palette = plot_colors,
               line.color = "black", line.size = 0.4, point.size = 2.5,
               ylab = "percent. of stem-associated (%)", xlab = FALSE) +
  labs(subtitle = sprintf("GLMM %s", fmt_p(p_stem))) + # 添加 GLMM P值
  theme_classic(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(color = "black"), plot.subtitle = element_text(hjust = 0.5, face = "bold"))

p2 <- ggpaired(df_plot, x = "grade", y = "Differentiated", id = "patient", 
               color = "grade", fill = "white", palette = plot_colors,
               line.color = "black", line.size = 0.4, point.size = 2.5,
               ylab = "percent. of differentiated-like (%)", xlab = FALSE) +
  labs(subtitle = sprintf("GLMM %s", fmt_p(p_diff))) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(color = "black"), plot.subtitle = element_text(hjust = 0.5, face = "bold"))

p3 <- ggpaired(df_plot, x = "grade", y = "Stem_and_Diff", id = "patient", 
               color = "grade", fill = "white", palette = plot_colors,
               line.color = "black", line.size = 0.4, point.size = 2.5,
               ylab = "percent. of stem-associated &\ndifferentiated-like (%)", xlab = FALSE) +
  labs(subtitle = sprintf("GLMM %s", fmt_p(p_both))) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(color = "black"), plot.subtitle = element_text(hjust = 0.5, face = "bold"))

final_plot <- p1 | p2 | p3

pdf("figures/Sup_Fig1D_review.pdf", width = 9, height = 4.5)
print(final_plot)
dev.off()

## Figure 1H: tumor cell evolution analysis--------
library(patchwork)
library(dplyr)
library(infercnv)

infercnv_obj <- readRDS("infercnv/all_malignant/run.final.infercnv_obj")
infercnv::plot_cnv(infercnv_obj = infercnv_obj,   #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = paste0('infercnv/all_malignant/infercnv.pdf')) #保存为pdf文件


## Figure Sup_Fig1F: heterogenity of tumor cells, defined by ROGUE--------
library(ROGUE)

expr = as.data.frame(malignant@assays$RNA@counts)
meta <- malignant@meta.data
ent.res <- SE_fun(expr)
SEplot(ent.res)
table(malignant$sample)

rogue.res <- rogue(expr, labels = meta$grade, samples = meta$sample, platform = "UMI", span = 0.6)
rogue.res
rogue.boxplot(rogue.res)

res1 <- rogue.res
res1$sample <- rownames(rogue.res)
res1 <- melt(res1)
res1 <- na.omit(res1)
res1$value <- as.numeric(res1$value)
res1$group <- sub("_.*", "", res1$sample)

# subtype
rogue.res <- rogue(expr, labels = meta$subtype, samples = meta$sample, platform = "UMI", span = 0.6)
rogue.res
rogue.boxplot(rogue.res)

res2 <- rogue.res
res2$sample <- rownames(rogue.res)
res2 <- melt(res2)
res2 <- na.omit(res2)
res2$value <- as.numeric(res2$value)

# the higher the ROGUE score, the higher consistent among cells in this group, so lower score means higher heterogenity
pdf("figures/Sup_Fig1F.pdf", height = 4, width = 4)                                                                                                                                                                                                                                                                                                                                                                             
ggplot(res1, aes(x = variable, y = value, color = variable)) +
  geom_boxplot() +
  geom_point(size = 2) +
  geom_line(aes(group = group), color = "gray") +
  ylab("ROGUE") +
  xlab("Grade") +
  scale_color_manual(values = GradeColor) +
  theme_minimal()
ggplot(res2, aes(x = variable, y = value, color = variable)) +
  geom_boxplot() +
  geom_point(size = 2) +
  ylab("ROGUE") +
  xlab("Subtype") +
  scale_color_manual(values = SubtypeColor) +
  theme_minimal()
dev.off()

write.table(merge@meta.data[,c("sample","grade","cell_type","tumor_ident")], file = "tables/Fig1C.tsv", quote = F, sep = "\t")
write.table(merge@reductions$umap@cell.embeddings, file = "tables/Fig1D.tsv", quote = F, sep = "\t")
write.table(malignant@meta.data[,c("sample","grade","cell_type","tumor_ident","subtype","AUC_MES","AUC_NES")], file = "tables/Fig1E.tsv", quote = F, sep = "\t")
write.table(malignant@reductions$umap@cell.embeddings, file = "tables/Fig1F.tsv", quote = F, sep = "\t")
write.table(gsea_results@result, file = "tables/Fig1G.tsv", quote = F, sep = "\t")
write.table(matrix, file = "tables/Fig1H.tsv", quote = F, sep = "\t")

