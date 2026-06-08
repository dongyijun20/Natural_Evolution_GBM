## Figure 5A: SMF genes regulation overview-------------
library(OmnipathR)

# net <- get_collectri(organism='human', split_complexes=FALSE)
load("tmp/net_omnipath.rdata")
intersection <- read.table("tables/TableS5_refined.tsv", header = 1)[,1]

net_filter <- net %>% 
  dplyr::filter(target %in% intersection) %>%
  distinct(source, target, .keep_all = TRUE)

smf_tf_table <- net_filter %>% 
  group_by(source) %>% 
  summarise(n_targets = n_distinct(target)) %>% 
  filter(n_targets >= 2) %>%
  arrange(desc(n_targets))

pdf("figures/Sup_Fig4A.pdf", height = 4, width = 4)
ggplot(smf_tf_table, aes(x = reorder(source, n_targets), y = n_targets)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Transcription Factor", y = "Number of SMF target genes",
       title = "Upstream TFs of the SMF gene signature") +
  theme_minimal()
dev.off()

interest_tf <- smf_tf_table$source
interest_tf
write.table(interest_tf, "tables/TableS6_review.tsv", quote=F, row.names=F, sep="\t")
                                         
library(igraph)
library(ggraph)
library(tidyverse)

net_plot <- net_filter %>%
  filter(source %in% interest_tf) %>% 
  mutate(reg_type = ifelse(mor == 1, "Activation",
                           ifelse(mor == -1, "Repression", "Unknown")))

g <- graph_from_data_frame(net_plot, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$node_type <- case_when(
  V(g)$name %in% interest_tf ~ "TF",
  V(g)$name %in% intersection ~ "SMF Gene",
    TRUE ~ "Other"
)

set.seed(123)
layout <- create_layout(g, layout = "kk")

pdf("figures/Fig5A_review.pdf", height = 5, width = 6)
ggraph(layout) +
  geom_edge_arc(aes(color = reg_type),
  strength = 0.1,
  arrow = arrow(length = unit(2, 'mm'), type = "closed"),
  end_cap = circle(3, 'mm'),
  width = 0.5, alpha = 0.5) +
  # Nodes: Size by degree (Hubs are bigger)
  geom_node_point(aes(color = node_type, size = degree)) +
  # Labels: Only label TFs and SMF genes
  geom_node_text(aes(label = name),
  repel = TRUE,
  size = 3.5) +
  # Scales & Colors
  scale_edge_color_manual(values = c("Activation" = "forestgreen", "Repression" = "firebrick")) +
  scale_color_manual(values = c("SMF Gene" = "orange",
  "TF" = "steelblue",
  "Other" = "grey30")) +
  # Fix node size range
  scale_size_continuous(range = c(2, 8)) +
  theme_void() +
  labs(title = "Core TF-Target Regulatory Network",
  edge_color = "Regulation",
  color = "Node Type",
  size = "Connectivity")
dev.off()


## Figure 5B: decoupleR of interested TFs, to ensure these TF have higher activiti in S2 BMDM------------
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(decoupleR)
library(ggplot2)
library(patchwork)
library(OmnipathR)

options(omnipathr.curl_verbose = TRUE)
omnipath_set_cachedir(tempdir())
options(future.globals.maxSize = 16 * 1024^3)

sce <- readRDS("result/BMDM_withXHsample_harmony.rds")

Idents(sce) <- "subcluster"

mat <- as.matrix(sce@assays$RNA@data)

plan("multisession",workers = 8)
acts <- run_ulm(mat, net, minsize = 5,
                .source='source', .target='target',.mor='mor')

sce[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = sce) <- "tfsulm"
sce <- ScaleData(sce)
sce@assays$tfsulm@data <- sce@assays$tfsulm@scale.data

df <- t(as.matrix(sce@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(sce)) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score), .groups = 'drop')

tfs_highest_in_S2 <- df %>%
  group_by(source) %>%
  slice_max(order_by = mean, n = 1) %>% 
  filter(cluster == "S2") %>%  
  pull(source)                  

interest_tf_final <- intersect(interest_tf, tfs_highest_in_S2)

net_plot_final <- net_plot %>% 
  filter(source %in% interest_tf_final)

g <- graph_from_data_frame(net_plot_final, directed = TRUE)

# 计算度 (Degree) 用来控制节点大小
V(g)$degree <- degree(g, mode = "all")

# 设置节点类型 (Node Type)
V(g)$node_type <- case_when(
  V(g)$name %in% interest_tf_final ~ "TF",
  V(g)$name %in% intersection ~ "SMF Gene",
  TRUE ~ "Other"
)

# 设置布局
set.seed(123)
layout <- create_layout(g, layout = "nicely") # "kk" 或 "nicely" 都不错

# 绘图输出
# pdf("figures/Fig5A_review.pdf", height = 6, width = 7) # 如果需要保存PDF请取消注释

ggraph(layout) +
  # 1. 画边 (Edges)
  geom_edge_arc(aes(color = reg_type),
                strength = 0.1, # 弧度大小
                arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                end_cap = circle(3, 'mm'), # 连线不直接插进点里，留点空隙
                width = 0.8, alpha = 0.6) +
  
  # 2. 画节点 (Nodes)
  geom_node_point(aes(color = node_type, size = degree)) +
  
  # 3. 画标签 (Labels)
  geom_node_text(aes(label = name),
                 repel = TRUE,
                 fontface = "bold", # 加粗字体
                 size = 3.5) +
  
  # 4. 颜色与标尺设置
  # 因为没有 Act/Rep，我们把 Association 设为深灰色
  scale_edge_color_manual(values = c("Association" = "grey50", 
                                     "Activation" = "forestgreen", 
                                     "Repression" = "firebrick")) +
  
  scale_color_manual(values = c("Intersection Gene" = "#E69F00", # 橙色
                                "TF" = "#56B4E9",              # 蓝色
                                "Other" = "grey30")) +
  
  scale_size_continuous(range = c(3, 10)) + # 调整点的大小范围
  
  # 5. 主题设置
  theme_void() +
  labs(title = "TF-Target Regulatory Network",
       subtitle = "Based on ChEA/ENCODE ChIP-seq data",
       edge_color = "Regulation",
       color = "Node Type",
       size = "Connectivity")

top_acts_mat <- df %>%
  filter(source %in% interest_tf_final) %>% # 条件2: 必须在 S2 中活性最高
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()


# 检查一下是否还有剩下的 TF，防止绘图报错
if(ncol(top_acts_mat) == 0) {
  stop("没有 TF 同时满足 '在 interest_tf 列表里' 且 '在 S2 中活性最高' 这两个条件！")
}

palette_length = 100
my_color = colorRampPalette(c("darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# 确保行序正确
top_acts_mat <- top_acts_mat[paste0("S", 0:4), , drop = FALSE] # drop=FALSE 防止只剩1列变成向量


palette_length = 100
my_color = colorRampPalette(c("darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

top_acts_mat <- top_acts_mat[paste0("S",0:4),]

library(pheatmap)
pdf("figures/Fig5B_review.pdf", height = 4, width = 6)
pheatmap(top_acts_mat,
         border_color = NA,
         color = my_color,
         breaks = my_breaks,
         cluster_rows = FALSE,
         cluster_cols = TRUE,  
         fontsize_row = 12,
         fontsize_col = 10,
         cellwidth = 15, 
         cellheight = 25,
         main = "TF Activity (Highest in S2)")
dev.off()

## Figure 5C: bulk ATAC-seq validation--------
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# TSS ±2kb
tss_gr <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), upstream = 2000, downstream = 2000)

library(Rsubread)

tss_saf <- data.frame(
  GeneID = names(tss_gr),
  Chr = as.character(seqnames(tss_gr)),
  Start = start(tss_gr),
  End = end(tss_gr),
  Strand = ifelse(strand(tss_gr) == "*", "+", as.character(strand(tss_gr)))  # replace '*' with '+'
)
tss_saf <- tss_saf %>% filter(Start >= 1)

# Save with tab-separated format and no row/col names
write.table(tss_saf, "tmp/tss_regions.saf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

bam_files <- c("data/bulk_ATAC/Cleandata/Erastin/Erastin.bam", "data/bulk_ATAC/Cleandata/NC/NC.bam")

counts <- featureCounts(files = bam_files,
                        annot.ext = "tmp/tss_regions.saf",
                        isGTFAnnotationFile = FALSE,
                        isPairedEnd = TRUE,
                        nthreads = 8)

saveRDS(counts, file = "result/bulk_ATAC_count_matrix.rds")
# counts <- readRDS("result/bulk_ATAC_count_matrix.rds")
count_matrix <- counts$counts
colnames(count_matrix) <- c("Erastin", "NC")

library(edgeR)

group <- factor(c("Erastin", "NC"))
dge <- DGEList(counts = count_matrix, group = group)

dge <- estimateCommonDisp(dge)
et <- exactTest(dge, dispersion = 0.1)

res_atac <- topTags(et, n = Inf)$table %>%
  tibble::rownames_to_column("region") %>%
  dplyr::select(region, logFC) %>%
  rename(logFC_accessibility = logFC)

library(org.Hs.eg.db)
symbol_map <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = res_atac$region,
                                    columns = c("SYMBOL","ENSEMBL"),
                                    keytype = "ENTREZID")
merged <- merge(res_atac, symbol_map, by.x = "region", by.y = "ENTREZID")

plot_df <- merged %>%
  mutate(
    group = case_when(
      SYMBOL %in% intersection ~ "SMF gene",
      SYMBOL %in% interest_tf ~ "TF",
      TRUE ~ "Background"
    )
  )

rna_logFC <- res %>%
  select(gene_name, logFC) %>%
  distinct() %>%
  filter(!is.na(gene_name)) %>%
  rename(SYMBOL = gene_name, logFC_expr = logFC)

plot_df <- plot_df %>%
  left_join(rna_logFC, by = "SYMBOL") %>%
  mutate(logFC_expr = ifelse(is.na(logFC_expr),
                             rna_logFC$logFC_expr[match(plot_df$ENSEMBL, rna_logFC$SYMBOL)],
                             logFC_expr))

color_map <- c(
  "SMF gene" = "orange",
  "TF" = "steelblue",
  "Background" = "gray80"
)

library(ggrepel)

pdf("figures/Sup_Fig4C.pdf", height = 5, width = 6)
ggplot(plot_df, aes(x = logFC_expr, y = logFC_accessibility)) +
  geom_point(aes(color = group), alpha = 0.3, size = 0.8) +
  
  geom_point(
    data = plot_df %>% filter(group != "Background"),
    aes(color = group),
    size = 3
  )  +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  
  geom_text_repel(
    data = plot_df %>%
      filter(SYMBOL %in% c(intersection, interest_tf)),
    aes(label = SYMBOL, color = group),
    size = 3, max.overlaps = 500
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_map) +
  theme_minimal() +
  labs(
    title = "Erastin treated THP-1 RNA & ATAC changes",
    x = "log2 FC (RNA-seq)",
    y = "log2 FC (ATAC-seq)",
    color = "Gene Set"
  )
dev.off()

## Figure 5D: fingerprint of ATAC-seq-------
library(ATACseqQC)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

genome <- Hsapiens

tf_motif <- query(MotifDb, "CEBPA")
pfm <- as.list(tf_motif)[[1]]

bamfile_A <- "data/bulk_ATAC/Cleandata/Erastin/Erastin.bam"
result_A <- factorFootprints(bamfile_A, pfm=pfm, genome=genome,
                             min.score="95%", upstream=100, downstream=100)

bamfile_B <- "data/bulk_ATAC/Cleandata/NC/NC.bam"
result_B <- factorFootprints(bamfile_B, pfm=pfm, genome=genome,
                             min.score="95%", upstream=100, downstream=100)

avg_A <- (colMeans(result_A$signal$`+`) + colMeans(result_A$signal$`-`)) / 2
avg_B <- (colMeans(result_B$signal$`+`) + colMeans(result_B$signal$`-`)) / 2

df <- data.frame(
  Position = -((ncol(result_A$signal$`+`) - 1) / 2):((ncol(result_A$signal$`+`) - 1) / 2),
  Erastin = avg_A,
  NC = avg_B
)
df_long1 <- pivot_longer(df, cols = -Position, names_to = "Group", values_to = "Signal")

pdf("figures/Sup_Fig4E.pdf", height = 4, width = 5)
ggplot(df_long1, aes(x = Position, y = Signal, color = Group)) +
  geom_line(size = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TF Footprint (CEBPA)",
    x = "Position relative to motif (bp)",
    y = "Insertion Signal"
  ) +
  scale_color_manual(values = c("Erastin" = "#D95F02", "NC" = "#1B9E77"))
dev.off()

## Figure 5E: UMAP of scATAC-seq data-------------
library(Signac)
library(ArchR)

integrated <- readRDS("result/XH02_signac_integration.rds")

integrated$group <- paste(integrated$grade, integrated$celltype)

group_pal <- c("older Malignant" = "#8C1217", "older Myeloid" = "#171C40", "older Oligodendrocyte" = "#14602E", "older Stromal" = "#5E1C62", "older T Cell" = "#A6531D",
               "younger Malignant" = "#F68A8E", "younger Myeloid" = "#8087B5", "younger Oligodendrocyte" = "#6FC98B", "younger Stromal" = "#C280C9", "younger T Cell" = "#FBBF94")

pdf("figures/Fig4D.pdf", width = 6, height = 4)
DimPlot(
  object = integrated,
  group.by = 'group', alpha = 0.8,
  cols = group_pal, shuffle = T,
  repel = TRUE) + 
  ggtitle('scATAC-seq')
dev.off()

DefaultAssay(integrated) <- "RNA"
library(UCell)
library(readxl)
score_list <- readRDS("result/gene_list_Fig3A.rds")
markers<-as.list(read_excel("reference/mmc2.xlsx", skip = 3))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1/S","G2/M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NES = c('S100A10','FOSL2','SPP1','CAV1','ANXA1','VIM','CD44','SERPINH1',
                      'LGALS3','CEBPB','ATF5','LGALS1'),
              Hypoxia = score_list$Hypoxia,
              BMDM = score_list$BMDM,
              Ferroptosis_neg = score_list$Ferroptosis_neg)
integrated <- AddModuleScore_UCell(integrated, markers, name = NULL)

pdf("figures/Sup_Fig4D.pdf", height = 4, width = 3)
VlnPlot(integrated[,integrated$celltype=="Malignant"], group.by = "grade", c("MES","NES"), stack = T, flip = T) + NoLegend() 
VlnPlot(integrated[,integrated$celltype=="Myeloid"], group.by = "grade", c("BMDM","Hypoxia","Ferroptosis_neg"), stack = T, flip = T) + NoLegend()
dev.off()

## Figure 5F: Motif Deviation Score of scATAC---------------
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DefaultAssay(integrated) <- "peaks"
integrated <- AddMotifs(
  object = integrated, assay = "peaks",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = integrated,
  ident.1 = 'older Myeloid',
  ident.2 = 'younger Myeloid',
  group.by = "group",
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2 & da_peaks$avg_log2FC > 1, ])
length(top.da.peak)

enriched.motifs <- FindMotifs(
  object = integrated,
  features = top.da.peak
)

MotifPlot(
  object = integrated,
  motifs = head(rownames(enriched.motifs))
)

library(GenomeInfoDb)
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
counts <- GetAssayData(integrated[["peaks"]], slot = "counts")
current_ranges <- granges(integrated[["peaks"]])
names(current_ranges) <- rownames(counts)
filtered_ranges <- current_ranges[seqnames(current_ranges) %in% standard_chromosomes]
counts_filtered <- counts[names(filtered_ranges), ]
new_peaks_assay <- CreateChromatinAssay(counts = counts_filtered, ranges = filtered_ranges)
integrated[["peaks"]] <- new_peaks_assay

integrated <- RunChromVAR(
  object = integrated,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(integrated, file = "result/XH02_signac_integration.rds")

DefaultAssay(integrated) <- 'chromvar'

differential.activity <- FindMarkers(
  object = integrated, 
  group.by = "group",
  ident.1 = 'older Myeloid',
  ident.2 = 'younger Myeloid',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

top.diff.peak <- rownames(differential.activity[differential.activity$p_val < 0.005 & differential.activity$pct.1 > 0.2 & 
                                                  differential.activity$avg_diff > 1, ])

MotifPlot(
  object = integrated,
  motifs = head(top.diff.peak),
  assay = 'peaks'
)

opts <- list(species = 9606)
pfm_list <- getMatrixSet(JASPAR2020, opts)

motif_df <- data.frame(
  motif_id = sapply(pfm_list, ID),
  tf_name = sapply(pfm_list, name),
  stringsAsFactors = FALSE
)

matched_tfs <- motif_df[motif_df$motif_id %in% top.diff.peak, ]

DefaultAssay(integrated) <- "chromvar"
chromvar_avg <- AverageExpression(integrated, assays = "chromvar", group.by = "group")$chromvar

gene_avg <- AverageExpression(integrated, assays = "RNA", group.by = "group")$RNA

matched_tfs$clean_tf <- gsub("::.*", "", matched_tfs$tf_name)
matched_tfs$clean_tf <- gsub("\\(.*\\)", "", matched_tfs$clean_tf)
matched_tfs$clean_tf <- gsub("var\\.\\d+", "", matched_tfs$clean_tf)
matched_tfs <- matched_tfs %>% filter(clean_tf %in% rownames(gene_avg) & tf_name %in% enriched.motifs$motif.name)

motif_matrix <- chromvar_avg[matched_tfs$motif_id, ]
gene_matrix <- gene_avg[matched_tfs$clean_tf, ]
rownames(gene_matrix) <- matched_tfs$tf_name

motif_matrix <- t(scale(t(motif_matrix)))
gene_matrix <- t(scale(t(gene_matrix)))

keep_motif <- which(apply(motif_matrix, 1, function(x) which.max(x) == 2))
keep_gene <- which(apply(gene_matrix, 1, function(x) x[2]>x[7]))
keep_row <- intersect(keep_motif, keep_gene)

motif_mat <- as.matrix(motif_matrix)[keep_row,]
gene_mat <- as.matrix(gene_matrix)[keep_row,]

dist_mat <- dist(gene_mat)
hc <- hclust(dist_mat)
ordered_TFs <- hc$order

motif_mat <- motif_mat[ordered_TFs,]
gene_mat <- gene_mat[ordered_TFs,]

SMF_tss <- annotations[annotations$gene_name %in% intersection]
promoter_regions <- promoters(SMF_tss, upstream = 2000, downstream = 500)

motif_object <- Motifs(integrated)
motif_matrix <- motif_object@data[,rownames(motif_mat)]

peak_gr <- granges(integrated[["peaks"]])
names(peak_gr) <- rownames(integrated[["peaks"]])
peaks_in_promoter <- subsetByOverlaps(peak_gr, SMF_promoter_region)

motif_promoter_matrix <- motif_matrix[rownames(motif_matrix) %in% names(peaks_in_promoter), ]

motif_sums <- Matrix::colSums(motif_promoter_matrix)
motif_keep <- as.vector(which(motif_sums > 0))

motif_mat_filtered <- motif_mat[motif_keep,]
gene_mat_filtered <- gene_mat[motif_keep,]

library(Matrix)
motif_avg_score <- motif_promoter_matrix[, motif_keep]
motif_score_df <- as.data.frame(as.matrix(motif_avg_score))
motif_score_df$peak <- rownames(motif_promoter_matrix)

peak_to_gene <- findOverlaps(peaks_in_promoter, SMF_promoter_region)
peak_gene_map <- data.frame(
  peak = names(peaks_in_promoter)[queryHits(peak_to_gene)],
  gene = SMF_promoter_region$gene_name[subjectHits(peak_to_gene)]
)

motif_score_long <- reshape2::melt(motif_score_df, id.vars = "peak", variable.name = "TF", value.name = "score")
motif_score_long <- merge(motif_score_long, peak_gene_map, by = "peak") %>% 
  distinct()
motif_score_summarized <- motif_score_long %>%
  dplyr::group_by(TF, gene) %>%
  dplyr::summarise(score = sum(score), .groups = 'drop')

motif_df <- melt(motif_mat_filtered)
colnames(motif_df) <- c("TF", "CellType", "Activity")
motif_df$Source <- "Motif"

gene_df <- melt(gene_mat_filtered)
colnames(gene_df) <- c("TF", "CellType", "Activity")
gene_df$Source <- "Gene"

library(RColorBrewer)
library(scales)

colors <- rev(brewer.pal(11, "RdYlBu"))
activity_range <- range(motif_df$Activity, na.rm = TRUE)
mid_val <- rescale(0, to = c(0, 1), from = activity_range)

p1 <- ggplot(motif_df, aes(x = CellType, y = TF, fill = Activity)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colors,
    values = rescale(c(activity_range[1], 0, activity_range[2])),
    limits = activity_range
  ) +
  theme_minimal() +
  labs(title = "chromVAR Motif Activity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(gene_df, aes(x = CellType, y = TF, fill = Activity)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  labs(title = "TF Gene Activity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(motif_score_summarized, aes(x = gene, y = TF, fill = score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  theme_minimal() +
  labs(title = "Motif binding to 14 genes' promoters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("figures/Fig4E.pdf", height = 16, width = 12)
(p1 | p2 | p3) + plot_layout(ncol = 3, widths = c(1, 1, 1))
dev.off()

interest_tf2 <- intersect(gene_df$TF, interest_tf)
interest_tf2

unique(as.vector(motif_df$TF[gene_df$TF=="CEBPA"]))
motif_score_summarized %>% filter(TF == "MA0102.4")
motif_df %>% filter(TF == "MA0102.4")

# motif footprinting-------
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# For sample 1
cells_1 <- grep("_1$", colnames(integrated), value = TRUE)
barcode_map_1 <- setNames(
  object = sub("_1$", "", cells_1),
  nm = cells_1
)

frags1 <- CreateFragmentObject(
  path = "data/scATAC_new/result/fragments.tsv.gz",
  cells = barcode_map_1
)

# For sample 2
cells_2 <- grep("_2$", colnames(integrated), value = TRUE)
barcode_map_2 <- setNames(
  object = sub("_2$", "", cells_2),
  nm = cells_2
)

frags2 <- CreateFragmentObject(
  path = "data/scATAC_new/result2/fragments.tsv.gz",
  cells = barcode_map_2
)

# Assign the fragment objects to the Seurat object
Fragments(integrated) <- list(frags1, frags2)

DefaultAssay(integrated) <- "peaks"

integrated <- Footprint(
  object = integrated,
  motif.name = interest_tf2,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

plot.data <- GetFootprintData(
  object = integrated,
  features = interest_tf2,
  group.by = "group"
) %>% na.omit() %>% 
  filter(!group %in% c("younger T Cell","older T Cell","younger Stromal","older Stromal")) %>% 
  mutate(group = factor(group, levels = rev(names(group_pal))))
  
plot.list <- lapply(interest_tf2, function(x){
  ggplot(plot.data[plot.data$feature==x,], aes(x = position, y = norm.value, color = group)) +
    geom_line(alpha = 0.8) +
    theme_classic() +
    scale_color_manual(values = group_pal) +
    labs(y = "Footprint signal", color = "Cell type") +
    ggtitle(x)
})

pdf("figures/Fig4F.pdf", height = 4, width = 6)
plot.list
dev.off()

pdf("figures/Fig4G.pdf", height = 5, width = 7)
MotifPlot(
  object = integrated,
  motifs = unique(as.vector(motif_df$TF[gene_df$TF %in% interest_tf2])),
  assay = 'peaks'
)
dev.off()

## Figure 5G: regulatory function of intersection genes------------
library(cicero)
library(monocle3)
library(SeuratWrappers)

DefaultAssay(integrated) <- "peaks"
# convert to CellDataSet format and make the cicero object
integrated.cds <- as.cell_data_set(x = integrated)
integrated.cicero <- make_cicero_cds(integrated.cds, reduced_coordinates = reducedDims(integrated.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(integrated.cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(integrated) <- links

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")

Annotation(integrated) <- annotations

nupr1_tss <- annotations[annotations$gene_name == "NUPR1"]
nupr1_tss <- resize(nupr1_tss, width = 1, fix = "start")
nupr1_promoter_region <- resize(nupr1_tss, width = 4000, fix = "center")
peaks <- granges(integrated[["peaks"]])
promoter_peaks_nupr1 <- subsetByOverlaps(peaks, nupr1_promoter_region)

link_start <- GRanges(
  seqnames = seqnames(links),
  ranges = IRanges(start = start(links), width = 1)
)
link_end <- GRanges(
  seqnames = seqnames(links),
  ranges = IRanges(start = end(links), width = 1)
)

start_hits <- findOverlaps(link_start, promoter_peaks)
end_hits <- findOverlaps(link_end, promoter_peaks)

linked_idx <- unique(c(queryHits(start_hits), queryHits(end_hits)))
nupr1_links <- links[linked_idx]
nupr1_links <- nupr1_links[nupr1_links$score > 0.1]
nupr1_links
linked_idx_filtered <- match(nupr1_links, links)

enhancer_peaks <- GRanges()
for (i in linked_idx_filtered) {
  this_link <- links[i]
  
  s <- GRanges(seqnames = seqnames(this_link),
               ranges = IRanges(start(this_link), width = 1))
  e <- GRanges(seqnames = seqnames(this_link),
               ranges = IRanges(end(this_link), width = 1))
  
  if (countOverlaps(s, promoter_peaks) > 0) {
    enhancer_peaks <- c(enhancer_peaks, e)
  } else {
    enhancer_peaks <- c(enhancer_peaks, s)
  }
}

enhancer_in_peakset_idx <- findOverlaps(enhancer_peaks, peaks)
enhancer_peak_names <- names(peaks)[subjectHits(enhancer_in_peakset_idx)]

peak_matrix <- GetAssayData(integrated, assay = "peaks", slot = "counts")
enhancer_matrix <- peak_matrix[enhancer_peak_names, ]

older_cells <- colnames(integrated)[integrated$group == "older Myeloid"]
younger_cells <- colnames(integrated)[integrated$group == "younger Myeloid"]

older_avg <- Matrix::rowMeans(enhancer_matrix[, older_cells])
younger_avg <- Matrix::rowMeans(enhancer_matrix[, younger_cells])

df <- data.frame(
  peak = enhancer_peak_names,
  older = older_avg,
  younger = younger_avg,
  log2FC = log2(older_avg + 1e-5) - log2(younger_avg + 1e-5)
)

df_up_in_older <- df[order(-df$log2FC), ]
df_up_in_older

motif.enrichment <- FindMotifs(
  object = integrated,
  features = enhancer_peak_names
)

top_motifs <- motif.enrichment[motif.enrichment$p.adjust < 0.05 & motif.enrichment$fold.enrichment > 1, ]
intersect(interest_tf2, top_motifs$motif.name)

unique(as.vector(motif_df$TF[gene_df$TF=="ETV4"]))

motif_matrix <- GetMotifData(integrated)
etv4_peaks <- rownames(motif_matrix)[which(motif_matrix[,"MA0764.2"] == 1)]
cebpa_peaks <- rownames(motif_matrix)[which(motif_matrix[,"MA0102.4"] == 1)]

etv4_enhancer_peaks <- intersect(etv4_peaks, enhancer_peak_names)
etv4_enhancer_gr <- granges(integrated[["peaks"]])[etv4_enhancer_peaks][1:4]
etv4_enhancer_gr$color <- "forestgreen"

promoter_peaks$color <- "red"
highlight_region <- c(etv4_enhancer_gr, promoter_peaks_nupr1)
selected_links <- nupr1_links[c(5,6,7,3)]

coverage_region <- range(selected_links)
start(coverage_region) <- start(coverage_region)-2000
end(coverage_region) <- end(coverage_region)+2000

Links(integrated) <- links

integrated_subset <- integrated[,!integrated$group %in% c("younger T Cell","older T Cell","younger Stromal","older Stromal")]
cov_plot <- CoveragePlot(
  object = integrated_subset,
  region = coverage_region,
  annotation = FALSE,
  peaks = FALSE,
  region.highlight = highlight_region,
  assay.scale = "common",
  group.by = "group"
) + scale_fill_manual(values = group_pal)

gene_plot <- AnnotationPlot(
  object = integrated_subset,
  region = coverage_region
)

peak_plot <- PeakPlot(
  object = integrated_subset,
  region = coverage_region
)

Links(integrated_subset) <- links
link_plot <- LinkPlot(
  object = integrated_subset,
  region = coverage_region
)

expr_plot_nupr1 <- ExpressionPlot(
  object = integrated_subset,
  features = "NUPR1",
  assay = "RNA",
  group.by = "group"
) + scale_fill_manual(values = group_pal)

expr_plot_etv4 <- ExpressionPlot(
  object = integrated_subset,
  features = "ETV4",
  assay = "RNA",
  group.by = "group"
) + scale_fill_manual(values = group_pal)

expr_plot_cebpa <- ExpressionPlot(
  object = integrated_subset,
  features = "CEBPA",
  assay = "RNA",
  group.by = "group"
) + scale_fill_manual(values = group_pal)

combined_expr_plot <- expr_plot_nupr1 + expr_plot_etv4 + expr_plot_cebpa + plot_layout(ncol = 3)

pdf("figures/Fig4H.pdf", height = 6, width = 10)
CombineTracks(
  plotlist = list(cov_plot, peak_plot, gene_plot, link_plot),
  expression.plot = combined_expr_plot,
  heights = c(6, 1, 1, 1),
  widths = c(8, 3)
)
dev.off()

## Figure 5: CEBPA Cut&Tag tracks (two genes from intersection) ---------
if (file.exists("data/CnT/NC1.bw")) {
  library(rtracklayer)
  cnT_genes <- c("GPX4", "NUPR1")
  bw_path <- "data/CnT/NC1.bw"
  peak_path <- "data/CnT/NC1_peaks.narrowPeak"
  cebpa_peaks_gr <- rtracklayer::import(peak_path)
  cebpa_peaks_gr <- keepStandardChromosomes(cebpa_peaks_gr, pruning.mode = "coarse")
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(cebpa_peaks_gr) <- "UCSC")
  region_width <- 8000L

  make_cnt_panel <- function(gene_name) {
    gene_anno <- annotations[annotations$gene_name == gene_name]
    if (length(gene_anno) == 0) return(NULL)
    tss <- resize(gene_anno, width = 1L, fix = "start")
    region <- resize(tss, width = region_width, fix = "center")
    region <- region[1]
    peaks_in_region <- subsetByOverlaps(cebpa_peaks_gr, region)
    reg_start <- start(region)
    reg_end <- end(region)

    peak_track <- if (length(peaks_in_region) > 0) {
      peak_df <- data.frame(
        xmin = start(peaks_in_region),
        xmax = end(peaks_in_region),
        ymin = 0, ymax = 1
      )
      ggplot(peak_df) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "#2166ac", color = NA) +
        scale_x_continuous(limits = c(reg_start, reg_end), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_classic(base_size = 10) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
        ) +
        labs(y = "CEBPA\npeaks", title = gene_name)
    } else {
      ggplot() +
        scale_x_continuous(limits = c(reg_start, reg_end), expand = c(0, 0)) +
        theme_classic(base_size = 10) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
        ) +
        labs(y = "CEBPA\npeaks", title = gene_name)
    }

    region_bw <- region
    bw_seqinfo <- rtracklayer::seqinfo(rtracklayer::BigWigFile(bw_path))
    if (!as.character(seqnames(region)) %in% seqnames(bw_seqinfo)) {
      suppressWarnings(GenomeInfoDb::seqlevelsStyle(region_bw) <- "NCBI")
    }
    bigwig_track <- Signac::BigwigTrack(
      region = region_bw,
      bigwig = list(CEBPA = bw_path),
      type = "coverage",
      smooth = 200,
      y_label = "CEBPA Cut&Tag",
      bigwig.scale = "common",
      ymax = "q95"
    ) + theme_classic(base_size = 10)

    gene_track <- Signac::AnnotationPlot(
      object = integrated_subset,
      region = region
    )

    Signac::CombineTracks(
      plotlist = list(peak_track, bigwig_track, gene_track),
      heights = c(0.8, 2, 1)
    )
  }

  p_gpx4 <- make_cnt_panel(cnT_genes[1])
  p_nupr1 <- make_cnt_panel(cnT_genes[2])
  if (!is.null(p_gpx4) && !is.null(p_nupr1)) {
    p_combined <- (p_gpx4 | p_nupr1) +
      plot_annotation(
        title = "CEBPA Cut&Tag at SMF gene loci",
        theme = theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        )
      )
    pdf("figures/Fig5_CnT.pdf", height = 5, width = 11)
    print(p_combined)
    dev.off()
  }
} else {
  message("data/CnT/NC1.bw not found, skipping Cut&Tag figure.")
}

## Figure 5H review: bulk ATAC accessibility around CEBPA Cut&Tag peaks---------
if (Sys.which("bamCoverage") == "" || Sys.which("computeMatrix") == "") {
  warning("deepTools not found; skipping Fig5H_review and Sup_Fig5H_review.")
} else {
  peak_path <- "data/CnT/NC1_peaks.narrowPeak"
  peak_bed <- "tmp/NC1_top20k_peaks.bed"
  erastin_bam <- "data/bulk_ATAC/Cleandata/Erastin/Erastin.bam"
  nc_bam <- "data/bulk_ATAC/Cleandata/NC/NC.bam"
  erastin_bw <- "tmp/Erastin.RPKM.bw"
  nc_bw <- "tmp/NC.RPKM.bw"
  matrix_path <- "tmp/CnT_peak_bulkATAC_profile.matrix.gz"
  
  up_bp <- 1500
  down_bp <- 1500
  bin_size <- 10
  center_window <- 250
  n_top_peaks <- 20000
  
  erastin_col <- "#D95F02"
  nc_col <- "#1B9E77"
  delta_col <- "#7B3294"
  
  dir.create("tmp", showWarnings = FALSE)
  dir.create("figures", showWarnings = FALSE)
  
  peak_tbl <- read.delim(peak_path, header = FALSE, stringsAsFactors = FALSE) %>%
    transmute(
      chr = ifelse(grepl("^chr", V1), V1, paste0("chr", V1)),
      start = V2,
      end = V3,
      score = V7
    ) %>%
    arrange(desc(score)) %>%
    slice_head(n = n_top_peaks)
  
  write.table(
    peak_tbl[, c("chr", "start", "end")],
    file = peak_bed,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if (!file.exists(erastin_bw)) {
    cmd <- sprintf(
      "bamCoverage -b %s -o %s --normalizeUsing RPKM --binSize 25 -p 8",
      shQuote(erastin_bam),
      shQuote(erastin_bw)
    )
    stopifnot(system(cmd) == 0)
  }
  
  if (!file.exists(nc_bw)) {
    cmd <- sprintf(
      "bamCoverage -b %s -o %s --normalizeUsing RPKM --binSize 25 -p 8",
      shQuote(nc_bam),
      shQuote(nc_bw)
    )
    stopifnot(system(cmd) == 0)
  }
  
  cmd <- sprintf(
    paste(
      "computeMatrix reference-point --referencePoint center",
      "-R %s -S %s %s -a %d -b %d --skipZeros -p 8 -o %s"
    ),
    shQuote(peak_bed),
    shQuote(erastin_bw),
    shQuote(nc_bw),
    up_bp,
    down_bp,
    shQuote(matrix_path)
  )
  stopifnot(system(cmd) == 0)
  
  mat_df <- read.table(
    gzfile(matrix_path),
    sep = "\t",
    header = FALSE,
    skip = 1,
    check.names = FALSE
  )
  
  signal_mat <- as.matrix(mat_df[, -(1:6), drop = FALSE])
  n_bins <- (up_bp + down_bp) / bin_size
  erastin_mat <- signal_mat[, seq_len(n_bins), drop = FALSE]
  nc_mat <- signal_mat[, n_bins + seq_len(n_bins), drop = FALSE]
  
  x <- seq(
    from = -up_bp + bin_size / 2,
    by = bin_size,
    length.out = n_bins
  )
  
  overlay_df <- bind_rows(
    data.frame(position = x, signal = colMeans(erastin_mat), group = "Erastin"),
    data.frame(position = x, signal = colMeans(nc_mat), group = "NC")
  )
  
  delta_mat <- erastin_mat - nc_mat
  mean_delta <- colMeans(delta_mat)
  
  bootstrap_ci <- function(xmat, n_boot = 80, seed = 123) {
    set.seed(seed)
    boot_mean <- vapply(
      seq_len(n_boot),
      function(i) {
        idx <- sample(seq_len(nrow(xmat)), replace = TRUE)
        colMeans(xmat[idx, , drop = FALSE])
      },
      numeric(ncol(xmat))
    )
    data.frame(
      lower = apply(boot_mean, 1, quantile, probs = 0.025),
      upper = apply(boot_mean, 1, quantile, probs = 0.975)
    )
  }
  
  ci_df <- bootstrap_ci(delta_mat)
  delta_df <- data.frame(
    position = x,
    delta = mean_delta,
    lower = ci_df$lower,
    upper = ci_df$upper
  )
  
  center_idx <- abs(x) <= center_window
  erastin_center <- rowMeans(erastin_mat[, center_idx, drop = FALSE])
  nc_center <- rowMeans(nc_mat[, center_idx, drop = FALSE])
  median_log2fc <- median(log2((erastin_center + 1e-6) / (nc_center + 1e-6)))
  wilcox_p <- wilcox.test(erastin_center, nc_center, paired = TRUE, alternative = "greater")$p.value
  
  p_delta <- ggplot(delta_df, aes(x = position, y = delta)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = delta_col, alpha = 0.16) +
    geom_area(data = subset(delta_df, delta >= 0), fill = erastin_col, alpha = 0.24) +
    geom_area(data = subset(delta_df, delta < 0), fill = nc_col, alpha = 0.24) +
    geom_line(color = delta_col, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    annotate(
      "text",
      x = min(delta_df$position),
      y = max(delta_df$upper, na.rm = TRUE),
      hjust = 0, vjust = 1.1, size = 3.1,
      label = paste0(
        "Erastin - NC\n",
        "median log2FC (center ±", center_window, " bp) = ", sprintf("%.3f", median_log2fc), "\n",
        "Wilcoxon p = ", format(wilcox_p, scientific = TRUE, digits = 2)
      )
    ) +
    theme_classic(base_size = 11) +
    labs(
      title = "Bulk ATAC accessibility shift at CEBPA Cut&Tag peaks",
      x = "Distance to CEBPA peak center (bp)",
      y = "Delta ATAC signal"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p_overlay <- ggplot(overlay_df, aes(x = position, y = signal, color = group)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Erastin" = erastin_col, "NC" = nc_col)) +
    theme_classic(base_size = 11) +
    labs(
      title = "Bulk ATAC metaprofile around CEBPA Cut&Tag peaks",
      x = "Distance to CEBPA peak center (bp)",
      y = "ATAC signal (RPKM)",
      color = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.82, 0.9)
    )
  
  pdf("figures/Fig5H_review.pdf", height = 3.4, width = 4.0)
  print(p_delta)
  dev.off()
  
  pdf("figures/Sup_Fig5H_review.pdf", height = 3.4, width = 4.0)
  print(p_overlay)
  dev.off()
}

write.table(net_plot, file = "tables/Fig5A.tsv", quote = F, sep = "\t")
write.table(top_acts_mat, file = "tables/Fig5B.tsv", quote = F, sep = "\t")
write.table(integrated@meta.data[,c("grade","celltype")], file = "tables/Fig5C.tsv", quote = F, sep = "\t")
write.table(integrated@reductions$umap@cell.embeddings, file = "tables/Fig5D.tsv", quote = F, sep = "\t")
write.table(df_long1, file = "tables/Fig5E.tsv", quote = F, sep = "\t")
write.table(df_long2, file = "tables/Fig5F.tsv", quote = F, sep = "\t")



