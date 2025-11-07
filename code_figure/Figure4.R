## Figure 4A: SSS genes regulation overview-------------
library(OmnipathR)

# net <- get_collectri(organism='human', split_complexes=FALSE)
load("tmp/net_omnipath.rdata")
intersection <- read.table("tables/TableS5.tsv", header = 1)[,1]
net_filter <- net %>% dplyr::filter(source %in% intersection | target %in% intersection)

sss_regulators <- unique(net$source)
sss_tf_table <- net_filter %>% group_by(source) %>% summarise(n_targets = n()) %>% 
  filter(n_targets >= 3)

pdf("figures/Sup_Fig4A.pdf", height = 4, width = 4)
ggplot(sss_tf_table, aes(x = reorder(source, n_targets), y = n_targets)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Transcription Factor", y = "Number of SSS target genes",
       title = "Upstream TFs of the SSS gene signature") +
  theme_minimal()
dev.off()

library(igraph)
library(ggraph)
library(tidyverse)

net_plot <- net_filter %>%
  filter(source %in% sss_tf_table$source) %>% 
  mutate(reg_type = ifelse(mor == 1, "Activation",
                           ifelse(mor == -1, "Repression", "Unknown")))

g <- graph_from_data_frame(net_plot, directed = TRUE)

node_df <- data.frame(name = V(g)$name)

node_df <- node_df %>%
  mutate(node_type = case_when(
    name %in% intersection ~ "SSS gene",
    name %in% net$source ~ "TF",
    TRUE ~ "Other target"
  ))

set.seed(663)
layout <- create_layout(g, layout = "kk") %>%
  left_join(node_df, by = "name")

pdf("figures/Fig4A.pdf", height = 6, width = 8)
ggraph(layout) +
  geom_edge_link(aes(color = reg_type),
                 arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                 end_cap = circle(3, 'mm'),
                 edge_width = 0.8, alpha = 0.6) +
  geom_node_point(aes(color = node_type), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, 
                 point.padding = unit(0.8, "lines")) +
  scale_edge_color_manual(values = c("Activation" = "forestgreen", "Repression" = "firebrick")) +
  scale_color_manual(values = c(
    "TF" = "steelblue",
    "Other target" = "grey30",
    "SSS gene" = "orange"
  )) +
  theme_void() +
  labs(title = "Known TF–Target Regulatory Network of SSS Genes",
       edge_color = "Regulation",
       color = "Node type")
dev.off()

interest_tf <- unique(node_df$name[node_df$node_type=="TF"])
write.table(interest_tf, file = "tables/TableS6.tsv", quote = F, row.names = F)

## Figure 4B: decoupleR of interested TFs, to ensure these TF have higher activiti in S2 BMDM------------
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
net <- get_collectri(organism='human', split_complexes=FALSE)

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
  summarise(mean = mean(score))

top_acts_mat <- df %>%
  filter(source %in% interest_tf) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

palette_length = 100
my_color = colorRampPalette(c("darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

top_acts_mat <- top_acts_mat[paste0("S",0:4),]

library(pheatmap)
pdf("figures/Fig4B.pdf", height = 4, width = 6)
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
         main = "TF activity across BMDM subclusters (decoupleR)")
dev.off()

## Figure 4B: bulk RNA-seq validation--------
library(edgeR)
library(tidyverse)
library(readxl)

# Step 1: Read count data
count_data <- read_excel("data/bulk_RNA/Quantification_default_0_0_none2025-05-15-07-52-23.xlsx")

# Step 2: Prepare count matrix
counts <- count_data %>%
  column_to_rownames("gene_id") %>%
  select(starts_with("NC_"), starts_with("Erastin_")) %>%
  as.matrix()

group <- factor(c(rep("Erastin", 3),rep("NC", 3)))
dge <- DGEList(counts = counts, group = group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

et <- exactTest(dge)

score_list <- readRDS("result/gene_list_Fig3A.rds")

res <- topTags(et, n = Inf)$table %>%
  rownames_to_column("gene_id") %>%
  left_join(count_data %>% select(gene_id, gene_name), by = "gene_id") %>% 
  mutate(group = case_when(
    gene_name %in% intersection ~ "SSS gene",
    gene_name %in% interest_tf ~ "TF",
    gene_name %in% c(score_list$Ferroptosis_neg, score_list$Ferroptosis_dual) ~ "Ferroptosis_related",
    !gene_name %in% c(intersection, interest_tf, score_list$Ferroptosis_neg, score_list$Ferroptosis_dual) ~ "Other"
  ))

res_intersection <- res %>%
  dplyr::filter(gene_name %in% c(intersection, interest_tf, score_list$Ferroptosis_neg, score_list$Ferroptosis_dual) &
                  FDR < 0.05)

pdf("figures/Fig4C.pdf", height = 5, width = 5)
ggplot(res_intersection, aes(x = reorder(gene_name, logFC), y = logFC, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TF" = "steelblue", "SSS gene" = "orange", "Ferroptosis_related" = "forestgreen")) +
  coord_flip() +
  labs(title = "Significant expression changes in erastin-induced THP-1",
       x = "Gene", y = "log2 Fold Change") +
  theme_minimal()
dev.off()

## Figure 4C: bulk ATAC-seq validation--------
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
      SYMBOL %in% intersection ~ "SSS gene",
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
  "SSS gene" = "orange",
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

## Figure 4D: fingerprint of ATAC-seq-------
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

tf_motif <- query(MotifDb, "ETV4")
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
df_long2 <- pivot_longer(df, cols = -Position, names_to = "Group", values_to = "Signal")


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
ggplot(df_long2, aes(x = Position, y = Signal, color = Group)) +
  geom_line(size = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TF Footprint (ETV4)",
    x = "Position relative to motif (bp)",
    y = "Insertion Signal"
  ) +
  scale_color_manual(values = c("Erastin" = "#D95F02", "NC" = "#1B9E77"))
dev.off()

## Figure 4E: UMAP of scATAC-seq data-------------
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

## Figure 4F: Motif Deviation Score of scATAC---------------
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

saveRDS(integrated, file = "result/XH02_signac_integration.rds")

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

sss_tss <- annotations[annotations$gene_name %in% intersection]
promoter_regions <- promoters(sss_tss, upstream = 2000, downstream = 500)

motif_object <- Motifs(integrated)
motif_matrix <- motif_object@data[,rownames(motif_mat)]

peak_gr <- granges(integrated[["peaks"]])
names(peak_gr) <- rownames(integrated[["peaks"]])
peaks_in_promoter <- subsetByOverlaps(peak_gr, sss_promoter_region)

motif_promoter_matrix <- motif_matrix[rownames(motif_matrix) %in% names(peaks_in_promoter), ]

motif_sums <- Matrix::colSums(motif_promoter_matrix)
motif_keep <- as.vector(which(motif_sums > 0))

motif_mat_filtered <- motif_mat[motif_keep,]
gene_mat_filtered <- gene_mat[motif_keep,]

library(Matrix)
motif_avg_score <- motif_promoter_matrix[, motif_keep]
motif_score_df <- as.data.frame(as.matrix(motif_avg_score))
motif_score_df$peak <- rownames(motif_promoter_matrix)

peak_to_gene <- findOverlaps(peaks_in_promoter, sss_promoter_region)
peak_gene_map <- data.frame(
  peak = names(peaks_in_promoter)[queryHits(peak_to_gene)],
  gene = sss_promoter_region$gene_name[subjectHits(peak_to_gene)]
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

## Figure 4G: regulatory function of intersection genes------------
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

write.table(net_plot, file = "tables/Fig5A.tsv", quote = F, sep = "\t")
write.table(top_acts_mat, file = "tables/Fig5B.tsv", quote = F, sep = "\t")
write.table(integrated@meta.data[,c("grade","celltype")], file = "tables/Fig5C.tsv", quote = F, sep = "\t")
write.table(integrated@reductions$umap@cell.embeddings, file = "tables/Fig5D.tsv", quote = F, sep = "\t")
write.table(df_long1, file = "tables/Fig5E.tsv", quote = F, sep = "\t")
write.table(df_long2, file = "tables/Fig5F.tsv", quote = F, sep = "\t")



