BMDM <- readRDS("result/BMDM.rds")
Idents(BMDM) <- "seurat_clusters"
markers_df_BMDM2 <- FindAllMarkers(BMDM, only.pos = TRUE, test.use = "wilcox")

markers_BMDM2 <- markers_df_BMDM2$gene[markers_df_BMDM2$cluster==2]
length(markers_BMDM2)

Idents(BMDM) <- "grade"
markers_df_grade <- FindAllMarkers(BMDM, only.pos = TRUE, test.use = "wilcox")
table(markers_df_grade$cluster)
summary(markers_df_grade$avg_log2FC)

markers_older <- markers_df_grade$gene[markers_df_grade$cluster=="older"]
length(markers_older)

markers_common <- intersect(markers_older, markers_BMDM2)
markers_common <- markers_common[-grep("MT",markers_common)]
write.table(markers_common,file="result/BMDM_NES_signature.csv", row.names = F, col.names = F, quote = F)

## SPRING
MNES <- c("ANXA1", "ANXA2", "ARF6", "FN1", "HK2", "LDHA", "LYZ", "MPC2", "PDIA3", "PGK1",
          "RAB1A", "RETN", "S100A4", "S100A9", "SDCBP", "SOCS3", "TGFBI", "TIMP1")

library(clusterProfiler)
wp2gene <- read.gmt("./reference/h.all.v2023.1.Hs.symbols.gmt")
wp2gene <- lapply(names(wp2gene), function(x) data.frame(term=x, gene=wp2gene[[x]]))
wp2gene <- Reduce(rbind,wp2gene)
head(wp2gene)
ewp <- enricher(MNES, TERM2GENE = wp2gene, pvalueCutoff=0.05)

pdf("plot_new/GO_MNES.pdf", height = 3, width = 5)
barplot(ewp, showCategory=30, font.size = 8) + ggtitle("BMDM natural evolution signature")
dev.off()

ewp@result[grep("HYPOXIA",ewp@result$Description),]

cluster_color <- readRDS("color/ClusterColor.rds")
c <- readRDS("color/GradeColor.rds")
BMDM = Seurat::AddModuleScore(BMDM, features = list(MNES), name='NES_BMDM')
df<-data.frame(NES_BMDM = BMDM$NES_BMDM1, grade = BMDM$grade, cluster = BMDM$seurat_clusters)

pdf("plot_new/MNES_score.pdf", height = 4, width = 4)
FeaturePlot(BMDM, features = 'NES_BMDM1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")
ggboxplot(df, x = "cluster", y = "NES_BMDM",
  color = "cluster", palette = cluster_color)+ ylab("NES_BMDM")+
  stat_compare_means(comparisons = split(t(combn(as.character(0:4),2)),1:10), 
                     label = "p.signif", hide.ns = T)+NoLegend()
ggboxplot(df, x = "grade", y = "NES_BMDM",
          color = "grade", palette = grade_color)+ ylab("NES_BMDM")+
  stat_compare_means(comparisons = list(c("older","younger")), 
                     label = "p.signif", hide.ns = T)+NoLegend()
dev.off()





