library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)

merge <- readRDS("result/merged_withXHsample.rds")
pbmc <- merge
#pbmc$batch <- pbmc$sample
pbmc$batch <- "NJ"
pbmc$batch[pbmc$sample=="TT01_1"|pbmc$sample=="TT01_2"|pbmc$sample=="TT02_1"|pbmc$sample=="TT02_2"] <- "TT"
pbmc$batch[pbmc$sample=="XH01_1"|pbmc$sample=="XH01_2"] <- "XH"

table(pbmc$batch)
set.seed(10)
pbmc <- pbmc %>% 
  RunHarmony("batch", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

pbmc <- pbmc %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 
pbmc <- pbmc %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

library(RColorBrewer)
names<-c("NJ02_1","NJ02_2","TT01_1","TT01_2","TT02_1","TT02_2","NJ01_1","NJ01_2","XH01_1","XH01_2")
SampleColor <- brewer.pal(8, 'Paired')[c(2,1,4,3,6,5,8,7,10,9)]
names(SampleColor)<-names

GradeColor <- brewer.pal(10, 'Paired')[9:10]
names(GradeColor) <- c("younger","older")

pdf("plot_new/merged_withXHsample_rmbatch.pdf")
DimPlot(pbmc, label = TRUE)
DimPlot(pbmc, group.by = "batch")
DimPlot(pbmc, group.by = "sample", cols = SampleColor)
DimPlot(pbmc, group.by = "grade", cols = GradeColor)
DimPlot(pbmc, group.by = "celltype_mid_main")
DimPlot(pbmc, group.by = "celltype_bped_main")
DimPlot(pbmc, group.by = "celltype_hpca_main")
DimPlot(pbmc, group.by = "celltype_iced_main")
dev.off()

saveRDS(pbmc,"result/merged_withXHsample_rmbatch.rds")
