library(Seurat)

NJ01_1<-readRDS("result/NJ01_1.rds")
NJ01_2<-readRDS("result/NJ01_2.rds")
NJ02_1<-readRDS("result/NJ02_1.rds")
NJ02_2<-readRDS("result/NJ02_2.rds")
TT01_1<-readRDS("result/TT01_1.rds")
TT01_2<-readRDS("result/TT01_2.rds")
TT02_1<-readRDS("result/TT02_1.rds")
TT02_2<-readRDS("result/TT02_2.rds")
XH01_1<-readRDS("result/XH01_1.rds")
XH01_2<-readRDS("result/XH01_2.rds")

library(RColorBrewer)
names<-c("NJ02_1","NJ02_2","TT01_1","TT01_2","TT02_1","TT02_2","NJ01_1","NJ01_2","XH01_1","XH01_2")
SampleColor <- brewer.pal(10, 'Paired')[c(2,1,4,3,6,5,8,7,10,9)]
names(SampleColor)<-names

GradeColor <- c("younger" = ,
                "older" = )

merge <- merge(x=NJ01_1,y=c(NJ01_2,NJ02_1,NJ02_2,TT01_1,TT01_2,TT02_1,TT02_2,XH01_1,XH01_2))
merge
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge)
merge <- ScaleData(merge)
merge <- RunPCA(merge, npcs = 50)
merge <- RunUMAP(merge, dims = 1:40)
merge <- FindNeighbors(merge, dims = 1:40)
merge <- FindClusters(merge, resolution = 0.8)
merge$sample <- sub(".*\\_(.*)","\\1",colnames(merge))
table(merge$sample)
merge$sample[merge$sample==1] <- "NJ01_1"
merge$sample[merge$sample==2] <- "NJ01_2"
merge$sample[merge$sample==3] <- "NJ02_1"
merge$sample[merge$sample==4] <- "NJ02_2"
merge$sample[merge$sample==5] <- "TT01_1"
merge$sample[merge$sample==6] <- "TT01_2"
merge$sample[merge$sample==7] <- "TT02_1"
merge$sample[merge$sample==8] <- "TT02_2"
merge$sample[merge$sample==9] <- "XH01_1"
merge$sample[merge$sample==10] <- "XH01_2"
saveRDS(merge,file = "result/merged_withXHsample.rds")

merge$grade <- "older"
merge$grade[sub(".*\\_(.*)","\\1",merge$sample)==2] <- "younger"

pdf("plot_new/merged_withXHsample.pdf")
DimPlot(merge, label = T)
DimPlot(merge, group.by = "sample", cols = SampleColor)
DimPlot(merge, group.by = "grade", cols = GradeColor)
ggplot(as.data.frame(table(merge$sample)))+geom_col(aes(Var1,Freq,fill=Var1))+
  scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")
DimPlot(merge, group.by = "celltype_mid_main")
DimPlot(merge, group.by = "celltype_bped_main")
DimPlot(merge, group.by = "celltype_hpca_main")
DimPlot(merge, group.by = "celltype_iced_main")
dev.off()


# ifnb.list<-list(NJ01_1,NJ01_2,NJ02_1,NJ02_2,TT01_1,TT01_2,TT02_1,TT02_2,high,low)
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# immune.combined <- IntegrateData(anchorset = immune.anchors)
# saveRDS(immune.combined, file = "result/combined.rds")
# 
# DefaultAssay(immune.combined) <- "integrated"
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# immune.combined$sample <- sub(".*\\_(.*)","\\1",colnames(immune.combined))
# table(immune.combined$sample)
# immune.combined$sample[immune.combined$sample==1] <- "NJ01_1"
# immune.combined$sample[immune.combined$sample==2] <- "NJ01_2"
# immune.combined$sample[immune.combined$sample==3] <- "NJ02_1"
# immune.combined$sample[immune.combined$sample==4] <- "NJ02_2"
# immune.combined$sample[immune.combined$sample==5] <- "TT01_1"
# immune.combined$sample[immune.combined$sample==6] <- "TT01_2"
# immune.combined$sample[immune.combined$sample==7] <- "TT02_1"
# immune.combined$sample[immune.combined$sample==8] <- "TT02_2"
# immune.combined$sample[immune.combined$sample==9] <- "high"
# immune.combined$sample[immune.combined$sample==10] <- "low"
# 
# pdf("plot_new/combined.pdf")
# DimPlot(immune.combined)
# DimPlot(immune.combined, group.by = "sample", cols = SampleColor)
# ggplot(as.data.frame(table(immune.combined$sample)))+geom_col(aes(Var1,Freq,fill=Var1))+
#   scale_fill_manual(values = SampleColor)+theme_bw()+xlab("Sample")+ylab("Number of cells")
# DimPlot(immune.combined, group.by = "celltype_mid_main")
# DimPlot(immune.combined, group.by = "celltype_iced_main")
# dev.off()
# 
# saveRDS(immune.combined,file = "result/combined.rds")

