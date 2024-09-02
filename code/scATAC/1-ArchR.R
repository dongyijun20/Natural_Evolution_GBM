library(ArchR)
library(parallel)
addArchRThreads(threads = 4)
addArchRGenome("hg38")
h5disableFileLocking()
sampleID_fileID<-"data/N20240290-01-ATAC/0_CellRanger/fragments.tsv.gz"
names(sampleID_fileID) <- c("XH01")
ArrowFiles_sample <- createArrowFiles(
  inputFiles = sampleID_fileID,
  outputNames = names(sampleID_fileID),
  sampleNames = names(sampleID_fileID),
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  minTSS =4, #Dont set this too high because you can always increase later
  minFrags =1000,
  maxFrags = 100000,
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE)
h5disableFileLocking()
doubScores <- addDoubletScores(
  input = ArrowFiles_sample,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

h5disableFileLocking()
proj_name <- ArchRProject(
  ArrowFiles = ArrowFiles_sample,
  outputDirectory = "ArchROutput1",
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.sdo
)
genescore<-getMatrixFromProject(proj_name,useMatrix = "GeneScoreMatrix")
save(genescore,file="result/genescore.rdata")

h5disableFileLocking()
proj_name_new <- filterDoublets(proj_name)

proj_name_new <- addCellColData(ArchRProj = proj_name, data = case_when(sub(".*\\-(.*)","\\1",proj_name$cellNames)==1 ~ "XH01_2",
                                                                        sub(".*\\-(.*)","\\1",proj_name$cellNames)==1 ~ "XH01_1"),
                            cells = proj_name$cellNames, name = "Sample", force = TRUE)

proj_name_new <- addIterativeLSI(
  ArchRProj = proj_name_new,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

# proj_name_new <- addHarmony(
#   ArchRProj = proj_name_new,
#   reducedDims = "IterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample",
#   force =TRUE
# )
h5disableFileLocking()
proj_name_new<-addClusters(
  input = proj_name_new,
  reducedDims = "IterativeLSI",
  method = "seurat",
  name = "Clusters",
  resolution =2,
  force =TRUE)

h5disableFileLocking()
proj_name_new<-addUMAP(
  ArchRProj = proj_name_new,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force =TRUE
)
h5disableFileLocking()
# proj_name_new<-addTSNE(
#   ArchRProj = proj_name_new,
#   reducedDims = "Harmony",
#   name = "TSNE",
#   perplexity = 100,
#   force =TRUE
# )
# h5disableFileLocking()

saveArchRProject(proj_name_new,outputDirectory ="result/XH01_scATAC", dropCells = T)

p1 <- ArchR::plotEmbedding(ArchRProj = proj_name_new, colorBy = "cellColData", name = "Sample", embedding = "UMAP", randomize = T)
p2 <- ArchR::plotEmbedding(ArchRProj = proj_name_new, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", randomize = T)
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj_name_new, addDOC = FALSE, width = 5, height = 5)

p1 <- plotGroups(
  ArchRProj = proj_name_new, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = proj_name_new, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj_name_new, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = proj_name_new, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_name_new, addDOC = FALSE, width = 4, height = 4)

proj_name_new<-loadArchRProject("result/all_proj_QC")

h5disableFileLocking()
markersGS <- getMarkerFeatures(
  ArchRProj = proj_name_new, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

proj_name_new <- addImputeWeights(proj_name_new)


markerGenes  <- c(
  "CLDN5","PECAM1", #endothelial
  "RGS5", "ACTA2", #vascular
  "ISLR","CTHRC1", # mesenchymal stromal cell
  "JCHAIN", "MZB1", #B-Cell
  "AIF1", "LYZ", #metastasis-associated macrophages
  "CD1C", "CLEC10A", #DC cell
  "GFAP", "S100B", #reactive astrocyte
  "CD3D", "IL7R" #T-Cells
)


p <- plotEmbedding(
  ArchRProj = proj_name_new, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj_name_new)
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = proj_name_new, 
        addDOC = FALSE, width = 5, height = 5)

markerGenes <- sapply(markerList, function(x) x$name[1:2])
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj_name_new, addDOC = FALSE)

seRNA <- readRDS("../neurosurgery/metastasis/landscape/GSM5645891_seurat.rds")
proj_name_new <- addGeneIntegrationMatrix(
  ArchRProj = proj_name_new, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
cM <- as.matrix(confusionMatrix(proj_name_new$Clusters, proj_name_new$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
df<-cbind(preClust, rownames(cM)) #Assignments
table(proj_name_new$Sample,proj_name_new$predictedGroup_Un)

proj_name_new$celltype<-sapply(proj_name_new$Clusters, function(x) df[df[,2]==x,1])
table(proj_name_new$celltype)
proj_name_new$celltype <- sub("(.*)\\-.*","\\1",proj_name_new$celltype)
proj_name_new$celltype <- sub("(.*)\\:.*","\\1",proj_name_new$celltype)

proj_name_new$celltype[proj_name_new$Clusters%in%paste0("C",c(1:15))]<-"Malignant"

pal <- paletteDiscrete(values = seRNA$cell_type)
p1 <- plotEmbedding(
  proj_name_new, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)
pal <- paletteDiscrete(values = proj_name_new$celltype)
p2 <- plotEmbedding(
  proj_name_new, 
  colorBy = "cellColData", 
  name = "celltype", 
  pal = pal
)
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj_name_new, addDOC = FALSE, width = 5, height = 5)

proj_name_new <- addGroupCoverages(ArchRProj = proj_name_new, groupBy = "celltype")

pathToMacs2 <- "/home/tianchen/.conda/envs/python3.6/bin/macs2"
proj_name_new <- addReproduciblePeakSet(
  ArchRProj = proj_name_new, 
  groupBy = "celltype", 
  pathToMacs2 = pathToMacs2
)

saveArchRProject(proj_name_new,outputDirectory ="result/all_proj_QC", dropCells = T)

projHeme5 <- addPeakMatrix(proj_name_new)

getAvailableMatrices(projHeme5)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

pv <- markerPlot(seMarker = markersPeaks, name = "Malignant", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv

markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Malignant",
  bgdGroups = "T:CD4+"
)

projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
plotPDF(ggUp, ggDo, name = "Malignant-vs-T:CD4-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
