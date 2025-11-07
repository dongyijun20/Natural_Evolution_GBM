## Figure 5A: calculation for the Venn plot--------
library(Seurat)
library(ArchR)
library(dplyr)

log2FC_filter = 1
p_val_filter = 0.05

# standard 1
marker_list <- readRDS("result/BMDM_markers_filter.rds")
BMDM_subcluster_markers <- marker_list$S2
length(BMDM_subcluster_markers) #46

# standard 2
TAM <- readRDS("result/TAM_WJS.rds")
Idents(TAM) <- "position"
TAM_necrotic_markers_raw <- FindMarkers(TAM, ident.1 = "HS", ident.2 = "SX", only.pos = T)
TAM_necrotic_markers <- rownames(TAM_necrotic_markers_raw)[TAM_necrotic_markers_raw$avg_log2FC>=log2FC_filter&
                                                             TAM_necrotic_markers_raw$p_val<p_val_filter]
length(TAM_necrotic_markers) #52

# standard 3
ferroptosis_suppressor <- read_csv("reference/ferroptosis_suppressor.csv")
length(unique(ferroptosis_suppressor$symbol)) #238

# Intersect with external marker lists
intersection <- Reduce(intersect, list(BMDM_subcluster_markers, TAM_necrotic_markers))
intersection_core <- Reduce(intersect, list(BMDM_subcluster_markers, ferroptosis_suppressor$symbol))

length(intersection) #14
length(intersection_core) #4
intersection
intersection_core

write.table(intersection, file = "tables/TableS5.tsv", quote = F, row.names = F)
intersection <- read.table("tables/TableS5.tsv", header = T)[,1]

library(VennDiagram)

pdf("figures/Fig5A.pdf", height = 3.5, width = 4)
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = length(BMDM_subcluster_markers),
  area2 = length(TAM_necrotic_markers),
  cross.area = length(intersection),
  category = c("S2 BMDM markers", "PC vs IT upregulated markers"),
  fill = c("coral", "lightblue"),
  cat.pos = c(0, 0),        # Positions in degrees (0 is center-top)
  cat.dist = c(0.02, 0.02), # Smaller distance = closer to center
  cat.cex = 1.2             # Optional: change label size
)
dev.off()
    
## Figure 5A: Survival Analysis with Expression-----------
library(survival)
library(survminer)
library(tidyverse)
library(TCGAbiolinks)
library(GSVA)

# Query for gene expression data (HTSeq counts or FPKM)
query_expression <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_expression)
expression_data <- GDCprepare(query_expression)
saveRDS(expression_data, "reference/TCGA/TCGA_GBM_expression_data.rds")

# Extract the expression matrix and rename rows
fpkm_data_TCGA <- assay(expression_data, "fpkm_unstrand")
rownames(fpkm_data_TCGA) <- rowData(expression_data)$gene_name

# Read gene set
intersection <- read.table("tables/TableS5.tsv", header = 1)[,1]
gene_set_list <- list(SSS = as.character(intersection))

# Run ssGSEA
ssgsea_scores <- gsva(fpkm_data_TCGA, gene_set_list, method = "ssgsea", verbose = FALSE)

# Add to clinical data
clinical_data_TCGA <- as.data.frame(colData(expression_data)) 
clinical_data_TCGA <- clinical_data_TCGA %>% 
  mutate(SSS = ssgsea_scores["SSS", clinical_data_TCGA$barcode])

# Group by quartiles
clinical_data_TCGA <- clinical_data_TCGA %>%
  mutate(group = case_when(
    SSS <= quantile(SSS, 1/3, na.rm = TRUE) ~ "Low",
    SSS >= quantile(SSS, 2/3, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_ 
  ))

# Prepare survival info
clinical_data_TCGA <- clinical_data_TCGA %>%
  mutate(OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
         OS_status = ifelse(vital_status == "Dead", 1, 0))

# Run survival analysis
surv_obj <- Surv(time = clinical_data_TCGA$OS_time, event = clinical_data_TCGA$OS_status)
fit <- survfit(surv_obj ~ group, data = clinical_data_TCGA)

# Plot survival curve
pdf("figures/Fig5A.pdf", height = 5, width = 4.5)
ggsurvplot(fit, data = clinical_data_TCGA,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "Survival Analysis with SSS ssGSEA scoring in TCGA GBM",
           palette = c("#8B235A","#005787"))
dev.off()

cox_data_TCGA <- clinical_data_TCGA %>% 
  mutate(gender = factor(gender),
         paper_IDH.status = factor(paper_IDH.status),
         paper_MGMT.promoter.status = factor(paper_MGMT.promoter.status),
         SSS_score_group = factor(group, levels = c("Low", "High")))

cox_model <- coxph(Surv(OS_time, OS_status) ~ 
                     SSS_score_group + age_at_diagnosis + gender + paper_IDH.status + paper_MGMT.promoter.status, 
                   data = cox_data_TCGA, control = coxph.control(iter.max = 50))
library(forestmodel)

pdf("figures/Sup_Fig5A.pdf", height = 5, width = 9)
forest_model(cox_model)
dev.off()

## Figure 5B: Correlation of immune score----------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Subset your expression data
high_samples <- na.omit(rownames(clinical_data_TCGA)[clinical_data_TCGA$group == "High"])
low_samples <- na.omit(rownames(clinical_data_TCGA)[clinical_data_TCGA$group == "Low"])

# Calculate average expression in high vs low group
expr_high <- rowMeans(fpkm_data_TCGA[, high_samples])
expr_low <- rowMeans(fpkm_data_TCGA[, low_samples])

# Compute log2 fold-change
logFC <- log2(expr_high + 1) - log2(expr_low + 1)

# Prepare ranked gene list (named vector: gene_symbol -> logFC)
gene_list <- sort(logFC, decreasing = TRUE)

gene_df <- bitr(names(gene_list), fromType = "SYMBOL", 
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Replace gene symbols with Entrez IDs
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

ego <- gseGO(geneList = gene_list,
             OrgDb = org.Hs.eg.db,
             ont = "BP",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             verbose = TRUE)

library(enrichplot)
grep("immune", ego@result$Description, value = TRUE)
ego@result[ego@result$Description == "negative regulation of immune response",]

pdf("figures/Sup_Fig5F.pdf", height = 4, width = 4)
gseaplot2(ego, geneSetID = ego@result$ID[ego@result$Description == "negative regulation of immune response"], 
          title = "negative regulation of immune response")
gseaplot2(ego, geneSetID = ego@result$ID[ego@result$Description == "immune response-inhibiting signal transduction"], 
          title = "immune response-inhibiting signal transduction")
dev.off()

## Figure 5C: Survival Analysis with CGGA-------------
fpkm_data1 <- read_table("reference/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
fpkm_data2 <- read_table("reference/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")

# Keep shared genes and merge expression matrices
identical_genes <- intersect(fpkm_data1$Gene_Name, fpkm_data2$Gene_Name)
fpkm_data1 <- fpkm_data1 %>% column_to_rownames("Gene_Name")
fpkm_data2 <- fpkm_data2 %>% column_to_rownames("Gene_Name")
fpkm_data <- cbind(fpkm_data1[identical_genes, ], fpkm_data2[identical_genes, ])
saveRDS(fpkm_data, file = "reference/CGGA/CGGA_GBM_expression_data.rds")

# Define your signature gene set
gene_set_list <- list(SSF = as.character(intersection))  # ensure `intersection` is defined

# Run ssGSEA
ssgsea_scores <- gsva(as.matrix(fpkm_data), gene_set_list, method = "ssgsea", verbose = FALSE)

# Load and merge clinical data
clinical_data1 <- read.table("reference/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data2 <- read.table("reference/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data <- rbind(clinical_data1, clinical_data2) %>%
  dplyr::filter(Histology == "GBM")

# Match ssGSEA score to clinical samples
clinical_data$signature_score <- ssgsea_scores["SSF", clinical_data$CGGA_ID]

# Group patients by quartiles of signature score
clinical_data <- clinical_data %>%
  mutate(group = case_when(
    signature_score <= quantile(signature_score, 1/3, na.rm = TRUE) ~ "Low",
    signature_score >= quantile(signature_score, 2/3, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_
  ))

# Prepare survival variables
surv_obj <- Surv(time = clinical_data$OS, event = clinical_data$Censor..alive.0..dead.1.)
fit <- survfit(surv_obj ~ group, data = clinical_data)

# Plot KM survival curves
pdf("figures/Fig5B.pdf", height = 5, width = 4.5)
ggsurvplot(fit, data = clinical_data,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "Survival Analysis with SSS ssGSEA scoring in CGGA GBM",
           palette = c("#FDDE00", "#1068B1"))
dev.off()

## Figure 5C: Survival Analysis with CIBERSORTx Deconvolution---------
# Save bulk expression data 
fpkm_data_TCGA <- cbind(GeneSymbol = rownames(fpkm_data_TCGA), fpkm_data_TCGA)
fpkm_data_TCGA <- fpkm_data_TCGA[!duplicated(rownames(fpkm_data_TCGA)),]
dim(fpkm_data_TCGA)

write.table(fpkm_data_TCGA, "reference/cybersort/bulk_expression.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save single-cell reference matrix 
malignant <- readRDS("result/malignant_withXHsample_harmony.rds")
myeloid <- readRDS("result/myeloid_withXHsample_harmony.rds")

reference_matrix1 <- as.matrix(GetAssayData(BMDM, assay = "RNA", slot = "data"))
reference_matrix2 <- as.matrix(GetAssayData(subset(myeloid, TAM_type != "BMDM"), assay = "RNA", slot = "data"))
reference_matrix3 <- as.matrix(GetAssayData(malignant, assay = "RNA", slot = "data"))
colnames(reference_matrix1) <- BMDM$subcluster
colnames(reference_matrix2) <- subset(myeloid, TAM_type != "BMDM")$TAM_type
colnames(reference_matrix3) <- malignant$subtype

intersected_genes <- intersect(rownames(reference_matrix1), rownames(reference_matrix3))

reference_matrix <- cbind(reference_matrix1[intersected_genes,],reference_matrix2[intersected_genes,],reference_matrix3[intersected_genes,])
reference_matrix <- cbind(GeneSymbol = rownames(reference_matrix), reference_matrix)

selected_cells <- sapply(unique(colnames(reference_matrix)), function(x) {
  cell_indices <- which(colnames(reference_matrix) == x)
  sample_size <- min(200, length(cell_indices))
  sample(cell_indices, sample_size)
})

reference_matrix <- reference_matrix[,unlist(selected_cells)]
dim(reference_matrix)

write.table(reference_matrix, "reference/cybersort/single_cell_reference_matrix.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

clinical_data
# Merge CIBERSORTx results with clinical data
cibersortx_results <- read.csv("result/CIBERSORTx_Job12_Results.csv")
merged_data_cibersort <- cibersortx_results %>%
  inner_join(clinical_data_TCGA, by = c("Mixture"="barcode"))

# Define high and low S2 BMDM groups
merged_data_cibersort <- merged_data_cibersort %>%
  mutate(S2_cibersort_group = case_when(
    S2 <= quantile(S2, 1/3, na.rm = TRUE) ~ "Low",
    S2 >= quantile(S2, 2/3, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_ 
  ))

# Kaplan-Meier survival analysis
surv_obj_cibersort <- Surv(time = merged_data_cibersort$OS_time, event = merged_data_cibersort$OS_status)
fit_cibersort <- survfit(surv_obj_cibersort ~ S2_cibersort_group, data = merged_data_cibersort)

pdf("figures/Fig5C.pdf", height = 5, width = 4.5)
ggsurvplot(fit_cibersort, data = merged_data_cibersort,
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "S2 BMDM Proportion",
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "S2 BMDM Proportion in TCGA GBM (deconvoluted)",
           palette = c("#D55E00","#0072B2"))
dev.off()

## Figure 5D: GLASS cohort for longitudinal validation---------------
# Load expression and metadata
gene_tpm_matrix <- read_tsv("reference/GLASS/gene_tpm_matrix_all_samples.tsv")
analysis_rnaseq_pairs <- read_csv("reference/GLASS/analysis_rnaseq_pairs.csv") %>%
  mutate(across(starts_with("tumor_barcode"), ~ gsub("-", ".", .)))

clinical_cases <- read_csv("reference/GLASS/clinical_cases.csv")
clinical_surgeries <- read_csv("reference/GLASS/clinical_surgeries.csv")

# Merge clinical info
clinical_GLASS <- clinical_cases %>%
  inner_join(clinical_surgeries, by = "case_barcode") %>%
  dplyr::filter(grade == "IV") %>%
  dplyr::select(case_barcode, case_project, case_source, case_sex, case_age_diagnosis_years,
         case_vital_status, case_overall_survival_mo, grade, idh_status, histology) %>%
  distinct()

merged_data <- analysis_rnaseq_pairs %>%
  inner_join(clinical_GLASS, by = "case_barcode")

# Ensure your gene set (intersection) is a character vector of gene symbols
# e.g., intersection <- read.table("tables/TableS5.tsv", header = FALSE)[,1]
gene_set_list <- list(signature = as.character(intersection))

# Prepare expression matrix: rows = genes, columns = samples
tpm_matrix <- gene_tpm_matrix %>%
  column_to_rownames("Gene_symbol") %>%
  as.matrix()

# Optional: log2 transform
tpm_matrix_log <- log2(tpm_matrix + 1)

# Compute ssGSEA score per sample
ssgsea_result <- gsva(tpm_matrix_log, gene_set_list, method = "ssgsea", verbose = FALSE)

# Reshape to long format
signature_scores <- data.frame(
  tumor_barcode = colnames(ssgsea_result),
  signature = as.numeric(ssgsea_result[1, ])
)

# Merge ssGSEA scores into paired data
final_data <- merged_data %>%
  inner_join(signature_scores, by = c("tumor_barcode_a" = "tumor_barcode")) %>%
  rename(signature_a = signature) %>%
  inner_join(signature_scores, by = c("tumor_barcode_b" = "tumor_barcode")) %>%
  rename(signature_b = signature) 

quartile_a <- quantile(final_data$signature_a, probs = c(0.5), na.rm = TRUE)
quartile_b <- quantile(final_data$signature_b, probs = c(0.5), na.rm = TRUE)

final_data <- final_data %>%
  mutate(
    tumor_type_a = ifelse(sample_type_a == "TP", "Primary", "Recurrent"),
    tumor_type_b = ifelse(sample_type_b == "TP", "Primary", "Recurrent"),
    survival_time = case_overall_survival_mo,
    survival_status = ifelse(case_vital_status == "dead", 1, 0),
    signature_group = case_when(
      signature_a <= signature_b ~ "Up",
      signature_a > signature_b ~ "Down",
      TRUE ~ NA
    ),
    diff = signature_b - signature_a,
    group_a = case_when(
      signature_a <= quartile_a[1] ~ "Low",
      signature_a >= quartile_a[1] ~ "High",
      TRUE ~ NA
    ),
    group_b = case_when(
      signature_b <= quartile_b[1] ~ "Low",
      signature_b >= quartile_b[1] ~ "High",
      TRUE ~ NA
    ),
    group_combined = paste0(group_a, "_to_", group_b)
  ) %>%
  arrange(diff) %>%
  mutate(case_barcode_ordered = paste0(case_barcode, "_", row_number())) %>%
  mutate(
    case_sex = factor(case_sex),
    idh_status = factor(idh_status)
  )

# Survival analysis
surv_obj <- Surv(time = final_data$survival_time, event = final_data$survival_status)
fit <- survfit(surv_obj ~ group_a + group_b, data = final_data)

pdf("figures/Fig5D.pdf", height = 6, width = 6)
ggsurvplot(fit, data = final_data,
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "Gene Set Change",
           xlab = "Survival Time (Months)",
           ylab = "Overall Survival Probability")
dev.off()

# Waterfall plot
pdf("figures/Sup_Fig5C.pdf", height = 4, width = 6)
ggplot(final_data, aes(x = factor(case_barcode_ordered, levels = case_barcode_ordered), y = diff, fill = signature_group)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Case Barcode",
    y = "Difference in Gene Set Score (Recurrent - Primary)",
    fill = "Change Group",
    title = "Waterfall Plot of Gene Set Score Differences"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
dev.off()


# Figure 5F: immune response from TCGA--------
TCIA_clinical <- read.table("reference/TCGA/TCIA-ClinicalData.tsv", sep = "\t")
colnames(TCIA_clinical) <- TCIA_clinical[1,]
TCIA_clinical <- TCIA_clinical[-1,]

merged_data <- merge(clinical_data_TCGA, TCIA_clinical,
                     by.x = "patient", by.y = "barcode",
                     all.x = TRUE)
merged_data$ips_ctla4_neg_pd1_pos <- as.numeric(merged_data$ips_ctla4_neg_pd1_pos)
merged_data$ips_ctla4_pos_pd1_pos <- as.numeric(merged_data$ips_ctla4_pos_pd1_pos)

ggboxplot(merged_data, x = "group", y = "ips_ctla4_pos_pd1_pos",
          color = "group") + 
  stat_compare_means(comparisons = list(c("High","Low")), label = "p.signif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

# Figure 5F: validation on larger cohorts--------------
GBM_intergrated_24 <- readRDS("reference/GBM.RNA.integrated.24.rds")
GBM_intergrated_24_TAM <- subset(GBM_intergrated_24, subset = anno_ident=="Macrophages")

seurat_workflow <- function(seu){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 50)
  seu <- RunUMAP(seu, dims = 1:40)
  seu <- FindNeighbors(seu, dims = 1:40)
  seu <- FindClusters(seu)
  seu
}

GBM_intergrated_24_TAM <- seurat_workflow(GBM_intergrated_24_TAM)

GBM_intergrated_24_TAM <- GBM_intergrated_24_TAM %>% 
  RunHarmony("Pt_number", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

GBM_intergrated_24_TAM <- GBM_intergrated_24_TAM %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.1) 

GBM_intergrated_24_TAM <- GBM_intergrated_24_TAM %>%
  RunUMAP(dims = 1:40, reduction = "harmony")

DimPlot(GBM_intergrated_24_TAM)
DimPlot(GBM_intergrated_24_TAM, group.by = "treatment_3")

BMDM <- readRDS("result/BMDM_withXHsample_harmony.rds")
BMDM <- seurat_workflow(BMDM)
anchors <- FindTransferAnchors(reference = BMDM, query = GBM_intergrated_24_TAM, dims = 1:30, reference.reduction = "pca")

# set UMAP models
umap_new_model <- list()
umap_new_model$n_epochs <- 500
umap_new_model$alpha <-1
umap_new_model$method <- "umap"
umap_new_model$negative_sample_rate <- 5
umap_new_model$gamma <- 1
umap_new_model$approx_pow <- 0
umap_new_model$n_neighbors <- 30
umap_new_model$metric$cosine <- list()
umap_new_model$embedding <- BMDM[["umap"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap_new_model$a <- ab_param["a"]
umap_new_model$b <- ab_param["b"]
BMDM[["umap"]]@misc$model <- umap_new_model

GBM_intergrated_24_TAM <- MapQuery(
  anchorset = anchors,
  query = GBM_intergrated_24_TAM,
  reference = BMDM,
  refdata = list(
    subcluster = "subcluster"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

SubclusterColor <- c('#ea9994','#f2c396','#9cd2ed','#86c7b4','#a992c0')
names(SubclusterColor) <- unique(BMDM$subcluster)

pdf("figures/Fig5E.pdf", height = 4.5, width = 5.5)
DimPlot(GBM_intergrated_24_TAM, group.by = "predicted.subcluster", label = T, cols = SubclusterColor)
dev.off()

saveRDS(GBM_intergrated_24_TAM, file = "result/BMDM_neoPD1_treatments.rds")

# Figure 5G: PD-1 response is related to S2 BMDM proportion
library(forcats)
library(ggplot2)
library(gridExtra)
GBM_intergrated_24_TAM@meta.data$treatment_3 <- factor(GBM_intergrated_24_TAM@meta.data$treatment_3, levels = c("responder","nonresponder","Rec","Untreated"))
CellInfo <- GBM_intergrated_24_TAM@meta.data
P1=CellInfo %>% ggplot(aes(x=treatment_3, fill=fct_rev(predicted.subcluster))) +scale_fill_manual(values = SubclusterColor)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(predicted.subcluster , fill=predicted.subcluster))+geom_bar(stat="count",colour = "black",width = 0.7)+  scale_fill_manual(values = SubclusterColor)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)

pdf("figures/Fig5F.pdf",width=7,height=6)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

GBM_intergrated_24_TAM <- AddModuleScore_UCell(GBM_intergrated_24_TAM, list("SSS"=intersection,"S2"=BMDM_subcluster_markers), name = NULL)
VlnPlot(GBM_intergrated_24_TAM, "SSS", group.by = "treatment_3")
VlnPlot(GBM_intergrated_24_TAM, "SSS", group.by = "predicted.subcluster")
VlnPlot(GBM_intergrated_24_TAM, "S2", group.by = "predicted.subcluster")

library(viridis)
library(scales)
vlnplot <- VlnPlot(GBM_intergrated_24_TAM, 
                   features = c(intersection, "SSS"), 
                   group.by = "treatment_3", 
                   stack = TRUE, 
                   cols = viridis_pal(option = "D")(length(c(intersection, "SSS")))) + 
  NoLegend()
pdf("figures/Fig5G.pdf", height = 4, width = 9)
vlnplot
dev.off()

## Figure 5H: correlation plot of res and non-res-----------
library(ggrepel)

# 1. Aggregate expression data (mean or median) per group
DefaultAssay(GBM_intergrated_24_TAM) <- "RNA"
GBM_intergrated_24_TAM <- NormalizeData(GBM_intergrated_24_TAM)
avg_exp <- AggregateExpression(
  GBM_intergrated_24_TAM,
  assays = "RNA",
  features = rownames(GBM_intergrated_24_TAM),
  group.by = c("treatment_3"),  # "treatment_3" = responder/non-responder
  slot = "data", 
  return.seurat = FALSE
)$RNA

group_sizes <- table(GBM_intergrated_24_TAM$treatment_3)
avg_exp <- sweep(avg_exp, 2, group_sizes[colnames(avg_exp)], FUN = "/")
range(avg_exp)

# 2. Reshape data for plotting (responder vs. non-responder)
library(tidyr)
lower_cutoff <- 2^-5

plot_data <- avg_exp %>% as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  mutate(
    responder_adj = log2(responder + 1), 
    nonresponder_adj = log2(nonresponder + 1), 
    to_label = ifelse(gene %in% intersection, gene, NA)
  )

# 3. Plot
library(ggrepel)
library(scales)

cor_plot <- ggplot(plot_data, aes(x = responder_adj, y = nonresponder_adj)) +
  # Plot non-intersection genes (light gray)
  geom_point(
    data = subset(plot_data, !gene %in% intersection),
    color = "gray90", 
    size = 1.5, 
    alpha = 0.4
  ) +
  # Plot intersection genes (highlighted)
  geom_point(
    data = subset(plot_data, gene %in% intersection),
    color = "#E64B35",
    size = 2.5,
    alpha = 0.8
  ) +
  # Label only intersection genes
  geom_text_repel(
    aes(label = to_label),
    color = "black",
    size = 3,
    max.overlaps = 100,
    box.padding = 0.5,
    segment.color = "grey50",
    min.segment.length = 0.2
  ) +
  # Reference line y = x
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "gray40",
    linewidth = 0.5
  ) +
  # Labels and theme
  labs(
    x = "log2(CPM + 1) in Responder",
    y = "log2(CPM + 1) in Non-responder"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save plot
pdf("figures/Fig5H.pdf", height = 5, width = 5)
cor_plot
dev.off()

# Figure 5I: drug response in SSS-high BMDM population--------------
BMDM_merged <- readRDS("result/BMDM_merged.rds")
BMDM_merged$subcluster[BMDM_merged$TAM_subtype=="other TAMs"] <- "other"
BMDM_merged$subcluster[BMDM_merged$TAM_subtype=="S2 BMDM"] <- "S2"
BMDM_merged$subcluster[BMDM_merged$subcluster!="S2"] <- "other"
BMDM_merged$subcluster[BMDM_merged$predicted.subcluster!="S2"] <- "other"
BMDM_merged$subcluster[BMDM_merged$predicted.subcluster=="S2"] <- "S2"
table(BMDM_merged$subcluster)

intersection <- read.table("tables/TableS5.tsv", header = T)[,1]

DimPlot(BMDM_merged, group.by = "subcluster", reduction = "tsne")

library(UCell)
BMDM_merged <-AddModuleScore_UCell(BMDM_merged, list(SSS=intersection), name = NULL)
FeaturePlot(BMDM_merged, "SSS", max.cutoff = "q90", min.cutoff = "q10", reduction = "tsne")

BMDM$label_SSS <- "Middle"
sample_list <- unique(BMDM$orig.ident)  

for (sample_id in sample_list) {
  cells_in_sample <- colnames(BMDM)[BMDM$orig.ident==sample_id]
  sss_scores <- BMDM$SSS[cells_in_sample]
  
  low_cutoff <- quantile(sss_scores, 1/3, na.rm = TRUE)
  high_cutoff <- quantile(sss_scores, 2/3, na.rm = TRUE)
  
  BMDM$label_SSS[cells_in_sample][sss_scores < low_cutoff] <- "Low"
  BMDM$label_SSS[cells_in_sample][sss_scores >= high_cutoff] <- "High"
}

DimPlot(BMDM, group.by="label_SSS")

saveRDS(BMDM_merged, file = "result/BMDM_merged.rds")

library(beyondcell, lib.loc = "/home/fengqian/miniconda3/envs/beyondcell/lib/R/library")

gs <- GetCollection(PSc)
#ss <- GetCollection(SSc)

sc <- BMDM_merged

bc <- bcScore(sc, gs, expr.thres = 0.1) 
#bc <- bcScore(ss, ss, expr.thres = 0.1)

# Run the UMAP reduction. 
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.2)

bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 0.5)

bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")
bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

bcRanks <- function(bc, idents = NULL, extended = TRUE, 
                    resm.cutoff = c(0.1, 0.9), 
                    sp.cutoff = c(0.1, 0.4, 0.6, 0.9)) {
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (idents %in% colnames(bc@meta.data)) {
      if (idents %in% names(bc@ranks)) {
        warning(paste0('$', idents, ' already exists in bc@ranks. ',
                       'Entry has been overwritten.'))
      }
      meta <- bc@meta.data[colnames(bc@normalized), idents, drop = TRUE]
    } else {
      stop('Idents not found.')
    }
  } else {
    stop("You must supply the name of a metadata column to group by.")
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # Check resm.cutoff.
  if (length(resm.cutoff) != 2 | !is.numeric(resm.cutoff)) {
    stop('resm.cutoff must be a numeric vector of length 2.')
  }
  if (resm.cutoff[2] < resm.cutoff[1]) {
    warning(paste('Upper residuals\' mean cut-off is smaller than lower', 
                  'residuals\' mean cut-off. Sorting residuals\' mean cut-offs', 
                  'in increasing order.'))
    resm.cutoff <- sort(resm.cutoff, decreasing = FALSE)
  }
  # Check sp.cutoff.
  if (length(sp.cutoff) != 4 | !is.numeric(sp.cutoff)) {
    stop('sp.cutoff must be a numeric vector of length 4.')
  }
  if (any(sp.cutoff < 0 | sp.cutoff > 1)) {
    stop('sp.cutoff must contain 4 switch point values between 0 and 1.')
  }
  sorted.sp.cutoff <- sort(sp.cutoff, decreasing = FALSE)
  if (!identical(sp.cutoff, sorted.sp.cutoff)) {
    warning(paste('Sorting switch point cut-offs in increasing order.'))
    sp.cutoff <- sorted.sp.cutoff
  }
  # --- Code ---
  # Progress bar.
  pb <- txtProgressBar(min = 0, max = 100, style = 3, file = stderr())
  bins <- 10
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # Keep column to group by.
  meta <- bc@meta.data %>%
    tibble::rownames_to_column("cells") %>%
    dplyr::select(cells, all_of(idents)) %>%
    dplyr::rename(group.var := !!idents) %>%
    dplyr::mutate(group.var = factor(group.var)) %>%
    unique()
  lvls <- levels(meta$group.var)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 5)
  # Column to order by.
  order.col <- paste0("rank.", levels(meta$group.var)[1])
  # Final column order.
  if (extended) {
    cols.additional <- c("median", "sd", "variance", "min", "max", "prop.na")
  } else {
    cols.additional <- NULL
  }
  cols.stats <- c("rank", "switch.point", "mean", cols.additional, 
                  "residuals.mean", "group")
  cols.stats.level <- tidyr::expand_grid(lvls, cols.stats) %>%
    dplyr::mutate(col.name = paste(cols.stats, lvls, sep = ".")) %>%
    dplyr::pull(col.name)
  # Get switch points.
  sp <- data.frame(switch.point = bc@switch.point) %>%
    tibble::rownames_to_column("IDs")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 10)
  # Compute long normalized BCS.
  normalized.long <- bc@normalized %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cells") %>%
    tidyr::pivot_longer(cols = all_of(sigs), names_to = "IDs", 
                        values_to = "enrichment", values_drop_na = FALSE)
  # Add grouping information and switch point.
  normalized.long <- normalized.long %>%
    dplyr::inner_join(sp, by = "IDs") %>%
    dplyr::inner_join(meta, by = "cells")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 25)
  # Compute mean BCS and residual's mean per signature.
  stats.long <- normalized.long %>%
    dplyr::group_by(IDs) %>%
    dplyr::mutate(mean = round(mean(enrichment, na.rm = TRUE), digits = 2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(resid = enrichment - mean) %>%
    dplyr::group_by(IDs, group.var) %>%
    dplyr::mutate(residuals.mean = round(mean(resid, na.rm = TRUE), digits = 2)) %>%
    dplyr::ungroup()
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 45)
  # If extended == TRUE, compute the median, standard deviation, variance, min, 
  # max and proportion of NaNs per signature.
  if (extended) {
    stats.long <- stats.long %>%
      dplyr::group_by(IDs) %>%
      dplyr::mutate(median = round(median(enrichment, na.rm = TRUE), digits = 2),
                    sd = round(sd(enrichment, na.rm = TRUE), digits = 2),
                    variance = round(var(enrichment, na.rm = TRUE), digits = 2),
                    min = round(min(enrichment, na.rm = TRUE), digits = 2),
                    max = round(max(enrichment, na.rm = TRUE), digits = 2),
                    prop.na = round(sum(is.na(enrichment))/length(cells), 
                                    digits = 2)) %>%
      dplyr::ungroup()
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 50)
  # Residual's deciles.
  res.decil <- stats.long %>%
    dplyr::group_by(group.var) %>%
    dplyr::group_modify(~as.data.frame(t(quantile(.$residuals.mean, resm.cutoff, 
                                                  na.rm = TRUE)))) %>%
    dplyr::ungroup()
  colnames(res.decil)[2:3] <- c("Pmin", "Pmax")
  stats.long <- stats.long %>%
    dplyr::select(-cells, -enrichment) %>%
    unique() %>%
    dplyr::inner_join(res.decil, by = "group.var")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 75)
  # Group annotation.
  stats.long.annotated <- stats.long %>%
    dplyr::mutate(
      group = dplyr::case_when(switch.point < sp.cutoff[1] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-HighSensitivity",
                               switch.point > sp.cutoff[4] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-LowSensitivity",
                               switch.point > sp.cutoff[2] & 
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-Differential-LowSensitivity",
                               switch.point > sp.cutoff[2] &
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-Differential-HighSensitivity",
                               TRUE ~ NA_character_))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 80)
  # Order.
  rank <- stats.long.annotated %>%
    dplyr::mutate(in.range = switch.point > sp.cutoff[2] & 
                    switch.point < sp.cutoff[3],
                  sp.rank = switch.point * as.numeric(in.range)) %>%
    dplyr::select(IDs, group.var, sp.rank, residuals.mean, in.range) %>%
    unique() %>%
    dplyr::group_split(group.var)
  rank <- lapply(rank, FUN = function(x) {
    dt <- data.table::as.data.table(x)
    dt[, rank := data.table::frank(dt, -sp.rank, -residuals.mean, 
                                   ties.method = "dense")]
    return(dt)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(rank = dplyr::if_else(in.range, rank, NA_integer_)) %>%
    dplyr::select(IDs, group.var, rank) %>%
    unique()
  stats.long.ranked <- stats.long.annotated %>%
    dplyr::inner_join(rank, by = c("IDs", "group.var"))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 85)
  # Pivot wider
  final.stats <- stats.long.ranked %>%
    dplyr::select(IDs, group.var, all_of(cols.stats)) %>%
    unique() %>%
    tidyr::pivot_wider(names_from = group.var, values_from = all_of(cols.stats),
                       names_sep = ".")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 90)
  # Add Drug name and MoA to final.stats.
  info <- drugInfo$IDs %>%
    dplyr::filter(IDs %in% final.stats$IDs) %>%
    dplyr::select(IDs, preferred.drug.names, studies) %>%
    dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Targets, by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Synonyms, by = "IDs",
                     relationship = "many-to-many")
  if (dim(info)[1] > 0) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(x) {
      paste(na.omit(unique(x)), collapse = "; ")
    })
    cols.druginfo <- c("drugs", "preferred.drug.names", "MoAs", "targets", 
                       "studies")
  } else {
    info <- data.frame(IDs = rownames(bc@normalized))
    cols.druginfo <- NULL
  }
  final.stats <- final.stats %>%
    dplyr::left_join(info, by = "IDs") %>%
    tibble::column_to_rownames("IDs") %>%
    unique()
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 95)
  # Order by rank and reorder columns.
  final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                             c(cols.druginfo, cols.stats.level)]
  # Add to beyondcell object.
  bc@ranks[[idents]] <- final.stats
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 100)
  return(bc)
}

bc <- bcRanks(bc, idents = "subcluster")
saveRDS(bc, file = "result/drug_prediction_bc_multifocal.rds")

# Explore the statistics table.
df_cor <- cor(t(bc@data), bc@meta.data$SSS)
colnames(df_cor) <- "cor"
df_cor <- df_cor %>% 
  as.data.frame() %>% 
  rownames_to_column("drugs") %>% 
  arrange(desc(cor))
head(df_cor)

# Step 1: Extract data from bc@ranks
df <- bc@ranks$subcluster

df <- df %>%
  dplyr::select(
    residuals.mean = residuals.mean.S2,
    switch.point = switch.point.S2,
    group = group.S2,
    drug = preferred.drug.names
  ) %>%
  tibble::rownames_to_column("sig_id") %>%
  dplyr::mutate(
    label = paste0(drug, " (", sig_id, ")")
  )

# Step 2: Label only top 5 from each group
df <- df %>%
  arrange(group, desc(abs(residuals.mean))) %>%
  group_by(group) %>%
  mutate(label_show = ifelse(
    (group == "TOP-Differential-HighSensitivity" & row_number() <= 10) | 
      (group == "TOP-Differential-LowSensitivity" & row_number() <= 3),
    label, 
    NA
  )) %>%
  ungroup()

# Step 3: Define color map
color_map <- c(
  "TOP-Differential-HighSensitivity" = "#F4A300",
  "TOP-Differential-LowSensitivity" = "#B086CC"
)

# Step 4: Set cutoffs
x_cutoff <- quantile(df$residuals.mean, probs = c(0.1, 0.9), na.rm = TRUE)
y_cutoff <- c(0.1, 0.4, 0.6, 0.9)

# Step 5: Plot
pdf("figures/Fig5I.pdf", height = 5, width = 8)
ggplot(df, aes(x = residuals.mean, y = switch.point)) +
  geom_point(aes(fill = group, color = group), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color_map, na.value = "lightgrey") +
  scale_color_manual(values = color_map, na.value = "lightgrey", guide = "none") +
  geom_vline(xintercept = x_cutoff, linetype = "dotted") +
  geom_hline(yintercept = y_cutoff, linetype = "dotted") +
  geom_text_repel(
    aes(label = label_show),
    size = 3,
    color = "black",
    nudge_x = ifelse(df$residuals.mean >= 0, 0.5, -0.5),
    direction = "y",      
    hjust = 0.5,           
    segment.size = 0.5,   
    segment.alpha = 0.6,  
    min.segment.length = 0, 
    box.padding = 0.5,  
    point.padding = 0.3,  
    force = 2,       
    max.time = 2,      
    max.iter = 20000,  
    max.overlaps = Inf   
  ) +
  theme_classic() +
  labs(
    title = "subcluster = High SSS",
    x = "Residuals' Mean",
    y = "Switch Point",
    fill = "Drug Group",
    caption = paste0("x cut-offs: first and last deciles; y cut-offs: ", paste(y_cutoff, collapse = ", "))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
dev.off()

library(patchwork)
FindDrugs(bc, "CEPHAELINE")
FindDrugs(bc, "EMETINE")
FindDrugs(bc, "CERAMIDE")

# UMAP with time points.
bcClusters(bc, UMAP = "beyondcell", idents = "label_SSS", pt.size = 0.5)

bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig-17814"), pt.size = 1.5)

pdf("figures/Sup_Fig5E.pdf")
bcHistogram(bc, signatures = "sig-7064", idents = "subcluster")
bcHistogram(bc, signatures = "sig-7353", idents = "label_SSS")
bcHistogram(bc, signatures = "sig-20501", idents = "label_SSS")
dev.off()

interest_mat <- df %>% 
  arrange(desc(residuals.mean)) %>% 
  filter(group=="TOP-Differential-HighSensitivity" & 
           switch.point < 0.6 & switch.point > 0.4 & 
           residuals.mean > x_cutoff[2] & sig_id %in% df_cor$drugs[df_cor$cor > 0.3]) %>% 
  dplyr::select(sig_id, drug)
View(interest_mat)

GBM_intergrated_24_TAM$subcluster_named <- GBM_intergrated_24_TAM$predicted.subcluster
GBM_intergrated_24_TAM$subcluster_named[GBM_intergrated_24_TAM$predicted.subcluster=="S0"] <- "FRMD4A+ MAC1"
GBM_intergrated_24_TAM$subcluster_named[GBM_intergrated_24_TAM$predicted.subcluster=="S3"] <- "HLA-hi MAC2"
GBM_intergrated_24_TAM$subcluster_named[GBM_intergrated_24_TAM$predicted.subcluster=="S2"] <- "NUPR1+ MDSC"
GBM_intergrated_24_TAM$subcluster_named[GBM_intergrated_24_TAM$predicted.subcluster=="S1"] <- "VCAN+ MDSC"
GBM_intergrated_24_TAM$subcluster_named[GBM_intergrated_24_TAM$predicted.subcluster=="S4"] <- "CXCL8+ MDSC"

write.table(clinical_data_TCGA[,c("patient","sample","SSS","group","OS_time", "OS_status")], file = "tables/Fig4B.tsv", quote = F, sep = ",")
write.table(clinical_data, file = "tables/Fig4C.tsv", quote = F, sep = ",")
write.table(merged_data_cibersort[,c("patient","sample","SSS","group","OS_time", "OS_status","S2_cibersort_group")], file = "tables/Fig4D.tsv", quote = F, sep = ",")
write.table(final_data, file = "tables/Fig4E.tsv", quote = F, sep = ",")
write.table(GBM_intergrated_24_TAM@meta.data[,c("diagnosis","PD1","treatment_3","subcluster_named")], file = "tables/Fig4F.tsv", quote = F, sep = "\t")
write.table(plot_data, file = "tables/Fig4G.tsv", quote = F, sep = "\t")
write.table(df, file = "tables/Fig4H.tsv", quote = F, sep = "\t")
write.table(vlnplot$data, file = "tables/Fig4I.tsv", quote = F, sep = "\t")


