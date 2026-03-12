## Figure 4A: calculation for the Venn plot--------
library(Seurat)
library(ArchR)
library(dplyr)
library(VennDiagram)


marker_list <- readRDS("result/BMDM_markers_filter.rds")
BMDM_subcluster_markers <- marker_list$S2
length(BMDM_subcluster_markers) #46

## bulk RNA-seq
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
    gene_name %in% intersection ~ "SMF gene",
    gene_name %in% interest_tf ~ "TF",
    gene_name %in% c(score_list$Ferroptosis_neg, score_list$Ferroptosis_dual) ~ "Ferroptosis_related",
    !gene_name %in% c(intersection, interest_tf, score_list$Ferroptosis_neg, score_list$Ferroptosis_dual) ~ "Other"
  ))

erastin_induced_genes <- res %>%
  filter(FDR < 0.05 & logFC > 0) %>%
  pull(gene_name)

intersection_refined <- intersect(BMDM_subcluster_markers, erastin_induced_genes)
write.table(intersection_refined, file = "tables/TableS5_refined.tsv", quote = F, row.names = F)

library(ggplot2)
library(ggrepel)

plot_data <- res %>%
  mutate(
    is_S2_marker = gene_name %in% BMDM_subcluster_markers,
    is_Selected = gene_name %in% intersection_refined,
    
    color_group = case_when(
      is_Selected ~ "SMF Functional Core", 
      is_S2_marker & !is_Selected ~ "S2 Repressed", 
      TRUE ~ "Background"
    ),
    
    label = ifelse(is_Selected, gene_name, NA)
  )

xlim_val <- 2.5
ylim_val <- 16

pdf("figures/Fig4A_review.pdf", height = 4, width = 5)
ggplot(plot_data, aes(x = logFC, y = -log10(FDR))) +
  
  # 1. Background Genes (Light Grey)
  # Keeping them faint to highlight the S2 markers
  geom_point(data = subset(plot_data, color_group == "Background"), 
             color = "grey90", size = 1, alpha = 0.5) +
  
  # 2. S2 Repressed Genes (Scientific Blue)
  # These are S2 markers that are downregulated by Erastin (Left side)
  geom_point(data = subset(plot_data, color_group == "S2 Repressed"), 
             color = "#3C5488", size = 2, alpha = 0.8) +
  
  # 3. SMF Functional Core (Scientific Red)
  # These are the adaptive S2 markers upregulated by Erastin (Right side)
  geom_point(data = subset(plot_data, color_group == "SMF Functional Core"), 
             color = "#E64B35", size = 3) +
  
  # 4. Threshold Lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = c(0), linetype = "dashed", color = "grey60") +
  
  # 5. Labels (Only for the Functional Core)
  geom_text_repel(aes(label = label), 
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  max.overlaps = 50,
                  size = 3.5, 
                  fontface = "bold", 
                  color = "#E64B35",
                  # Ensure labels stay within the zoomed area
                  xlim = c(-xlim_val, xlim_val),
                  ylim = c(0, ylim_val)) +
  
  # 6. Theme and Aesthetics
  theme_classic() +
  labs(x = "Log2 Fold Change (Erastin vs Control)",
       y = "-Log10 FDR") +
  
  # 7. Zoom in without removing data points
  coord_cartesian(xlim = c(-xlim_val, xlim_val), 
                  ylim = c(0, ylim_val))
dev.off()

# pdf("figures/Fig4A_review.pdf", height = 4, width = 4)
# grid.newpage()
# venn.plot <- draw.pairwise.venn(
#   area1 = length(BMDM_subcluster_markers),
#   area2 = length(erastin_induced_genes),
#   cross.area = length(intersection_refined),
#   category = c("S2 BMDM markers", "PC vs IT upregulated markers"),
#   fill = c("coral", "lightblue"),
#   cat.pos = c(0, 0),        # Positions in degrees (0 is center-top)
#   cat.dist = c(0.02, 0.02), # Smaller distance = closer to center
#   cat.cex = 1.2             # Optional: change label size
# )
# dev.off()
    
## Figure 4B: Survival Analysis with Expression-----------
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
gene_set_list <- list(SMF = as.character(intersection_refined))

# Run ssGSEA
ssgsea_scores <- gsva(fpkm_data_TCGA, gene_set_list, method = "ssgsea", verbose = FALSE)

# Add to clinical data
clinical_data_TCGA <- as.data.frame(colData(expression_data)) 
clinical_data_TCGA <- clinical_data_TCGA %>% 
  mutate(SMF = ssgsea_scores["SMF", clinical_data_TCGA$barcode])

# Group by quartiles
clinical_data_TCGA <- clinical_data_TCGA %>%
  mutate(group = case_when(
    SMF <= median(SMF) ~ "Low",
    SMF >= median(SMF) ~ "High",
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
pdf("figures/Fig4B_review.pdf", height = 5, width = 4.5)
ggsurvplot(fit, data = clinical_data_TCGA,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "Survival Analysis with SMF ssGSEA scoring in TCGA GBM",
           palette = c("#8B235A","#005787"))
dev.off()

## Figure 4C: Survival Analysis with CGGA-------------
fpkm_data1 <- read_table("reference/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
fpkm_data2 <- read_table("reference/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")

# Keep shared genes and merge expression matrices
identical_genes <- intersect(fpkm_data1$Gene_Name, fpkm_data2$Gene_Name)
fpkm_data1 <- fpkm_data1 %>% column_to_rownames("Gene_Name")
fpkm_data2 <- fpkm_data2 %>% column_to_rownames("Gene_Name")
fpkm_data <- cbind(fpkm_data1[identical_genes, ], fpkm_data2[identical_genes, ])
saveRDS(fpkm_data, file = "reference/CGGA/CGGA_GBM_expression_data.rds")

# Define your signature gene set
gene_set_list <- list(SMF = as.character(intersection_refined))

# Run ssGSEA
ssgsea_scores <- gsva(as.matrix(fpkm_data), gene_set_list, method = "ssgsea", verbose = FALSE)

# Load and merge clinical data
clinical_data1 <- read.table("reference/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data2 <- read.table("reference/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data_CGGA <- rbind(clinical_data1, clinical_data2) %>%
  dplyr::filter(Histology == "GBM")

# Match ssGSEA score to clinical samples
clinical_data_CGGA$signature_score <- ssgsea_scores["SMF", clinical_data_CGGA$CGGA_ID]

# Group patients by quartiles of signature score
clinical_data_CGGA <- clinical_data_CGGA %>%
  mutate(group = case_when(
    signature_score < median(signature_score) ~ "Low",
    signature_score >= median(signature_score) ~ "High",
    TRUE ~ NA_character_
  ))

# Prepare survival variables
surv_obj <- Surv(time = clinical_data_CGGA$OS, event = clinical_data_CGGA$Censor..alive.0..dead.1.)
fit <- survfit(surv_obj ~ group, data = clinical_data_CGGA)

# Plot KM survival curves
pdf("figures/Fig4C_review.pdf", height = 5, width = 4.5)
ggsurvplot(fit, data = clinical_data_CGGA,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "Survival Analysis with SMF ssGSEA scoring in CGGA GBM",
           palette = c("#FDDE00", "#1068B1"))
dev.off()

## Figure 4D: Correlation of immune score----------
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

pdf("figures/Fig4D_review.pdf", height = 4, width = 4)
gseaplot2(ego, geneSetID = ego@result$ID[ego@result$Description == "negative regulation of immune response"], 
          title = "negative regulation of immune response")
dev.off()

## Figure 4D: Survival Analysis with CIBERSORTx Deconvolution---------
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

# Define high and low NUPR1+ BMDM groups
merged_data_cibersort <- merged_data_cibersort %>%
  mutate(S2_cibersort_group = case_when(
    S2 <= median(S2) ~ "Low",
    S2 >= median(S2) ~ "High",
    TRUE ~ NA_character_ 
  ))

# Kaplan-Meier survival analysis
surv_obj_cibersort <- Surv(time = merged_data_cibersort$OS_time, event = merged_data_cibersort$OS_status)
fit_cibersort <- survfit(surv_obj_cibersort ~ S2_cibersort_group, data = merged_data_cibersort)

pdf("figures/Fig4E_review.pdf", height = 5, width = 4.5)
ggsurvplot(fit_cibersort, data = merged_data_cibersort,
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "NUPR1+ MDSC Proportion",
           legend.labs = c("High", "Low"),
           xlab = "Days",
           ylab = "OS Rate (%)",
           title = "NUPR1+ MDSC Proportion in TCGA GBM (deconvoluted)",
           palette = c("#D55E00","#0072B2"))
dev.off()

# ## Figure 4E: GLASS cohort for longitudinal validation---------------
# # Load expression and metadata
# gene_tpm_matrix <- read_tsv("reference/GLASS/gene_tpm_matrix_all_samples.tsv")
# analysis_rnaseq_pairs <- read_csv("reference/GLASS/analysis_rnaseq_pairs.csv") %>%
#   mutate(across(starts_with("tumor_barcode"), ~ gsub("-", ".", .)))
# 
# clinical_cases <- read_csv("reference/GLASS/clinical_cases.csv")
# clinical_surgeries <- read_csv("reference/GLASS/clinical_surgeries.csv")
# 
# # Merge clinical info
# clinical_GLASS <- clinical_cases %>%
#   inner_join(clinical_surgeries, by = "case_barcode") %>%
#   dplyr::filter(grade == "IV") %>%
#   dplyr::select(case_barcode, case_project, case_source, case_sex, case_age_diagnosis_years,
#          case_vital_status, case_overall_survival_mo, grade, idh_status, histology) %>%
#   distinct()
# 
# merged_data <- analysis_rnaseq_pairs %>%
#   inner_join(clinical_GLASS, by = "case_barcode")
# 
# # Ensure your gene set is a character vector of gene symbols
# gene_set_list <- list(signature = as.character(intersection_refined))
# 
# # Prepare expression matrix: rows = genes, columns = samples
# tpm_matrix <- gene_tpm_matrix %>%
#   column_to_rownames("Gene_symbol") %>%
#   as.matrix()
# 
# # Optional: log2 transform
# tpm_matrix_log <- log2(tpm_matrix + 1)
# 
# # Compute ssGSEA score per sample
# ssgsea_result <- gsva(tpm_matrix_log, gene_set_list, method = "ssgsea", verbose = FALSE)
# 
# # Reshape to long format
# signature_scores <- data.frame(
#   tumor_barcode = colnames(ssgsea_result),
#   signature = as.numeric(ssgsea_result[1, ])
# )
# 
# # Merge ssGSEA scores into paired data
# final_data <- merged_data %>%
#   inner_join(signature_scores, by = c("tumor_barcode_a" = "tumor_barcode")) %>%
#   rename(signature_a = signature) %>%
#   inner_join(signature_scores, by = c("tumor_barcode_b" = "tumor_barcode")) %>%
#   rename(signature_b = signature) 
# 
# quartile_a <- quantile(final_data$signature_a, probs = c(0.5), na.rm = TRUE)
# quartile_b <- quantile(final_data$signature_b, probs = c(0.5), na.rm = TRUE)
# 
# final_data <- final_data %>%
#   mutate(
#     tumor_type_a = ifelse(sample_type_a == "TP", "Primary", "Recurrent"),
#     tumor_type_b = ifelse(sample_type_b == "TP", "Primary", "Recurrent"),
#     survival_time = case_overall_survival_mo,
#     survival_status = ifelse(case_vital_status == "dead", 1, 0),
#     signature_group = case_when(
#       signature_a <= signature_b ~ "Up",
#       signature_a > signature_b ~ "Down",
#       TRUE ~ NA
#     ),
#     diff = signature_b - signature_a,
#     group_a = case_when(
#       signature_a <= quartile_a[1] ~ "Low",
#       signature_a >= quartile_a[1] ~ "High",
#       TRUE ~ NA
#     ),
#     group_b = case_when(
#       signature_b <= quartile_b[1] ~ "Low",
#       signature_b >= quartile_b[1] ~ "High",
#       TRUE ~ NA
#     ),
#     group_combined = paste0(group_a, "_to_", group_b)
#   ) %>%
#   arrange(diff) %>%
#   mutate(case_barcode_ordered = paste0(case_barcode, "_", row_number())) %>%
#   mutate(
#     case_sex = factor(case_sex),
#     idh_status = factor(idh_status)
#   )
# 
# # Survival analysis
# surv_obj <- Surv(time = final_data$survival_time, event = final_data$survival_status)
# fit <- survfit(surv_obj ~ group_a + group_b, data = final_data)
# 
# pdf("figures/Fig4E_review.pdf", height = 6, width = 6)
# ggsurvplot(fit, data = final_data,
#            pval = TRUE,
#            risk.table = TRUE,
#            legend.title = "Gene Set Change",
#            xlab = "Survival Time (Months)",
#            ylab = "Overall Survival Probability")
# dev.off()
# 
# # Waterfall plot
# pdf("figures/Sup_Fig4C.pdf", height = 4, width = 6)
# ggplot(final_data, aes(x = factor(case_barcode_ordered, levels = case_barcode_ordered), y = diff, fill = signature_group)) +
#   geom_bar(stat = "identity") +
#   labs(
#     x = "Case Barcode",
#     y = "Difference in Gene Set Score (Recurrent - Primary)",
#     fill = "Change Group",
#     title = "Waterfall Plot of Gene Set Score Differences"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
# dev.off()

## Figure 4F: cox model------
h_hypoxia <- score_list$Hypoxia

macrophage_sig <- score_list$BMDM

gene_set_list <- list(
  SMF = as.character(intersection_refined),
  Hypoxia_Score = h_hypoxia,
  Macrophage_Abundance = macrophage_sig
)

# Run ssGSEA
ssgsea_scores <- gsva(fpkm_data_TCGA, gene_set_list, method = "ssgsea", verbose = FALSE)

# Prepare Clinical Data
clinical_data_TCGA <- as.data.frame(colData(expression_data)) 

# Merge Scores into Clinical Data
clinical_data_TCGA <- clinical_data_TCGA %>% 
  mutate(
    SMF = ssgsea_scores["SMF", clinical_data_TCGA$barcode],
    Hypoxia_Score = ssgsea_scores["Hypoxia_Score", clinical_data_TCGA$barcode],
    Macrophage_Abundance = ssgsea_scores["Macrophage_Abundance", clinical_data_TCGA$barcode]
  )

cox_data_stratified <- clinical_data_TCGA %>%
  filter(paper_IDH.status == "WT") %>% 
  mutate(
    # Handle Survival Data
    OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS_status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS_time) & !is.na(OS_status)) %>% # Remove missing survival data
  mutate(
    Hypoxia_Group = ifelse(Hypoxia_Score > median(Hypoxia_Score, na.rm=TRUE), "High Hypoxia", "Low Hypoxia"),
    Hypoxia = scale(Hypoxia_Score),
    MGMT = factor(paper_MGMT.promoter.status),
    SMF = scale(SMF),
    SMF_Group = factor(ifelse(SMF > median(SMF, na.rm=TRUE), "High SMF", "Low SMF"), levels = c("Low SMF","High SMF")),
    BMDM = scale(Macrophage_Abundance),
    Age = scale(age_at_diagnosis))

fit_all <- coxph(Surv(OS_time, OS_status) ~ 
                    SMF_Group + 
                    Hypoxia +
                    BMDM +
                    Age +
                    MGMT, 
                  data = cox_data_stratified)
summary(fit_all)

pdf("figures/Fig4F_review.pdf", height = 4, width = 7)
forest_model(fit_all)
dev.off()

## Figure 4G: SMF correlates with MES----------
library(readxl)
markers<-as.list(read_excel("reference/mmc2.xlsx", skip = 3))
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1/S","G2/M")]
markers<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)
mes_genes <- markers$MES

gene_set_list_mes <- list(MES_Score = mes_genes, SMF_Score = intersection_refined)
scores <- gsva(fpkm_data_TCGA, gene_set_list_mes, method = "ssgsea", verbose = F)

df_cor <- as.data.frame(t(scores))
library(ggpubr)

pdf("figures/Fig4G_review.pdf", height = 4, width = 4)
ggscatter(df_cor, x = "SMF_Score", y = "MES_Score", 
          color = "grey70",    
          shape = 16,          
          size = 1.5,         
          alpha = 0.6,        
          
          add = "reg.line", 
          
          add.params = list(color = "#B22222", fill = "#CCCCCC"), 
          
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          
          xlab = "SMF Signature Score", 
          ylab = "Mesenchymal (MES) Score") +
  
  ggtitle("SMF correlates with MES state") +
  theme_classic()
dev.off()

## Figure 4H: immunotherapy response-----------
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

GBM_intergrated_24 <- readRDS("reference/GBM.RNA.integrated.24.rds")
GBM_intergrated_24_TAM <- subset(GBM_intergrated_24, subset = anno_ident == "Macrophages")

# Define standard workflow (Looks good)
seurat_workflow <- function(seu){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  # Note: Do not run UMAP/Neighbors here if you run Harmony immediately after
  seu
}

# Process Query
GBM_intergrated_24_TAM <- seurat_workflow(GBM_intergrated_24_TAM)

# Run Harmony to correct patient batch effects in the Query
GBM_intergrated_24_TAM <- GBM_intergrated_24_TAM %>% 
  RunHarmony("Pt_number", plot_convergence = FALSE, nclust = 50, max_iter = 10, early_stop = TRUE)

# Calculate Query's own UMAP (based on Harmony)
GBM_intergrated_24_TAM <- GBM_intergrated_24_TAM %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.2) %>% # Increased resolution slightly (0.1 might be too low for large cohorts)
  RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony") 

BMDM <- readRDS("result/BMDM_withXHsample_harmony.rds")
BMDM <- RunUMAP(BMDM, dims = 1:30, reduction = "pca", return.model = TRUE) 

# Find anchors
anchors <- FindTransferAnchors(
  reference = BMDM, 
  query = GBM_intergrated_24_TAM, 
  dims = 1:30, 
  reference.reduction = "pca"
)

GBM_intergrated_24_TAM <- MapQuery(
  anchorset = anchors,
  query = GBM_intergrated_24_TAM,
  reference = BMDM,
  refdata = list(subcluster = "subcluster"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

SubclusterColor <- c('#ea9994','#f2c396','#9cd2ed','#86c7b4','#a992c0')
names(SubclusterColor) <- unique(BMDM$subcluster)

pdf("figures/Sup_Fig4A_review.pdf", height = 5, width = 5)
DimPlot(GBM_intergrated_24_TAM, 
        reduction = "umap.harmony", 
        group.by = "predicted.subcluster", 
        cols = SubclusterColor,
        label = TRUE, label.size = 3) + 
  ggtitle("Predicted Labels on Query's Structure")
dev.off()

# Save the object
saveRDS(GBM_intergrated_24_TAM, file = "result/BMDM_neoPD1_treatments.rds")

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

pdf("figures/Sup_Fig4B_review.pdf", height = 3, width = 5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2)))
dev.off()

GBM_intergrated_24_TAM <- AddModuleScore_UCell(GBM_intergrated_24_TAM, list("SMF"=intersection_refined,
                                                                            "S2"=BMDM_subcluster_markers), name = NULL)

pdf("figures/Sup_Fig4C_review.pdf", height = 3, width = 6)
DotPlot(GBM_intergrated_24_TAM, 
                  features = c(intersection_refined, "SMF"), 
                  group.by = "treatment_3",
                  cols = c("lightgrey", "#E64B35"),
                  dot.scale = 6) + 
  RotatedAxis() +
  ggtitle("Refined SMF Signature Genes") +
  theme(axis.text.x = element_text(size = 9),
        axis.title = element_blank())
dev.off()

# 1. Calculate Average Expression
Idents(GBM_intergrated_24_TAM) <- "treatment_3"
avg_data <- AverageExpression(GBM_intergrated_24_TAM, 
                              features = rownames(GBM_intergrated_24_TAM),
                              assays = "RNA", 
                              slot = "data")$RNA %>% as.data.frame()

# 2. Prepare Data for Plotting
plot_data <- avg_exp %>% as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  mutate(
    x_val = log2(responder + 1),
    y_val = log2(nonresponder + 1),
    logFC = y_val - x_val,
    
    is_target = gene %in% intersection_refined,
    
    to_label = ifelse(is_target & logFC > 0, gene, NA) 
  )

cor_plot <- ggplot(plot_data, aes(x = x_val, y = y_val)) +
  
  geom_point(data = subset(plot_data, !is_target), 
             color = "grey90", size = 1, alpha = 0.3) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey60") +
  
  geom_point(data = subset(plot_data, is_target), 
             aes(fill = logFC), 
             shape = 21,        
             size = 3, 
             color = "black",   
             stroke = 0.3) +
  
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#E64B35", 
                       midpoint = 0, limits = c(-1, 1), oob = scales::squish) +
  
  geom_text_repel(aes(label = to_label), 
                  size = 3, 
                  fontface = "bold",
                  color = "#E64B35",
                  box.padding = 0.5, 
                  max.overlaps = 50,
                  segment.color = "grey50",
                  min.segment.length = 0) +
  

  labs(
    x = "Expression in Responders (log2 CPM)",
    y = "Expression in Non-Responders (log2 CPM)",
    title = "Refined SMF Signature",
    subtitle = "Color indicates bias towards Non-Responders"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "right",
    aspect.ratio = 1 
  ) +
  coord_cartesian(xlim = c(-0.1, 6), 
                  ylim = c(-0.1, 6))

pdf("figures/Sup_Fig4D_review.pdf", height = 5, width = 6)
cor_plot
dev.off()

# Figure 4K: drug response in SMF-high BMDM population--------------
BMDM_merged <- readRDS("result/BMDM_merged.rds")
BMDM_merged$subcluster[BMDM_merged$TAM_subtype=="other TAMs"] <- "other"
BMDM_merged$subcluster[BMDM_merged$TAM_subtype=="S2 BMDM"] <- "S2"
BMDM_merged$subcluster[BMDM_merged$subcluster!="S2"] <- "other"
BMDM_merged$subcluster[BMDM_merged$predicted.subcluster!="S2"] <- "other"
BMDM_merged$subcluster[BMDM_merged$predicted.subcluster=="S2"] <- "S2"
table(BMDM_merged$subcluster)

intersection <- read.table("tables/TableS5_refined.tsv", header = T)[,1]

DimPlot(BMDM_merged, group.by = "subcluster", reduction = "tsne")

library(UCell)
BMDM_merged <-AddModuleScore_UCell(BMDM_merged, list(SMF=intersection), name = NULL)
FeaturePlot(BMDM_merged, "SMF", max.cutoff = "q90", min.cutoff = "q10", reduction = "tsne")

BMDM$label_SMF <- "Middle"
sample_list <- unique(BMDM$orig.ident)  

for (sample_id in sample_list) {
  cells_in_sample <- colnames(BMDM)[BMDM$orig.ident==sample_id]
  SMF_scores <- BMDM$SMF[cells_in_sample]
  
  low_cutoff <- quantile(SMF_scores, 1/3, na.rm = TRUE)
  high_cutoff <- quantile(SMF_scores, 2/3, na.rm = TRUE)
  
  BMDM$label_SMF[cells_in_sample][SMF_scores < low_cutoff] <- "Low"
  BMDM$label_SMF[cells_in_sample][SMF_scores >= high_cutoff] <- "High"
}

DimPlot(BMDM, group.by="label_SMF")

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
df_cor <- cor(t(bc@data), bc@meta.data$SMF)
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
pdf("figures/Fig4I.pdf", height = 5, width = 8)
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
    title = "subcluster = High SMF",
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
bcClusters(bc, UMAP = "beyondcell", idents = "label_SMF", pt.size = 0.5)

bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig-17814"), pt.size = 1.5)

pdf("figures/Sup_Fig5D.pdf")
bcHistogram(bc, signatures = "sig-7064", idents = "subcluster")
bcHistogram(bc, signatures = "sig-7353", idents = "label_SMF")
bcHistogram(bc, signatures = "sig-20501", idents = "label_SMF")
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

write.table(clinical_data_TCGA[,c("patient","sample","SMF","group","OS_time", "OS_status")], file = "tables/Fig4B.tsv", quote = F, sep = ",")
write.table(clinical_data, file = "tables/Fig4C.tsv", quote = F, sep = ",")
write.table(merged_data_cibersort[,c("patient","sample","SMF","group","OS_time", "OS_status","S2_cibersort_group")], file = "tables/Fig4D.tsv", quote = F, sep = ",")
write.table(final_data, file = "tables/Fig4E.tsv", quote = F, sep = ",")
write.table(GBM_intergrated_24_TAM@meta.data[,c("diagnosis","PD1","treatment_3","subcluster_named")], file = "tables/Fig4F.tsv", quote = F, sep = "\t")
write.table(plot_data, file = "tables/Fig4G.tsv", quote = F, sep = "\t")
write.table(df, file = "tables/Fig4H.tsv", quote = F, sep = "\t")
write.table(vlnplot$data, file = "tables/Fig4I.tsv", quote = F, sep = "\t")


