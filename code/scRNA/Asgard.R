library('Asgard')

#Please replace Your_local_path with your real local folder

PrepareReference(cell.info="reference/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="reference/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "reference/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "reference/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "reference/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "reference/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "DrugReference/"
)

#Note: the file names here maybe different after unzipping.
#Please note that it takes more than one hour to produce drug references in a standard computer with RAM>64GB.

library('Seurat')

# Load cells' cell type annotations for GSE113197
cell_types_file <- paste0(
  "https://raw.githubusercontent.com/lanagarmire/"
  "Single-cell-drug-repositioning/master/Drug/Normal_celltype.txt"
)
cell_types <- read.table(file=celltypes, header=TRUE, check.names=FALSE)

# Cell type of interest
cell_types_names <- c(
  "Luminal_L2_epithelial_cells", "Luminal_L1.1_epithelial_cells", 
  "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"
)

# Load normal sample Ind5 from GSE113197 dataset 
data <- read.table(file="GSM3099847_Ind5_Expression_Matrix.txt", 
                   header=TRUE, check.names=FALSE)
row.names(data) <- data[, 1]
data <- data[, -1]
ind5_cells <- subset(cell_type, sample=="Ind5" & celltype %in% celltypes_names)
common <- intersect(colnames(data), rownames(ind5_cells))
data <- data[, common]

metadata = data.frame(
  ind5_celltypes,
  cell = colnames(data),
  type = "normal"
)
epithelial2 <- CreateSeuratObject(counts=data, project="Epithelial", min.cells=3, 
                                  min.features=200, meta.data=metada)

#Load normal sample Ind6 from GSE113197 dataset
data <- read.table(file="GSM3099848_Ind6_Expression_Matrix.txt", header=TRUE,
                   check.names=FALSE)
row.names(data) <- data[, 1]
data <- data[, -1]
ind6_cells <- subset(celltype,sample=="Ind6" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype3))
data<-data[,common]
Epithelial3 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype3,cell=colnames(data),type="Normal"))

#Load normal sample Ind7 from GSE113197 dataset
data<-read.table(file="GSM3099849_Ind7_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype4<-subset(celltype,sample=="Ind7" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype4))
data<-data[,common]
Epithelial4 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype4,cell=colnames(data),type="Normal"))

#Load cancer sample PDX110 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516947_PDX110")
TNBC.PDX2 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-110",type="TNBC.PDX"))

#Load cancer sample PDX322 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516948_PDX322")
TNBC.PDX3 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-332",type="TNBC.PDX"))




