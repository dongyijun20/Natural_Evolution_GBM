library(Seurat)

malignant <- readRDS("result/malignant_new.rds")
myeloid <- readRDS("result/myeloid.rds")

counts <- t(as.data.frame(malignant@assays$RNA@counts))
write.table(counts[malignant$sample=="NJ01_1",], file = "cNMF/count_data/NJ01_1.count.txt", quote = F)
write.table(counts[malignant$sample=="NJ01_2",], file = "cNMF/count_data/NJ01_2.count.txt", quote = F)
write.table(counts[malignant$sample=="NJ02_1",], file = "cNMF/count_data/NJ02_1.count.txt", quote = F)
write.table(counts[malignant$sample=="NJ02_2",], file = "cNMF/count_data/NJ02_2.count.txt", quote = F)
write.table(counts[malignant$sample=="TT01_1",], file = "cNMF/count_data/TT01_1.count.txt", quote = F)
write.table(counts[malignant$sample=="TT01_2",], file = "cNMF/count_data/TT01_2.count.txt", quote = F)
write.table(counts[malignant$sample=="TT02_1",], file = "cNMF/count_data/TT02_1.count.txt", quote = F)
write.table(counts[malignant$sample=="TT02_2",], file = "cNMF/count_data/TT02_2.count.txt", quote = F)

library(reticulate)

use_condaenv(condaenv = "cnmf_env", required = T,conda = "/home/hsy/miniconda3/bin/conda")
py_config() #如果显示cnmf_env环境里面的python就OK

source("1.R")
step1(dir_input = "count_data",dir_output = "res1",k=3:5,iteration = 50) #这里为了演示方便，取值都比较小

source("2.R")
step2(dir_input = "res1",dir_output = "res2",dir_count = "count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)



