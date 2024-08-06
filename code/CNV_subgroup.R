library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(scales)

mat=read.table("./infercnv/older_malignant/infercnv.observations.txt",header=T,row.names=1,sep=" ",stringsAsFactors=F)
colnames(mat)=str_replace(colnames(mat),"^X","")
colnames(mat)=str_replace(colnames(mat),"\\.","-")

forCNVanno=read.table("./uphyloplot2-2.3/Inputs/older.cell_groupings",header = T,sep = "\t",stringsAsFactors = F)
forCNVanno=forCNVanno[!str_detect(forCNVanno$cell_group_name,"non-malignant"),]
forCNVanno$cell_group_name=str_replace(forCNVanno$cell_group_name,"malignant\\.malignant\\.","")
colnames(forCNVanno)=c("subgroup","CB")
rownames(forCNVanno)=forCNVanno$CB
forCNVanno$CB=NULL
forCNVanno$subgroup=factor(forCNVanno$subgroup,levels = sort(unique(forCNVanno$subgroup)))
mat=mat[,rownames(forCNVanno)]

gn <- rownames(mat)
geneFile <- read.table("geneLocate.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(geneFile$V1,gn),]
mat=mat[intersect(geneFile$V1,gn),]

##定义注释
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
len_c=length(names(table(forCNVanno$subgroup)))
color_c=RColorBrewer::brewer.pal(12, "Paired")[5:(len_c+4)]
names(color_c)=names(table(forCNVanno$subgroup))
left_anno <- rowAnnotation(df = forCNVanno,col=list(subgroup=color_c),border = F)

#重要的癌基因，随便找了一些，不同癌症不一样
#cosmic=read.table("Census_allMon.small.tsv",header = F,sep = "\t",stringsAsFactors = F)
# subset <- subset(pbmc, downsample=3000)
# Idents(subset) <- "tumor_ident"
# markers <- FindMarkers(subset, ident.1 = "malignant", ident.2 = "non-malignant")
# markers_df = markers[markers$avg_log2FC>0,] %>% top_n(n = 50, wt = avg_log2FC) %>% arrange(avg_log2FC)
# write.table(markers_df,file="result/malignant_markers50.txt",quote=F,sep="\t",row.names=T,col.names=T)

markers_df <- read.table("result/malignant_markers50.txt", header = T)
tail(markers_df,20)

# pdf("plot_new/malignant_featureplot.pdf",width = 12,height = 10)
# FeaturePlot(pbmc, features=c("SOX2","SOX4","PTPRZ1","PTN","CLU","CRYAB","TUBB2B","S100B","TUBA1A"), min.cutoff = "q1", max.cutoff = "q99",cols = c("lightgrey" ,"#DE1F1F"))
# dev.off()

#发生了CNV的基因
pre_stat=read.table("./infercnv/older_malignant/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_genes.dat",header = T,sep = "\t",stringsAsFactors = F)
pre_stat$cell_group_name=str_replace(pre_stat$cell_group_name,"all.*observations\\.","")
pre_stat=pre_stat %>% filter(state != 3)
####同一基因 不同CNV事件 要区别开来
pre_stat$cnv_type=""
chr_pq=read.table("reference/chr_pq.new.txt",header = T,sep = "\t",stringsAsFactors = F) 
#该文件由band文件转换得到，参考我的简书文章：染色体区段6p21.31在哪？https://www.jianshu.com/p/477a07192fe6
chr_pq$chr=paste("chr",chr_pq$chr,sep = "")

for (i in 1:nrow(pre_stat)) {
  tmp_chr_pq=chr_pq%>%filter(chr==pre_stat[i,"chr"])
  if(pre_stat[i,"start"] <= tmp_chr_pq[1,"cutoff"] & pre_stat[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
    pre_stat[i,"cnv_type"]=paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
  }else if (pre_stat[i,"start"] >= tmp_chr_pq[2,"cutoff"] & pre_stat[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
    pre_stat[i,"cnv_type"]=paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
  }else {
    pre_stat[i,"cnv_type"]=paste(pre_stat[i,"chr"],"p,",pre_stat[i,"chr"],"q",sep = "")
  }
  if (pre_stat[i,"state"] < 3) {
    pre_stat[i,"cnv_type"]=paste(pre_stat[i,"cnv_type"],"_loss",sep = "")
  }
  if (pre_stat[i,"state"] > 3) {
    pre_stat[i,"cnv_type"]=paste(pre_stat[i,"cnv_type"],"_gain",sep = "")
  }
  if (str_detect(pre_stat[i,"cnv_type"],",")) {
    tmp1=strsplit(pre_stat[i,"cnv_type"],",")[[1]][1]
    tmp2=strsplit(strsplit(pre_stat[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
    tmp3=strsplit(strsplit(pre_stat[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
    pre_stat[i,"cnv_type"]=paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
    rm(list = c("tmp1","tmp2","tmp3"))
  }
}
pre_stat$gene_cnv_type=paste(pre_stat$gene,pre_stat$cnv_type,sep = "_")
pre_stat$cell_group_name=str_replace(pre_stat$cell_group_name,"malignant\\.malignant\\.","")

##发生在极少数细胞的CNV事件——不考虑
tmpdf=as.data.frame(table(forCNVanno$subgroup))
colnames(tmpdf)=c("cell_group_name","cellnum")
pre_stat=pre_stat%>%inner_join(tmpdf,by="cell_group_name")
tmpdf2=aggregate(pre_stat$cellnum,by=list(pre_stat$gene_cnv_type),FUN=sum)
colnames(tmpdf2)=c("gene_cnv_type","num")
tmpdf2=tmpdf2%>%filter(num >= dim(mat)[2] * 0.05)
pre_stat=pre_stat[pre_stat$gene_cnv_type %in% tmpdf2$gene_cnv_type,]
#交集
key_gene=intersect(rownames(mat),intersect(pre_stat$gene,rownames(markers_df)))
ha=columnAnnotation(foo=anno_mark(at=which(rownames(mat) %in% key_gene),labels = rownames(mat)[which(rownames(mat) %in% key_gene)],
                                  which = "column",side = "bottom"))

pdf("plot_new/CNV_heatmap.clone_older.pdf",width = 15,height = 11)
ht_tp = Heatmap(t(mat),
                col = colorRamp2(c(0.6,1,1.4), c("#377EB8","#F0F0F0","#E41A1C")),
                cluster_rows = F,cluster_columns = F,
                show_column_names = F,show_row_names = F,
                column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),
                column_gap = unit(2, "mm"),
                row_split = forCNVanno$subgroup,
                
                heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
                
                top_annotation = top_anno,
                left_annotation = left_anno,
                bottom_annotation = ha,
                
                row_title = NULL,
                column_title = "patient ID",
                column_title_side = c("top"),
                column_title_gp = gpar(fontsize = 20),
)

draw(ht_tp, heatmap_legend_side = "right")
dev.off()
