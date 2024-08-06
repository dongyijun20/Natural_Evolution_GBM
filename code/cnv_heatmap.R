library(tidyverse)

cell_group=read.table("./infercnv/younger_malignant/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header = T,sep = "\t",stringsAsFactors = F)
cell_group=cell_group[!str_detect(cell_group$cell_group_name,"non-malignant"),]
cell_group$cell_group_name=str_replace(cell_group$cell_group_name,"malignant\\.malignant\\.","")

group_cellcount=as.data.frame(table(cell_group$cell_group_name))
colnames(group_cellcount)=c("cell_group_name","cellcount")
group_cellcount$cellratio=group_cellcount$cellcount / sum(group_cellcount$cellcount)

group_cnvtype=read.table("./infercnv/younger_malignant/CNVgroup_and_CNVtype_in_sampleA.txt",header = T,sep = "\t",stringsAsFactors = F)
group_cnvtype=group_cnvtype%>%inner_join(group_cellcount,by="cell_group_name")

alltype=c()
for (i in 1:22) {
  for (j in c("p","q")) {
    for (k in c("gain","loss")) {
      alltype=append(alltype,paste("chr",i,j,"_",k,sep = ""))
    }
  }
}

cellpercent=c()
for (i in alltype) {
  if(i %in% unique(group_cnvtype$cnv_type)){
    tmp=sum(group_cnvtype[group_cnvtype$cnv_type == i,"cellratio"])
    cellpercent=append(cellpercent,tmp)
  }else{
    cellpercent=append(cellpercent,0)
  }
}
names(cellpercent)=alltype
one.sample.stat=as.data.frame(cellpercent)
head(one.sample.stat)
# cellpercent
# chr1p_gain   0.0000000
# chr1p_loss   1.0000000
# chr1q_gain   0.0000000
# chr1q_loss   0.1197183
# chr2p_gain   0.0000000
# chr2p_loss   0.0000000

#### 当有多个样本时，直接合并就可以了。这里我为了演示，人为多加了几列 ####
some.sample.stat=one.sample.stat
colnames(some.sample.stat)="sampleA"
# some.sample.stat$sampleB=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleC=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleD=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleE=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleF=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleG=sample(some.sample.stat$sampleA,88,replace = T)
# some.sample.stat$sampleH=sample(some.sample.stat$sampleA,88,replace = T)

library(RColorBrewer)
library(scales)
library(pheatmap)

pheatmap(some.sample.stat,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         #rev(brewer.pal(n = 7, name ="RdYlBu"))
         #brewer.pal(9,"PuRd")
         #brewer.pal(9,"RdPu"),
         border_color = "grey",
         filename = "plot_new/malignant_younger_CNV.heatmap.pdf",width = 4,height = 11
)

