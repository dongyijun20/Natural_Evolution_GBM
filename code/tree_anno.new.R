library(tidyverse)
chr_pq=read.table("reference/chr_pq.new.txt",header = T,sep = "\t",stringsAsFactors = F) 
chr_pq$chr=paste("chr",chr_pq$chr,sep = "")

cnv_regions=read.table("./infercnv/older_malignant/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_regions.dat",header = T,sep = "\t",stringsAsFactors = F)
cnv_regions=cnv_regions%>%filter(state!=3)
cnv_regions$cell_group_name=str_replace(cnv_regions$cell_group_name,"malignant\\.malignant\\.","")
cnv_regions$cnv_type=""

for (i in 1:nrow(cnv_regions)) {
  tmp_chr_pq=chr_pq%>%filter(chr==cnv_regions[i,"chr"])
  if(cnv_regions[i,"start"] <= tmp_chr_pq[1,"cutoff"] & cnv_regions[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
    cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
  }else if (cnv_regions[i,"start"] >= tmp_chr_pq[2,"cutoff"] & cnv_regions[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
    cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
  }else {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"chr"],"p,",cnv_regions[i,"chr"],"q",sep = "")
  }
  if (cnv_regions[i,"state"] < 3) {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_loss",sep = "")
  }
  if (cnv_regions[i,"state"] > 3) {
    cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_gain",sep = "")
  }
  if (str_detect(cnv_regions[i,"cnv_type"],",")) {
    tmp1=strsplit(cnv_regions[i,"cnv_type"],",")[[1]][1]
    tmp2=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
    tmp3=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
    cnv_regions[i,"cnv_type"]=paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
    rm(list = c("tmp1","tmp2","tmp3"))
  }
}

#到这儿，我们得到的是一个比较完整的数据框，里面存储的是region对应的染色体位置，即位于长臂或者短臂。
#后续给肿瘤进化树做注释就是基于这些信息。例如, chr1p_loss说明的是在chr1p这条臂上面，存在loss的区域(区域大小未知)。
#我们能不能将chr1p_loss所表示的信息 更进一层，比如当出现chr1p_loss，我们就能知道chr1的短臂上有CNV，而且还是arm level的(就算不是arm level的，CNV区域也比较大)
#这两种解释都可以，本次代码更新针对第二种解释。
#这种解释之下，如果某个region的范围很小，仅因为这个region就把chr1p_loss标上去，就有点误导了。
#所以接下来，会对region的长度做些限制。
##########################################################
#我这里以染色体臂长度的30%作为阈值，臂内CNV的region之和超过这个长度才会保留下来
cnv_regions_copy=cnv_regions

cnv_regions_part1=cnv_regions[str_detect(cnv_regions$cnv_type,","),]
cnv_regions_part2=cnv_regions[!str_detect(cnv_regions$cnv_type,","),]
cnv_regions_part1_new=as.data.frame(matrix(NA,ncol = ncol(cnv_regions_part1), nrow = nrow(cnv_regions_part1)*2 ))
colnames(cnv_regions_part1_new)=colnames(cnv_regions_part1)
for (i in 1:dim(cnv_regions_part1)[1]) {
  cnv_regions_part1_new[2*i-1,]=cnv_regions_part1[i,]
  cnv_regions_part1_new[2*i,]  =cnv_regions_part1[i,]
  
  tmp_chr_pq=chr_pq[chr_pq$chr == cnv_regions_part1[i,"chr"],]
  cnv_regions_part1_new[2*i-1,"end"] = tmp_chr_pq[1,"cutoff"]
  cnv_regions_part1_new[2*i,"start"] = tmp_chr_pq[2,"cutoff"]
  
  cnv_regions_part1_new[2*i-1,"cnv_type"] = strsplit(cnv_regions_part1_new[2*i-1,"cnv_type"],",")[[1]][1]
  cnv_regions_part1_new[2*i,"cnv_type"]   = strsplit(cnv_regions_part1_new[2*i,  "cnv_type"],",")[[1]][2]
}
cnv_regions=cnv_regions_part2 %>% rbind(cnv_regions_part1_new)
cnv_regions$chr=factor(cnv_regions$chr,levels = paste("chr",1:22,sep = ""))
cnv_regions=cnv_regions%>%arrange(cell_group_name,chr,start)
cnv_regions$region_len=cnv_regions$end-cnv_regions$start+1
cnv_regions=cnv_regions[,c("cell_group_name","cnv_type","region_len")]
cnv_regions=as.data.frame(cnv_regions%>%group_by(cell_group_name,cnv_type)%>%dplyr::summarize(all_region_len=sum(region_len)))

####
cnv_regions$final_label=""
ratio=0.93

for (i in 1:dim(cnv_regions)[1]) {
  arm=strsplit(cnv_regions[i,"cnv_type"],"_")[[1]][1]
  ref_len=chr_pq[paste(chr_pq$chr,chr_pq$arm,sep = "") == arm,"arm_len"]
  
  if(cnv_regions[i,"all_region_len"] >= ref_len*ratio){
    cnv_regions[i,"final_label"]=cnv_regions[i,"cnv_type"]
  }else{
    cnv_regions[i,"final_label"]=""
  }
}

cnv_regions$cnv_type=NULL
colnames(cnv_regions)[3]="cnv_type"
cnv_regions=cnv_regions%>%filter(cnv_type != "")

write.table(cnv_regions[,c("cell_group_name","cnv_type")],file = "infercnv/older_malignant/CNVgroup_and_CNVtype_in_sampleA.txt",
            quote = F,sep = "\t",row.names = F,col.names = T)
#这个文件可以用来画热图，见另一个脚本
##########################################################

#CNV 有8类的代码类似
cnv_1.1.1.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.1"))[,"cnv_type"]
cnv_1.1.1.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.2"))[,"cnv_type"]
cnv_1.1.2.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.2.2"))[,"cnv_type"]
cnv_1.1.2.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.2.1"))[,"cnv_type"]
cnv_1.2.1.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.1"))[,"cnv_type"]
cnv_1.2.1.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.2"))[,"cnv_type"]
cnv_1.2.2.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.1"))[,"cnv_type"]
cnv_1.2.2.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.2"))[,"cnv_type"]

###cnv_1.1.1
cnv_1.1.1=intersect(cnv_1.1.1.1,cnv_1.1.1.2)
cnv_1.1.1.1_uniq=setdiff(cnv_1.1.1.1,cnv_1.1.1)
cnv_1.1.1.2_uniq=setdiff(cnv_1.1.1.2,cnv_1.1.1)
###cnv_1.1.2
cnv_1.1.2=intersect(cnv_1.1.2.2,cnv_1.1.2.1)
cnv_1.1.2.2_uniq=setdiff(cnv_1.1.2.2,cnv_1.1.2)
cnv_1.1.2.1_uniq=setdiff(cnv_1.1.2.1,cnv_1.1.2)
###cnv_1.2.1
cnv_1.2.1=intersect(cnv_1.2.1.1,cnv_1.2.1.2)
cnv_1.2.1.1_uniq=setdiff(cnv_1.2.1.1,cnv_1.2.1)
cnv_1.2.1.2_uniq=setdiff(cnv_1.2.1.2,cnv_1.2.1)
###cnv_1.2.2
cnv_1.2.2=intersect(cnv_1.2.2.1,cnv_1.2.2.2)
cnv_1.2.2.1_uniq=setdiff(cnv_1.2.2.1,cnv_1.2.2)
cnv_1.2.2.2_uniq=setdiff(cnv_1.2.2.2,cnv_1.2.2)
###cnv_1.1
cnv_1.1=intersect(cnv_1.1.1,cnv_1.1.2)
cnv_1.1.1_uniq=setdiff(cnv_1.1.1,cnv_1.1)
cnv_1.1.2_uniq=setdiff(cnv_1.1.2,cnv_1.1)
###cnv_1.2
cnv_1.2=intersect(cnv_1.2.1,cnv_1.2.2)
cnv_1.2.1_uniq=setdiff(cnv_1.2.1,cnv_1.2)
cnv_1.2.2_uniq=setdiff(cnv_1.2.2,cnv_1.2)
###cnv_1
cnv_1=intersect(cnv_1.1,cnv_1.2)
cnv_1.1_uniq=setdiff(cnv_1.1,cnv_1)
cnv_1.2_uniq=setdiff(cnv_1.2,cnv_1)

cat("cnv_1:",cnv_1,"\n",file="out.txt",append = T)
cat("cnv_1.1_uniq:",cnv_1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2_uniq:",cnv_1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1_uniq:",cnv_1.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.2_uniq:",cnv_1.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1_uniq:",cnv_1.2.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2_uniq:",cnv_1.2.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1.1_uniq:",cnv_1.1.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.1.1.2_uniq:",cnv_1.1.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1.1_uniq:",cnv_1.2.1.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.1.2_uniq:",cnv_1.2.1.2_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2.1_uniq:",cnv_1.2.2.1_uniq,"\n",file="out.txt",append = T)
cat("cnv_1.2.2.2_uniq:",cnv_1.2.2.2_uniq,"\n",file="out.txt",append = T)

